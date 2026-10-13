import importlib.util
import subprocess
import sys
from pathlib import Path

import pandas as pd
import polars as pl

import epykit as ep


def make_store(root):
    rows = []
    for j in range(6):
        sample = f"s{j}"
        rows.append(dict(sample_id=sample, group="case" if j < 3 else "control", treatment=int(j < 3)))
        positions = [0, 30, 70] if j != 5 else [0, 70]
        directory = root / f"sample={sample}" / "chrom=chr1"
        directory.mkdir(parents=True)
        pl.DataFrame(dict(pos=positions, strand=["*"] * len(positions),
                          N_meth=[12 if j < 3 else 5] * len(positions),
                          N_unmeth=[8 if j < 3 else 15] * len(positions),
                          coverage=[20] * len(positions))).write_parquet(directory / "part-0.parquet")
    return pl.DataFrame(rows)


def test_effective_counts_describe_observed_samples_not_group_size(tmp_path):
    obs = make_store(tmp_path / "counts")
    md = ep.MethylData(obs=obs, store=str(tmp_path / "counts"), analysis_root=str(tmp_path / "analysis"))
    ep.pp.set_unite_type(md, type="union")
    ep.tl.dmc(md, test="lr", min_samples_treatment=2, min_samples_control=2)
    result = md.dmc.sort("pos")
    assert result["n_case"].to_list() == [3, 3, 3]
    assert result["n_control"].to_list() == [3, 3, 3]
    assert "n_valid_case" in result.columns
    assert result["n_valid_case"].to_list() == [3, 3, 3]
    assert result["n_valid_control"].to_list() == [3, 2, 3]


def test_external_count_store_never_receives_output_cache(tmp_path):
    obs = make_store(tmp_path / "counts")
    md = ep.MethylData(obs=obs, store=str(tmp_path / "counts"), analysis_root=str(tmp_path / "analysis"))
    ep.pp.set_unite_type(md, type="intersect")
    ep.tl.dmc(md, test="lr", materialize=False)
    assert Path(md.uns["dmc"]["store_path"]).is_relative_to(tmp_path / "analysis")
    assert not (tmp_path / ".cache/dmc").exists()


def test_workflow_preserves_intersection_and_fresh_output(tmp_path):
    make_store(tmp_path / "counts")
    sheet = tmp_path / "samples.csv"
    pd.DataFrame([dict(sample_id=f"s{j}", group="case" if j < 3 else "control", path="unused")
                  for j in range(6)]).to_csv(sheet, index=False)
    script = Path(__file__).parents[1] / "run_power_baseline.py"
    output = tmp_path / "run"
    args = [sys.executable, str(script), "--samplesheet", str(sheet), "--treatment-group", "case",
            "--control-group", "control", "--input-store", str(tmp_path / "counts"),
            "--output-dir", str(output), "--no-report"]
    run = subprocess.run(args, capture_output=True, text=True, cwd=tmp_path)
    assert run.returncode == 0, run.stderr
    result = pl.read_csv(output / "tables/dmc_lr.tsv", separator="\t")
    assert result["pos"].to_list() == [0, 70]
    assert result["n_valid_control"].to_list() == [3, 3]
    second = subprocess.run(args, capture_output=True, text=True, cwd=tmp_path)
    assert second.returncode != 0
    assert "already exists" in second.stderr


def test_paper_coordinates_do_not_create_boundary_overlap():
    path = Path(__file__).parents[1] / "compare_dmrs.py"
    spec = importlib.util.spec_from_file_location("comparison", path)
    comparison = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(comparison)
    paper = pd.DataFrame(dict(chr=["chr1"], start=[11], end=[20], **{"diff.meth_mean": [0.2]}))
    calls = pd.DataFrame(dict(chrom=["chr1"], start=[20], end=[30], mean_meth_diff=[0.2]))
    p, e = comparison.standardize_columns(paper, calls)
    assert comparison.find_overlaps(p, e, min_overlap_bp=1).empty


def test_smyth_eb_dispersion_moderation():
    import numpy as np
    from epykit.dmc import _smyth_eb_moderation

    rng = np.random.default_rng(42)
    n = 1000
    df = 4.0
    true_phi = 1.8
    # Raw Pearson statistics under quasi-binomial with df=4
    chi_sum = true_phi * rng.chisquare(df, size=n)
    df_arr = np.full(n, df)
    valid = np.ones(n, dtype=bool)

    phi_eff, df_phi, d0, phi0 = _smyth_eb_moderation(chi_sum, df_arr, valid, min_dispersion=1.0)
    assert d0 >= 20.0
    assert 1.0 <= phi0 <= 3.0
    assert np.all(df_phi >= df + 20.0)
    assert np.all(phi_eff >= 1.0)

    # Fallback with few sites
    phi_eff_few, df_phi_few, d0_few, phi0_few = _smyth_eb_moderation(
        chi_sum[:20], df_arr[:20], valid[:20], min_dispersion=1.0, min_disp_sites=100
    )
    assert d0_few == 0.0
    assert np.all(df_phi_few == df)


def test_spatial_stouffer_variance_deflation_and_length_scale():
    import numpy as np
    from epykit.dmr import (
        _estimate_chrom_z_length_scale,
        _spatial_variance_sum,
        _stouffer_combine_signed,
        call_dmr_chain_merge,
        DMR_PRESETS,
    )

    # 1. Test length scale estimation on simulated AR(1) spatial Z-scores
    rng = np.random.default_rng(123)
    n = 1000
    positions = np.sort(rng.choice(100_000, size=n, replace=False))
    # AR(1) with decay lambda = 100 bp
    z = np.zeros(n)
    for i in range(1, n):
        d = positions[i] - positions[i - 1]
        rho = np.exp(-d / 100.0)
        z[i] = rho * z[i - 1] + np.sqrt(max(1.0 - rho * rho, 0.0)) * rng.standard_normal()

    lam_est = _estimate_chrom_z_length_scale(positions, z)
    assert 40.0 <= lam_est <= 250.0

    # 2. Test spatial variance deflation on dense cluster
    dense_pos = np.arange(8, dtype=float) * 20.0  # 8 CpGs separated by 20 bp
    w = np.ones(8)
    var_indep = 8.0
    var_spatial = _spatial_variance_sum(dense_pos, w, length_scale=100.0)
    # Variance should be significantly inflated due to pairwise positive correlations
    assert var_spatial > 2.0 * var_indep

    # Combined p-value should be properly deflated (larger p-value under spatial correction)
    pvals = np.full(8, 0.01)
    diffs = np.full(8, 0.2)
    p_uncorr = _stouffer_combine_signed(pvals, diffs)
    p_corr = _stouffer_combine_signed(pvals, diffs, positions=dense_pos, length_scale=100.0)
    assert p_corr > p_uncorr
    assert p_corr > 1e-8  # uncorrected gives ~3e-13

    # 3. Test DMR_PRESETS["default"] configuration
    assert DMR_PRESETS["default"]["alpha"] == 0.005
    assert DMR_PRESETS["default"]["min_cpgs"] == 5

    # 4. End-to-end chain merge with preset="default"
    dmc_df = pl.DataFrame({
        "chrom": ["chr1"] * 10,
        "pos": [100, 120, 140, 160, 180, 200, 220, 240, 260, 280],
        "meth_diff": [0.25] * 10,
        "pvalue": [0.001] * 10,
    })
    dmrs = call_dmr_chain_merge(dmc_df, preset="default")
    assert len(dmrs) == 1
    assert dmrs["n_cpgs"][0] == 10
    assert dmrs["combined_pvalue"][0] > 0.0

