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
