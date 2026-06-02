"""Layer 6a: CLI smoke + end-to-end pipeline tests.

These run the ``epykit`` console script as a real subprocess so we catch:

* argparse wiring (flag names, defaults, choices)
* the module-level no-side-effects rule on import
* exit codes
* the n=1 guard at the CLI level (B6)
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest


def _epykit(*args, check: bool = True, **kwargs) -> subprocess.CompletedProcess:
    """Invoke epykit as a module so the test works without a console
    script symlink installed on PATH (``pip install -e .`` may be missing
    in dev shells, but ``python -m epykit.cli`` always works)."""
    cmd = [sys.executable, "-m", "epykit.cli", *map(str, args)]
    return subprocess.run(
        cmd, capture_output=True, text=True, check=check, **kwargs,
    )



# Help text + defaults


def test_cli_top_level_help():
    """Top-level ``--help`` lists every subcommand."""
    res = _epykit("--help", check=False)
    # argparse returns 0 on --help.
    assert res.returncode == 0, res.stderr
    for sub in ("convert", "filter", "dmc", "dmr", "annotate", "qc-report", "smooth"):
        assert sub in res.stdout, f"subcommand {sub!r} missing from --help"


def test_cli_dmc_help_shows_lr_default_and_allow_n1():
    """B2 + B6: the dmc help must document --test default=lr and --allow-n1."""
    res = _epykit("dmc", "--help", check=False)
    assert res.returncode == 0, res.stderr
    out = res.stdout.lower()
    assert "--allow-n1" in res.stdout, "--allow-n1 missing from dmc help"
    # lr should be mentioned as the (default) choice; argparse prints
    # `default: lr` or includes it in the choices line.
    assert "lr" in out, "lr should appear in dmc --test help"


def test_cli_dmr_help_shows_method_and_allow_n1():
    res = _epykit("dmr", "--help", check=False)
    assert res.returncode == 0, res.stderr
    assert "--method" in res.stdout
    assert "tile" in res.stdout
    assert "sliding_window" in res.stdout
    assert "--allow-n1" in res.stdout
    # "hmm" must NOT appear in the choices list after the 1.0 removal.
    assert "hmm" not in res.stdout, "'hmm' still appears in dmr --help choices"


def test_cli_dmr_method_hmm_rejected():
    """CLI rejects --method=hmm now that the deprecation shim is gone.

    argparse enforces the choices list and exits non-zero with an 'invalid
    choice' message when an unlisted value is passed.  The library-level
    tl.dmr(method='hmm') raises ValueError as of 1.0; this test confirms the
    CLI surface is consistent (no shim translating hmm -> segment at the
    argparse layer).
    """
    res = _epykit("dmr", "--method", "hmm", "--output", "out.parquet", check=False)
    assert res.returncode != 0, (
        "CLI should reject --method=hmm with a non-zero exit code"
    )
    # argparse prints 'invalid choice' on stderr.
    assert "invalid choice" in res.stderr.lower() or "hmm" in res.stderr, (
        f"Expected an argparse rejection message mentioning 'hmm'; got: {res.stderr[:300]!r}"
    )


def test_cli_dmr_help_shows_chain_merge_default():
    """chain_merge is the tl.dmr library default; the CLI must expose it as a
    --method choice, surface --preset, and document it as the default."""
    res = _epykit("dmr", "--help", check=False)
    assert res.returncode == 0, res.stderr
    assert "chain_merge" in res.stdout, "chain_merge missing from dmr --method choices"
    assert "--preset" in res.stdout, "--preset missing from dmr help"
    assert "default" in res.stdout.lower()


def test_cli_dmr_chain_merge_requires_dmc_results():
    """chain_merge without --dmc-results must fail with a non-zero exit."""
    res = _epykit(
        "dmr", "--method", "chain_merge", "--output", "out.parquet", check=False
    )
    assert res.returncode != 0, "chain_merge without --dmc-results should fail"
    assert "dmc-results" in (res.stderr + res.stdout).lower()


def test_cli_dmr_chain_merge_end_to_end(tmp_path):
    """chain_merge CLI path turns a synthetic DMC parquet into DMRs.

    Fast: builds a tiny in-memory DMC table (the only columns chain_merge
    reads are chrom/pos/meth_diff/pvalue) rather than walking a methylstore.
    """
    import polars as pl

    dmc = pl.DataFrame(
        {
            "chrom": ["chr1"] * 5,
            "pos": [100, 150, 200, 250, 300],
            "meth_diff": [0.3, 0.3, 0.3, 0.3, 0.3],
            "pvalue": [1e-4] * 5,
        }
    )
    dmc_path = tmp_path / "dmc.parquet"
    dmc.write_parquet(str(dmc_path))

    out = tmp_path / "dmr.parquet"
    res = _epykit(
        "dmr", "--method", "chain_merge",
        "--dmc-results", str(dmc_path),
        "--output", str(out),
        check=False,
    )
    assert res.returncode == 0, f"stdout={res.stdout}\nstderr={res.stderr}"
    assert out.exists()

    df = pl.read_parquet(str(out))
    for col in ("chrom", "start", "end", "dmr_type"):
        assert col in df.columns, f"DMR output missing {col!r} column"
    assert len(df) >= 1, "chain_merge should call at least one DMR on this cluster"


def test_cli_smooth_help_says_gaussian_not_bsmooth():
    """B5: the smooth subcommand help should no longer mislead users into
    thinking they're getting BSmooth LOESS."""
    res = _epykit("smooth", "--help", check=False)
    assert res.returncode == 0, res.stderr
    assert "Gaussian" in res.stdout or "gaussian" in res.stdout



# Import-time side effects


def test_importing_epykit_does_not_clobber_root_logging(tmp_path):
    """Importing the package after configuring logging must NOT change
    the root logger's handlers or formatters (S2: ``logging.basicConfig``
    moved into ``main()``)."""
    script = tmp_path / "probe.py"
    script.write_text(
        "import logging\n"
        "import sys\n"
        "logging.basicConfig(level=logging.DEBUG, format='SENTINEL %(message)s')\n"
        "import epykit  # noqa\n"
        "h = logging.getLogger().handlers[0]\n"
        "sys.stdout.write(h.formatter._fmt + '\\n')\n"
    )
    res = subprocess.run(
        [sys.executable, str(script)], capture_output=True, text=True, check=False
    )
    assert "SENTINEL %(message)s" in res.stdout, (
        f"import epykit clobbered root logging config\n"
        f"stdout: {res.stdout!r}\nstderr: {res.stderr!r}"
    )



# End-to-end pipeline through the CLI (slow)


@pytest.mark.slow
def test_cli_end_to_end_dmc_pipeline(tmp_path, synth_bundle):
    """Full pipeline via the CLI: convert is implicit through epykit's
    read paths, but ``epykit dmc`` runs DMC against the methylstore that
    the session fixture has already filled.

    This is slow (~30 s) because it walks every chromosome -- gated behind
    the ``slow`` marker.
    """
    # Build a fresh methylstore via the Python API; the CLI dmc command
    # operates on the filtered store.
    import epykit as ep
    md = ep.read_bismark(
        synth_bundle.samplesheet,
        treatment_group="treatment",
        control_group="control",
        store_dir=str(tmp_path / "store"),
    )
    ep.pp.filter_coverage(md, lo_count=5, hi_perc=99.9)
    filtered_store = md.store

    output = tmp_path / "dmc.parquet"
    res = _epykit(
        "dmc",
        "--methylstore", filtered_store,
        "--samplesheet", synth_bundle.samplesheet,
        "--treatment-group", "treatment",
        "--control-group",   "control",
        "--output",          str(output),
        "--test",            "lr",
        check=False,
    )
    assert res.returncode == 0, (
        f"epykit dmc failed\nstdout: {res.stdout}\nstderr: {res.stderr}"
    )
    assert output.exists()

    # Sanity check: the output Parquet contains the expected schema.
    import polars as pl
    df = pl.read_parquet(str(output))
    assert "pvalue" in df.columns
    assert "qvalue" in df.columns
    assert len(df) > 0


@pytest.mark.slow
def test_cli_dmc_refuses_n1_without_allow_n1_flag(tmp_path):
    """The CLI mirror of the Python guard: n=1 without --allow-n1 must
    exit non-zero with a helpful error."""
    import epykit as ep
    from tests.fixtures.synth import SimConfig, generate

    cfg = SimConfig(
        n_per_group=1,
        chromosomes=("chr1",),
        cpgs_per_chrom=100,
        n_scattered_dmcs=5,
        n_dmrs=1,
        dmr_size_cpgs=3,
        seed=11,
    )
    result = generate(cfg, tmp_path / "n1")
    md = ep.read_bismark(
        result["samplesheet"],
        treatment_group="treatment",
        control_group="control",
        store_dir=str(tmp_path / "store"),
    )
    ep.pp.filter_coverage(md, lo_count=2, hi_perc=99.9)
    store = md.store

    res = _epykit(
        "dmc",
        "--methylstore", store,
        "--samplesheet", result["samplesheet"],
        "--treatment-group", "treatment",
        "--control-group",   "control",
        "--output",          str(tmp_path / "dmc_n1.parquet"),
        "--test", "lr",
        check=False,
    )
    assert res.returncode != 0, "CLI should refuse n=1 without --allow-n1"
    assert "at least 2 replicates" in (res.stderr + res.stdout)
