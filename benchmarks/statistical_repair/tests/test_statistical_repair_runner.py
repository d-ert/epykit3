"""Benchmark resumability must not mix sources or trust altered outputs."""

import json
from pathlib import Path

import pytest


def test_resume_rejects_changed_source_identity(tmp_path):
    from scripts.run_statistical_repair_benchmarks import completed_run

    (tmp_path / "run.json").write_text(json.dumps({"status": "complete", "identity": {"source": "old"}, "outputs": {}}))
    with pytest.raises(ValueError, match="identity"):
        completed_run(tmp_path, {"source": "new"})


def test_resume_checks_output_hashes(tmp_path):
    from scripts.run_statistical_repair_benchmarks import completed_run, sha256

    output = tmp_path / "candidates.tsv"
    output.write_text("original")
    identity = {"source": "same", "coverage": 10}
    (tmp_path / "run.json").write_text(json.dumps({"status": "complete", "identity": identity,
                                                "outputs": {output.name: sha256(output)}}))
    assert completed_run(tmp_path, identity)
    output.write_text("altered")
    with pytest.raises(ValueError, match="hash"):
        completed_run(tmp_path, identity)


def test_incomplete_runs_are_not_reported_as_complete(tmp_path):
    from scripts.run_statistical_repair_benchmarks import completed_run

    identity = {"source": "same"}
    (tmp_path / "run.json").write_text(json.dumps({"status": "running", "identity": identity}))
    assert not completed_run(tmp_path, identity)


def test_aggregate_keeps_null_error_separate_from_signal_recovery():
    from scripts.run_statistical_repair_benchmarks import aggregate_metrics

    rows = [
        dict(kind="pilot", config_id="a", dataset="empirical_signal_seed1", is_null=False,
             recall_cpg50=.5, any_signal_overlap_fdp=.1, n_calls=10, elapsed_seconds=2),
        dict(kind="pilot", config_id="a", dataset="empirical_signal_seed2", is_null=False,
             recall_cpg50=.25, any_signal_overlap_fdp=.2, n_calls=8, elapsed_seconds=3),
        dict(kind="pilot", config_id="a", dataset="empirical_null_seed1", is_null=True,
             recall_cpg50=0, any_signal_overlap_fdp=1, n_calls=2, elapsed_seconds=1),
    ]
    result = aggregate_metrics(rows)
    assert result["a"]["empirical_main_recall"] == .375
    assert result["a"]["empirical_main_signal_overlap_fdp"] == pytest.approx(.15)
    assert result["a"]["null_discovery_counts"] == [2]


@pytest.fixture
def provenance_fixture(tmp_path, monkeypatch):
    import subprocess
    from scripts import run_statistical_repair_benchmarks as runner

    repo = tmp_path / "repo"
    package = repo / "src/epykit"
    package.mkdir(parents=True)
    (package / "__init__.py").write_text("VERSION=1\n")
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    subprocess.run(["git", "-C", str(repo), "add", "."], check=True)
    subprocess.run(["git", "-C", str(repo), "-c", "user.name=Fixture", "-c", "user.email=fixture@example.test",
                    "commit", "-qm", "fixture"], check=True)
    root = tmp_path / "workspace"
    for name in ["scripts/benchmark.py", "scripts/validation.py", "scripts/helper.py", "wgbs_sim/evaluate.py"]:
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("original\n")
    monkeypatch.setattr(runner, "ROOT", root)
    return runner, repo, root, tmp_path / "out"


@pytest.mark.parametrize("changed", ["scripts/benchmark.py", "scripts/helper.py", "wgbs_sim/evaluate.py"])
def test_resume_rejects_changed_evaluation_dependency(provenance_fixture, changed):
    runner, repo, root, output = provenance_fixture
    runner.snapshot_source(repo, output, entrypoint=root / "scripts/benchmark.py")
    (root / changed).write_text("changed\n")
    with pytest.raises(ValueError, match="evaluation"):
        runner.snapshot_source(repo, output, entrypoint=root / "scripts/benchmark.py")


def test_validation_provenance_records_actual_entrypoint_and_checks_frozen_dependencies(provenance_fixture):
    runner, repo, root, output = provenance_fixture
    entrypoint = root / "scripts/validation.py"
    info = runner.snapshot_source(repo, output, entrypoint=entrypoint)
    assert info["runner_sha256"] == runner.sha256(entrypoint)
    assert info["entrypoint"] == "scripts/validation.py"
    (output / "evaluation_source/scripts/helper.py").write_text("altered frozen helper\n")
    with pytest.raises(ValueError, match="evaluation"):
        runner.snapshot_source(repo, output, entrypoint=entrypoint)


def test_pilot_evaluation_identity_includes_both_truth_files(tmp_path):
    from scripts.run_statistical_repair_benchmarks import pilot_truth_identity

    for name in ["truth.cpg.tsv.gz", "truth.regions.tsv"]:
        (tmp_path / name).write_bytes(b"original")
    before = pilot_truth_identity(tmp_path)
    (tmp_path / "truth.regions.tsv").write_bytes(b"changed")
    after = pilot_truth_identity(tmp_path)
    assert before != after
    assert before["truth.cpg.tsv.gz"] == after["truth.cpg.tsv.gz"]
    assert before["truth.regions.tsv"] != after["truth.regions.tsv"]


def test_dmc_scoring_preserves_zero_to_one_based_coordinates(tmp_path):
    import gzip
    pl = pytest.importorskip("polars", reason="Epykit benchmark scoring requires its optional environment")
    from scripts.run_statistical_repair_benchmarks import dmc_metrics

    with gzip.open(tmp_path / 'truth.cpg.tsv.gz', 'wt') as f:
        f.write('chrom\tpos\tis_nonzero\tis_relevant\nchr1\t1\t1\t1\nchr1\t21\t0\t0\nchr1\t41\t1\t1\n')
    table = pl.DataFrame(dict(chrom=['chr1']*3, pos=[0,20,40], qvalue=[.001,.01,.8]))
    score = dmc_metrics(table, tmp_path)
    assert score['n_dmc_calls'] == 2
    assert score['dmc_nonzero_recall'] == .5
    assert score['dmc_nonzero_fdp'] == .5
