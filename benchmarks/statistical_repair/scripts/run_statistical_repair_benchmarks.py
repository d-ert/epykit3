#!/usr/bin/env python3
"""Resumable, source-frozen GLS/DMC repair benchmarks with preserved baselines."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import math
import os
import shutil
import subprocess
import sys
import time
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPAIR_REPO = ROOT / "development/epykit_statistical_repair"
CONFIGS = {
    "cov10_boundNone": dict(coverage_min=10, boundary_min_effect=None),
    "cov10_bound05": dict(coverage_min=10, boundary_min_effect=.05),
    "cov5_boundNone": dict(coverage_min=5, boundary_min_effect=None),
    "cov5_bound05": dict(coverage_min=5, boundary_min_effect=.05),
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _safe_json(value):
    if isinstance(value, dict):
        return {str(k): _safe_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_safe_json(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if hasattr(value, "item"):
        return _safe_json(value.item())
    return value


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(_safe_json(value), indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def completed_run(directory, identity):
    manifest = Path(directory) / "run.json"
    if not manifest.is_file():
        return False
    run = json.loads(manifest.read_text())
    if run.get("identity") != identity:
        raise ValueError(f"Benchmark identity changed: {directory}")
    if run.get("status") != "complete":
        return False
    for name, expected in run.get("outputs", {}).items():
        path = Path(directory) / name
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f"Completed benchmark output hash changed: {path}")
    return True


def aggregate_metrics(records):
    result = {}
    for config in sorted({r["config_id"] for r in records}):
        selected = [r for r in records if r["config_id"] == config and r["kind"] == "pilot"]
        summary = {}
        for track in ("empirical", "generative"):
            main = [r for r in selected if r["dataset"].startswith(f"{track}_signal_seed")]
            for label, key in [("recall", "recall_cpg50"), ("signal_overlap_fdp", "any_signal_overlap_fdp")]:
                values = [r[key] for r in main if r.get(key) is not None]
                summary[f"{track}_main_{label}"] = sum(values)/len(values) if values else None
        null = [r for r in selected if r["is_null"]]
        summary["null_discovery_counts"] = [r["n_calls"] for r in null]
        summary["null_datasets"] = [r["dataset"] for r in null]
        summary["any_null_discovery_fraction"] = sum(r["n_calls"] > 0 for r in null)/len(null) if null else None
        result[config] = summary
    return result


def _evaluation_identity(entrypoint):
    files = set(ROOT.glob("scripts/*.py")) | set(ROOT.glob("wgbs_sim/**/*.py"))
    files.add(entrypoint)
    hashes = {str(p.relative_to(ROOT)): sha256(p) for p in sorted(files)}
    return hashes, hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest()


def verify_evaluation(info, output=None):
    """Prevent code changes during a run from producing a mixed completed set."""
    for name, expected in info["evaluation_sha256_by_file"].items():
        if sha256(ROOT / name) != expected:
            raise ValueError(f"Live evaluation dependency changed: {name}")
        if output is not None and sha256(output / "evaluation_source" / name) != expected:
            raise ValueError(f"Frozen evaluation dependency changed: {name}")


def pilot_truth_identity(directory):
    return {name: sha256(directory / name) for name in ("truth.cpg.tsv.gz", "truth.regions.tsv")}


def snapshot_source(repo, output, *, entrypoint=None):
    entrypoint = Path(__file__).resolve() if entrypoint is None else Path(entrypoint).resolve()
    package = repo / "src/epykit"
    dirty = subprocess.check_output(["git", "-C", str(repo), "status", "--porcelain", "--", "src"], text=True)
    if dirty:
        raise ValueError("Commit source changes before benchmarking")
    hashes = {str(p.relative_to(package)): sha256(p) for p in sorted(package.rglob("*"))
              if p.is_file() and "__pycache__" not in p.parts and p.name != ".DS_Store"}
    signature = hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest()
    evaluation_hashes, evaluation_signature = _evaluation_identity(entrypoint)
    manifest = output / "source_provenance.json"
    if manifest.exists():
        info = json.loads(manifest.read_text())
        if info["source_signature"] != signature:
            raise ValueError("This benchmark directory belongs to another source revision")
        if (info.get("evaluation_signature") != evaluation_signature
                or info.get("entrypoint") != str(entrypoint.relative_to(ROOT))
                or info.get("runner_sha256") != sha256(entrypoint)
                or info.get("evaluation_sha256") != sha256(ROOT / "wgbs_sim/evaluate.py")):
            raise ValueError("Benchmark evaluation identity changed; use a new versioned output directory")
        for name, expected in hashes.items():
            if sha256(output / "source/epykit" / name) != expected:
                raise ValueError("Frozen source snapshot hash mismatch")
        verify_evaluation(info, output)
        return info
    output.mkdir(parents=True, exist_ok=True)
    shutil.copytree(package, output / "source/epykit", ignore=shutil.ignore_patterns("__pycache__", ".DS_Store"))
    for name in evaluation_hashes:
        target = output / "evaluation_source" / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / name, target)
    info = dict(repo=str(repo), commit=subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip(),
                source_signature=signature, source_sha256=hashes, started_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                entrypoint=str(entrypoint.relative_to(ROOT)), runner_sha256=sha256(entrypoint),
                evaluation_sha256=sha256(ROOT / "wgbs_sim/evaluate.py"), evaluation_signature=evaluation_signature,
                evaluation_sha256_by_file=evaluation_hashes,
                python=sys.version, inference="regional permutation FDP estimates remain experimental")
    write_json(manifest, info)
    return info


def store_identity(store, manifest):
    hashes = {str(p.relative_to(store)): sha256(p) for p in sorted(store.rglob("*.parquet"))}
    if not hashes:
        raise ValueError(f"Count store has no Parquet files: {store}")
    signature = hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest()
    write_json(manifest, dict(store=str(store), sha256=hashes, signature=signature))
    return signature


def _write_tsv(path, records):
    if not records:
        return
    fields = list(dict.fromkeys(k for row in records for k in row))
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(records)


def pilot_metrics(candidates, dataset, config, elapsed, params, out):
    from scripts.gls_sensitivity import score_selection, select_regions, truth_arrays

    truth, universe, signal = truth_arrays(ROOT / "data/pilot", dataset["dataset"])
    calls = select_regions(candidates.to_dicts(), "combined_qvalue", .05)
    scores = score_selection(calls, truth, universe, signal)
    row = dict(kind="pilot", dataset=dataset["dataset"], track=dataset["track"], is_null=dataset["all_null"] == "True",
               config_id=config, **CONFIGS[config], elapsed_seconds=elapsed, n_candidates=len(candidates),
               common_cpgs=sum(params["common_cpgs"].values()), union_cpgs=sum(params["union_cpgs"].values()),
               null_regions=params["null_regions"], n_calls=scores["n_calls"], n_matches=scores["n_matches"],
               n_true=scores["n_true"], recall_cpg50=scores["recall"], precision=scores["precision"],
               recovery_fdp=scores["recovery_fdp"], any_signal_overlap_fdp=scores["any_signal_overlap_fdp"],
               recall_bp50=scores["recall_bp"], boundary_error_bp=scores["boundary_error_bp"],
               fragmented_truth_fraction=scores["fragmented_truth_fraction"], merged_call_fraction=scores["merged_call_fraction"])
    if "maxT_pvalue" in candidates.columns:
        import polars as pl
        row["n_maxT_calls"] = len(candidates.filter(pl.col("maxT_pvalue") <= .05))
    for shape in ("sharp", "tapered", "irregular", "nearby", "mixed_sign"):
        row[f"{shape}_recall"] = scores["scenario_recovery"].get(shape, {}).get("recall", 0.)
    write_json(out / "scenario_recovery.json", scores["scenario_recovery"])
    return row


def dmc_metrics(dmc, truth_directory=None):
    import polars as pl
    from scripts.run_epykit_pilot import normalize_dmc, table
    from wgbs_sim.evaluate import evaluate_sites

    selected = dmc.filter(pl.col("qvalue").is_finite() & (pl.col("qvalue") < .05))
    result = dict(n_dmc_calls=len(selected), n_dmc_tested=int(dmc["qvalue"].is_finite().sum()))
    if truth_directory is not None:
        scores = evaluate_sites(normalize_dmc(selected.to_dicts()), table(truth_directory / "truth.cpg.tsv.gz"))
        result.update(dmc_nonzero_recall=scores["nonzero"]["recall"], dmc_nonzero_fdp=scores["nonzero"]["fdp"],
                      dmc_relevant_recall=scores["relevant"]["recall"], dmc_relevant_fdp=scores["relevant"]["fdp"])
    return result


def real_metrics(ep, obs, store, candidates, config, elapsed, params, out):
    import polars as pl
    from scripts.compare_annotated_gls import METHODS, associations, set_comparison
    from scripts.compare_real_dmrs import WORKBOOK, read_workbook
    from scripts.run_real_gls import INPUT
    from wgbs_sim.evaluate import _overlap_pairs, evaluate_regions

    paper, _ = read_workbook(WORKBOOK)
    selected = candidates.filter(pl.col("combined_qvalue") <= .05)
    raw = selected.to_dicts()
    pairs = set(_overlap_pairs(raw, paper))
    reciprocal = evaluate_regions(raw, paper, {}, .5, "bp")
    refgene = INPUT.parent / "refseq/refGene.txt.gz"
    islands = INPUT.parent / "hg38_cpg_islands.bed"
    annotated = ep.MethylData(obs=obs, store=str(store), assembly="GRCh38", analysis_root=str(out / "annotation"))
    annotated.uns["dmr"] = selected
    if len(selected):
        ep.tl.annotate(annotated, refgene=str(refgene), cpg_islands=str(islands), tsv=False)
    annotated.uns["dmr"].write_parquet(out / "gls_annotated.parquet")
    genes = annotated.uns["dmr"].to_dicts()
    reported = associations(paper, "reported")
    nearest = associations(genes, "nearest")
    shared = set(nearest) & set(reported)
    gene_rows = []
    for method in METHODS:
        mapping = associations(genes, method)
        shared_method = set(mapping) & set(reported)
        coordinate_shared = {gene for gene in shared_method
                             if any((ci, pi) in pairs for ci in mapping[gene] for pi in reported[gene])}
        gene_rows.append(dict(config_id=config, method=method, reference="paper_reported", **set_comparison(mapping, reported),
                              shared_genes_with_coordinate_overlap=len(coordinate_shared)))
    _write_tsv(out / "gene_comparison.tsv", gene_rows)
    _write_tsv(out / "overlap_pairs.tsv", [dict(call_index=a, paper_index=b) for a, b in sorted(pairs)])
    shared_coords = {gene for gene in shared if any((ci, pi) in pairs for ci in nearest[gene] for pi in reported[gene])}
    return dict(kind="real", dataset="GSE263850", config_id=config, **CONFIGS[config], elapsed_seconds=elapsed,
                common_cpgs=sum(params["common_cpgs"].values()), union_cpgs=sum(params["union_cpgs"].values()),
                n_candidates=len(candidates), n_dmr_calls=len(selected), gls_overlapping_paper=len({a for a, _ in pairs}),
                paper_overlapping_gls=len({b for _, b in pairs}), overlap_pairs=len(pairs),
                same_direction_overlap_pairs=sum(raw[a]["mean_meth_diff"]*paper[b]["paper_effect"] > 0 for a, b in pairs),
                one_to_one_50pct_bp_matches=len(reciprocal["matches"]), shared_nearest_genes=len(shared),
                paper_shared_gene_fraction=len(shared)/len(reported), shared_genes_with_coordinate_overlap=len(shared_coords),
                n_maxT_calls=len(candidates.filter(pl.col("maxT_pvalue") <= .05)) if "maxT_pvalue" in candidates.columns else None)


def summarize(output):
    records = []
    for manifest in sorted(output.glob("runs/*/*/*/run.json")):
        run = json.loads(manifest.read_text())
        if run.get("status") == "complete":
            completed_run(manifest.parent, run["identity"])
            records.append(json.loads((manifest.parent / "metrics.json").read_text()))
    _write_tsv(output / "per_run.tsv", records)
    write_json(output / "metrics.json", records)
    write_json(output / "aggregate.json", aggregate_metrics(records))
    return records


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repo", type=Path, default=REPAIR_REPO)
    parser.add_argument("--kind", choices=["pilot", "real", "both"], default="both")
    parser.add_argument("--datasets", nargs="*")
    parser.add_argument("--configs", nargs="+", choices=list(CONFIGS), default=list(CONFIGS))
    parser.add_argument("--chromosomes", nargs="+")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--summarize-only", action="store_true")
    parser.add_argument("--include-dmc", action="store_true", help="Also benchmark depth-aware LR EB at boundary=None settings")
    args = parser.parse_args(argv)
    output = args.output.resolve()
    if args.summarize_only:
        summarize(output)
        return 0
    source = snapshot_source(args.repo.resolve(), output)
    os.environ.update(NUMBA_NUM_THREADS=str(args.workers), POLARS_MAX_THREADS=str(args.workers), OPENBLAS_NUM_THREADS="1",
                      OMP_NUM_THREADS="1", MPLBACKEND="Agg", PYTHONDONTWRITEBYTECODE="1",
                      NUMBA_CACHE_DIR=str(ROOT / "development/runtime/statistical_repair_numba"),
                      MPLCONFIGDIR=str(output / "runtime/mpl"))
    sys.path.insert(0, str(output / "source"))
    sys.path.insert(1, str(ROOT))
    import epykit as ep
    import polars as pl
    from epykit._regional_gls import run_regional_gls
    from scripts.run_epykit_pilot import table
    from scripts.compare_real_dmrs import WORKBOOK
    from scripts.run_real_gls import INPUT

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    ep.set_tmp_dir(output / "tmp")
    datasets = table(ROOT / "data/pilot/datasets.tsv") if args.kind in {"pilot", "both"} else []
    if args.datasets:
        known = {r["dataset"] for r in datasets}
        if set(args.datasets)-known:
            raise ValueError("Unknown selected pilot dataset")
        datasets = [r for r in datasets if r["dataset"] in args.datasets]
    jobs = [("pilot", row) for row in datasets]
    if args.kind in {"real", "both"}:
        jobs.append(("real", dict(dataset="GSE263850")))
    failed = False
    for kind, row in jobs:
        name = row["dataset"]
        previous = ROOT / "results" / ("epykit_real_gls_default" if kind == "real" else f"epykit_minimal_pilot/{name}")
        sheet = previous / "samplesheet.csv"
        obs = pl.read_csv(sheet).with_columns((pl.col("group") == ("KO" if kind == "real" else "treatment")).cast(pl.Int8).alias("treatment"))
        stores = {10: previous / "methylstore/.cache/filtered"}
        if kind == "real":
            stores[5] = ROOT / "results/gls_boundary_trial_real_stage2/cov5_filter/.cache/filtered"
        else:
            stores[5] = ROOT / f"results/epykit_gls_sensitivity/{name}/coverage_5/.cache/filtered"
        signatures = {}
        for config in args.configs:
            cfg = CONFIGS[config]
            cov = cfg["coverage_min"]
            store = stores[cov]
            if not store.is_dir():
                filtered = ep.MethylData(obs=obs, store=str(previous / "methylstore/.cache/raw"), assembly="GRCh38",
                                          analysis_root=str(output / f"derived_counts/{name}/coverage_{cov}"))
                ep.pp.filter_coverage(filtered, lo_count=cov, hi_perc=99.9)
                store = stores[cov] = Path(filtered.store)
            if cov not in signatures:
                signatures[cov] = store_identity(store, output / f"inputs/{name}_cov{cov}.json")
            verify_evaluation(source, output)
            identity = dict(evaluation_signature=source["evaluation_signature"], source_signature=source["source_signature"], dataset=name, kind=kind, config_id=config,
                            config=cfg, count_store_signature=signatures[cov], samplesheet_sha256=sha256(sheet),
                            chromosomes=args.chromosomes, n_perm=252 if kind == "pilot" else 100, seed=42, q_cutoff=.05,
                            include_dmc=args.include_dmc and cfg["boundary_min_effect"] is None)
            if kind == "pilot":
                identity["truth_sha256"] = pilot_truth_identity(ROOT / "data/pilot" / name)
            if kind == "real":
                identity.update(paper_sha256=sha256(WORKBOOK), refgene_sha256=sha256(INPUT.parent / "refseq/refGene.txt.gz"),
                                islands_sha256=sha256(INPUT.parent / "hg38_cpg_islands.bed"))
            run_dir = output / f"runs/{kind}/{name}/{config}"
            run_dir.mkdir(parents=True, exist_ok=True)
            if completed_run(run_dir, identity):
                logging.info("REUSE verified completed run %s/%s", name, config)
                continue
            state = dict(status="running", identity=identity, started_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
            write_json(run_dir / "run.json", state)
            logging.info("START %s %s %s; source %s", kind, name, config, source["commit"][:7])
            try:
                md = ep.MethylData(obs=obs, store=str(store), assembly="GRCh38", analysis_root=str(run_dir / "analysis"))
                started = time.perf_counter()
                _, candidates, params = run_regional_gls(
                    md, chromosomes=args.chromosomes, window_bp=500, min_cpgs=5, max_gap=500, min_length=50,
                    min_effect=.1, boundary_min_effect=cfg["boundary_min_effect"], n_perm=identity["n_perm"],
                    seed=42, q_cutoff=None, workers=args.workers, cache=False,
                    record_scan_scores=True,
                )
                elapsed = time.perf_counter()-started
                selected = candidates.filter(pl.col("combined_qvalue") <= .05)
                candidates.write_parquet(run_dir / "candidates.parquet")
                selected.write_parquet(run_dir / "selected.parquet")
                selected.write_csv(run_dir / "selected.tsv", separator="\t")
                write_json(run_dir / "params.json", params)
                metrics = (pilot_metrics(candidates, row, config, elapsed, params, run_dir) if kind == "pilot"
                           else real_metrics(ep, obs, store, candidates, config, elapsed, params, run_dir))
                if identity["include_dmc"]:
                    ep.pp.set_unite_type(md, type="union")
                    dmc_started = time.perf_counter()
                    ep.tl.dmc(md, test="lr", dispersion="eb", min_samples_treatment=2,
                              min_samples_control=2, chromosomes=args.chromosomes, tsv=False)
                    dmc = md.dmc
                    dmc.filter(pl.col("qvalue").is_finite() & (pl.col("qvalue") < .05)).write_parquet(run_dir / "dmc_selected.parquet")
                    metrics.update(dmc_metrics(dmc, ROOT / "data/pilot" / name if kind == "pilot" else None),
                                   dmc_seconds=time.perf_counter()-dmc_started)
                metrics.update(source_commit=source["commit"], source_signature=source["source_signature"],
                               chromosome_subset=args.chromosomes)
                verify_evaluation(source, output)
                if kind == "pilot" and identity["truth_sha256"] != pilot_truth_identity(ROOT / "data/pilot" / name):
                    raise ValueError("Pilot truth changed during evaluation")
                write_json(run_dir / "metrics.json", metrics)
                outputs = {p.name: sha256(p) for p in run_dir.iterdir() if p.is_file() and p.name != "run.json"}
                state.update(status="complete", outputs=outputs, completed_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
                logging.info("COMPLETE %s %s: %d candidates/%d selected in %.2fs", name, config, len(candidates), len(selected), elapsed)
            except Exception:
                failed = True
                state.update(status="failed", error=traceback.format_exc())
                logging.exception("FAILED %s %s", name, config)
            write_json(run_dir / "run.json", state)
            summarize(output)
    summarize(output)
    return int(failed)


if __name__ == "__main__":
    raise SystemExit(main())
