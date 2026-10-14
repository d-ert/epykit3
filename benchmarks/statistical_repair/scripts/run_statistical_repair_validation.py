#!/usr/bin/env python3
"""Independent null GLS/DMC validation on a generated simulator suite."""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import time
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.run_statistical_repair_benchmarks import (
    CONFIGS, REPAIR_REPO, completed_run, sha256, snapshot_source, store_identity, verify_evaluation, write_json,
)


def null_validation_summary(records, expected):
    from wgbs_sim.evaluate import null_summary

    if len(records) != expected or any(r["status"] != "complete" for r in records):
        raise ValueError("Every declared null run must be complete before reporting validation")
    if any(not r["all_null"] for r in records):
        raise ValueError("Null validation contains nonnull evidence")
    result = {}
    for config, track in sorted({(r["config_id"], r["track"]) for r in records}):
        selected = [r for r in records if (r["config_id"], r["track"]) == (config, track)]
        engines = ["gls", "dmc"]
        if all(r.get("gls_maxT_calls") is not None for r in selected):
            engines.append("gls_maxT")
        result[f"{config}/{track}"] = {
            engine: {**null_summary([r[f"{engine}_calls"] for r in selected]),
                     "discovery_counts": [r[f"{engine}_calls"] for r in selected]}
            for engine in engines
        }
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repo", type=Path, default=REPAIR_REPO)
    parser.add_argument("--configs", nargs="+", choices=list(CONFIGS), default=["cov10_boundNone", "cov5_bound05"])
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args(argv)
    datasets = json.loads((args.dataset_root / "suite.json").read_text())["runs"]
    if any(not r["all_null"] or not r["validated"] for r in datasets):
        raise ValueError("Need a complete, validated all-null suite")
    output = args.output.resolve()
    source = snapshot_source(args.repo.resolve(), output, entrypoint=Path(__file__))
    os.environ.update(NUMBA_NUM_THREADS=str(args.workers), POLARS_MAX_THREADS=str(args.workers), OPENBLAS_NUM_THREADS="1",
                      OMP_NUM_THREADS="1", MPLBACKEND="Agg", PYTHONDONTWRITEBYTECODE="1",
                      NUMBA_CACHE_DIR=str(ROOT / "development/runtime/statistical_repair_numba"),
                      MPLCONFIGDIR=str(output / "runtime/mpl"))
    sys.path.insert(0, str(output / "source"))
    import epykit as ep
    import polars as pl
    from epykit._regional_gls import run_regional_gls
    from scripts.run_epykit_pilot import write_samplesheet

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    ep.set_tmp_dir(output / "tmp")
    records = []
    for dataset in datasets:
        name = dataset["dataset"]
        data = args.dataset_root / name
        derived = output / f"derived_counts/{name}"
        derived.mkdir(parents=True, exist_ok=True)
        sheet = derived / "samplesheet.csv"
        if not sheet.exists():
            write_samplesheet(data, sheet)
        # Raw ingestion and filtering stay in this run's own derived directory.
        raw = ep.read_bismark(str(sheet), treatment_group="treatment", control_group="control", assembly="GRCh38",
                              coordinate_base="one_based", store_dir=str(derived / "methylstore"))
        for config in args.configs:
            cfg = CONFIGS[config]
            analysis = output / f"runs/{name}/{config}"
            analysis.mkdir(parents=True, exist_ok=True)
            md = ep.MethylData(obs=raw.obs, store=raw.store, assembly="GRCh38", analysis_root=str(analysis / "analysis"))
            ep.pp.filter_coverage(md, lo_count=cfg["coverage_min"], hi_perc=99.9)
            verify_evaluation(source, output)
            identity = dict(evaluation_signature=source["evaluation_signature"], source_signature=source["source_signature"], dataset=name, config=cfg,
                            manifest_sha256=sha256(data / "manifest.json"), samplesheet_sha256=sha256(sheet),
                            suite_sha256=sha256(args.dataset_root / "suite.json"), n_perm=252, seed=42,
                            store_signature=store_identity(Path(md.store), analysis / "input_store.json"),
                            dmc=dict(test="lr", dispersion="eb", min_samples_treatment=2, min_samples_control=2))
            if completed_run(analysis, identity):
                records.append(json.loads((analysis / "metrics.json").read_text()))
                continue
            state = dict(status="running", identity=identity)
            write_json(analysis / "run.json", state)
            logging.info("START null validation %s %s", name, config)
            metrics = dict(status="failed", config_id=config, track=dataset["track"], seed=dataset["seed"],
                           dataset=name, all_null=True, source_commit=source["commit"])
            try:
                started = time.perf_counter()
                _, candidates, params = run_regional_gls(
                    md, window_bp=500, min_cpgs=5, max_gap=500, min_length=50, min_effect=.1,
                    boundary_min_effect=cfg["boundary_min_effect"], n_perm=252, seed=42, q_cutoff=None,
                    workers=args.workers, cache=False,
                    record_scan_scores=True,
                )
                selected = candidates.filter(pl.col("combined_qvalue") <= .05)
                candidates.write_parquet(analysis / "candidates.parquet")
                selected.write_parquet(analysis / "selected.parquet")
                write_json(analysis / "params.json", params)
                gls_elapsed = time.perf_counter()-started
                ep.pp.set_unite_type(md, type="union")
                started = time.perf_counter()
                ep.tl.dmc(md, test="lr", dispersion="eb", min_samples_treatment=2, min_samples_control=2, tsv=False)
                dmc = md.dmc
                dmc_calls = len(dmc.filter(pl.col("qvalue").is_finite() & (pl.col("qvalue") < .05)))
                dmc.write_parquet(analysis / "dmc.parquet")
                metrics.update(status="complete", gls_calls=len(selected), dmc_calls=dmc_calls,
                               gls_maxT_calls=len(candidates.filter(pl.col("maxT_pvalue") <= .05)) if "maxT_pvalue" in candidates.columns else None,
                               gls_seconds=gls_elapsed, dmc_seconds=time.perf_counter()-started,
                               n_candidates=len(candidates), n_null_regions=params["null_regions"])
                verify_evaluation(source, output)
                write_json(analysis / "metrics.json", metrics)
                state.update(status="complete", outputs={p.name: sha256(p) for p in analysis.iterdir()
                                                            if p.is_file() and p.name != "run.json"})
                logging.info("COMPLETE null %s %s: GLS %d, DMC %d", name, config, len(selected), dmc_calls)
            except Exception:
                state.update(status="failed", error=traceback.format_exc())
                logging.exception("FAILED null validation %s %s", name, config)
            write_json(analysis / "run.json", state)
            records.append(metrics)
            write_json(output / "per_run.json", records)
    summary = null_validation_summary(records, len(datasets)*len(args.configs))
    write_json(output / "summary.json", summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
