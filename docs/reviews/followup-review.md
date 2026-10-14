# Follow-up review: d91378f and workspace provenance fixes

Verdict: **the three original Important findings are resolved; one additional Important null-validator failure-path issue remains**. No Critical findings. No source or active benchmark files were modified during this review.

Reviewed the supplied `review-followup.diff` first, the changed source and regression tests, current workspace benchmark/null-validation runners, and their provenance tests. Accepted the reported 96 passing source tests and focused harness evidence rather than rerunning the unchanged suite. The final 72 pilot configurations, four real configurations, and holdout validation remain live; this review makes no performance or calibration claim.

## Original findings resolved

1. **Exact-fit symmetry:** the shared observed-design residual check is mathematically appropriate: covariance changes do not change the fixed-effect column space, and explicit OLS residuals avoid the prior subtractive-cancellation issue. Its tolerance scales with response energy and uses squared floating-point precision rather than imposing a biological variance floor. Complete, masked, fixed-covariance, profiled, and batch entry paths now reject exact fits consistently. Regression coverage includes the original k=3/5/8/20 count examples and complements; small nonzero residuals are preserved across response scales.
2. **Pooled covariance/maxT:** the selector rejects the incompatible combination before inference. The pooled diagnostic is now NaN, with `maxT_available=False` and an explicit unavailable scope; the default replay path retains its original semantics. Cache source hashing and the revision bump prevent reuse of the old behavior.
3. **Evaluation provenance:** resume now checks the actual entry point, runner/evaluator hashes, all workspace Python evaluation/helper dependencies, and frozen copies. Per-run identity carries the evaluation signature. Both pilot truth files are hashed and checked again at completion. Null-validation provenance now identifies its own entry point, suite, and sample sheet. These changes close the prior resume-mixing defect.

## Important — Final null-validation verification can fail while the summary reports completion

Location: workspace `scripts/run_statistical_repair_validation.py:119-135`.

`metrics.update(status="complete", ...)` runs before the final `verify_evaluation`, metrics write, and output hashing. If any of those operations raises, the exception handler changes only `state` to failed. It leaves `metrics["status"]` complete, appends that record, and calls `null_validation_summary`, which checks the metric status rather than the persisted run state. A failed provenance check can therefore produce a successful exit and an accepted null-error summary.

This is especially relevant to the new provenance guarantee: detecting a code change is insufficient if the failed run still enters the aggregate. Successful runs are not invalidated by the existence of this failure path. Existing final outputs can be accepted only after independently checking every manifest is complete and its hashes match.

### Confirmed isolated actual-main reproduction

The following probe uses temporary data and mocks the expensive inference/ingestion, but executes the real validator's orchestration, exception handling, manifest writing, and summary functions. It does not modify root scripts or running benchmarks.

```python
import json, tempfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import polars as pl
import epykit as ep
from epykit import _regional_gls
from scripts import run_statistical_repair_validation as validator

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    data, output = root / "data", root / "out"
    name = "empirical_null_seed1"
    (data / name).mkdir(parents=True)
    (data / "suite.json").write_text(json.dumps({"runs": [dict(
        dataset=name, track="empirical", seed=1, all_null=True, validated=True
    )]}))
    (data / name / "manifest.json").write_text("{}")
    sheet = output / "derived_counts" / name / "samplesheet.csv"
    sheet.parent.mkdir(parents=True)
    sheet.write_text("fixture")
    fake_md = SimpleNamespace(
        obs=None, store=str(root / "store"),
        dmc=pl.DataFrame({"qvalue": [.9]}),
    )
    candidate = pl.DataFrame({"combined_qvalue": [.9], "maxT_pvalue": [.9]})
    source = dict(source_signature="fixture", evaluation_signature="fixture", commit="fixture")
    with (
        patch.object(validator, "snapshot_source", return_value=source),
        patch.object(validator, "store_identity", return_value="fixture"),
        patch.object(validator, "verify_evaluation", side_effect=[
            None, ValueError("Live evaluation dependency changed: fixture")
        ]),
        patch.object(ep, "read_bismark", return_value=fake_md),
        patch.object(ep, "MethylData", return_value=fake_md),
        patch.object(ep, "set_tmp_dir"),
        patch.object(ep.pp, "filter_coverage"),
        patch.object(ep.pp, "set_unite_type"),
        patch.object(ep.tl, "dmc"),
        patch.object(_regional_gls, "run_regional_gls", return_value=(
            candidate, candidate, {"null_regions": 0}
        )),
    ):
        result = validator.main([
            "--dataset-root", str(data), "--output", str(output),
            "--configs", "cov10_boundNone",
        ])
    print("RETURN_CODE", result)
    manifest = output / "runs" / name / "cov10_boundNone" / "run.json"
    print("RUN_STATE", json.loads(manifest.read_text())["status"])
    print("METRIC_STATE", json.loads((output / "per_run.json").read_text())[0]["status"])
    print("SUMMARY_EXISTS", (output / "summary.json").exists())
```

Observed output after the validator logs the injected provenance failure:

```text
RETURN_CODE 0
RUN_STATE failed
METRIC_STATE complete
SUMMARY_EXISTS True
```

Required repair: synchronize metric status with failed run status in the exception path, or delay marking metrics complete until all final verification/writes have succeeded. A regression should inject a failure at the final verification call and assert the metric record is failed, final completeness validation rejects the set, and no successful summary is created. Apply the repair in an isolated archived harness while active processes retain their frozen scripts; synchronize the live harness only after those processes terminate, preserving exact executed-script provenance.

## Acceptance status

Inference fixes in `d91378f`: accepted for the three reviewed concerns. Workspace provenance fix: accepted for identity enforcement; null-validator failure handling requires the additional repair above. Final result completeness/hash verification, final performance interpretation, and bundle acceptance remain pending terminal benchmark evidence.
