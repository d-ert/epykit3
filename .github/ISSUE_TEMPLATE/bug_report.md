---
name: Bug report
about: Something produced an error or a result you believe is wrong
labels: bug
---

## What happened

<!-- One or two sentences. Paste the full traceback in the block below. -->

```text

```

## Minimal reproduction

<!-- The smallest script or CLI call that shows the problem. If it needs data,
     say which input format (Bismark .cov, MethylDackel .bedGraph, combined-
     strand BED) and roughly how many samples and CpGs. -->

```python

```

## Expected

## Environment

- epykit version (`python -c "import epykit; print(epykit.__version__)"`):
- Python version and OS:
- Install method (`pip`, `uv`, conda) and extras installed:
- `polars`, `numpy`, `scipy`, `statsmodels` versions (`uv pip list` or `pip list`):

## Checklist

- [ ] I ran with `-v` (CLI) or `logging.basicConfig(level="INFO")` (API) and included relevant log lines above.
- [ ] If coordinates look off by one: my input is Bismark `.cov` (1-based, `start == end`) or MethylDackel/BED (0-based). Stated which.
