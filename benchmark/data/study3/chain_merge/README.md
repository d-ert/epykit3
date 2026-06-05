# Chain-merge replication of GSE263850 (paper-matched parameters)

This folder contains a clean, single-run reproduction of the differential
methylation analysis in Farhangdoost et al. (2024) on GSE263850, using
epykit's `dmr_chain_merge` engine with parameters mapped one-for-one to
the paper's DSS::callDMR call. See `parameters.json` for the full mapping.

## Files

| File | Description |
|---|---|
| `dmr_chain_merge.parquet` | All DMRs (chain_merge, alpha=1e-5, BH q < 0.05), HOMER-style annotated. Polars/Pandas-friendly. |
| `dmr_chain_merge.csv` | Same as above, CSV. |
| `dmr_gene_links_100kb.csv` | Long-form (DMR, gene) pairs where the gene's canonical TSS is within 100 kb of the DMR midpoint. Matches the paper's exact gene-linkage rule. |
| `parameters.json` | Exact parameters used, with DSS <-> epykit mapping. |
| `summary.md` | Headline numbers: DMR counts, morphology, paper-coord overlap, gene linkage. |
| `run_log.txt` | INFO log of the run. |
| `_store/` | Per-run epykit MethylStore (CpG counts, DMC parquet). Safe to delete; will be regenerated. |
| `_tmp/` | Scratch space; safe to delete. |

## Reproducing

```
py benchmark/scripts/run_chain_merge_replication.py
```

Single command. Deterministic (modulo floating-point summation order).
Total wall time: ~8-10 minutes on the local box.

## Headline numbers

See `summary.md` for the full table. Highlights:

- DMRs (chain_merge, alpha=1e-5): **852**
- Paper DMRs (DSS, p.threshold=1e-5): **NA**
- Coordinate recall of paper DMRs: **0.0%**
- Genes within 100 kb of any DMR: **1,290**
