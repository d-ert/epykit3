# methylKit multicore timing (cores=8)

Host nproc: 24.  Re-run of calculateDiffMeth with mc.cores=8.
Single-thread baseline = committed <seed>/methylkit.tsv.timing.tsv.

| seed | phase | wall_s (cores=8) | cpu_s | wall_s (1 core) |
|------|-------|----------------------|-------|-----------------|
| 2026000 | read | 1.014 | 1.045 | 1.172 |
| 2026000 | unite | 0.0960000000000001 | 0.64 | 0.548999999999999 |
| 2026000 | diffmeth | 12.318 | 3.711 | 72.627 |
| 2026001 | read | 0.944 | 0.975 | 1.173 |
| 2026001 | unite | 0.0890000000000004 | 0.599 | 0.520000000000001 |
| 2026001 | diffmeth | 12.724 | 3.648 | 71.832 |
| 2026002 | read | 0.955 | 0.986 | 0.994999999999999 |
| 2026002 | unite | 0.0909999999999993 | 0.598000000000001 | 0.683 |
| 2026002 | diffmeth | 12.32 | 3.76 | 76.412 |
