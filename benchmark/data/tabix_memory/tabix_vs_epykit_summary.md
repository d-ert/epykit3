# methylKit tabix vs RAM vs epykit lr — summary

Out-of-core head-to-head. tabix changes only storage/IO, not statistics (per-CpG p/q identical to RAM), so the axes are wall time and peak memory. USS (unique set size) is the fair memory metric under `mc.cores>1` because summed RSS double-counts copy-on-write shared pages across forked workers.

### gse263850

| tool | backend | cores | wall_s | RSS MB | USS MB | n_sites | speedup vs epykit | RSS× vs epykit | USS× vs epykit | tabix/ram RSS | tabix/ram wall |
|---|---|---|---|---|---|---|---|---|---|---|---|
| epykit_lr | parquet_ooc | 1 | 80.32 | 5890.7 | 5833.3 | 22020939 | 1.0× | 1.0× | 1.0× | — | — |

### simulator

| tool | backend | cores | wall_s | RSS MB | USS MB | n_sites | speedup vs epykit | RSS× vs epykit | USS× vs epykit | tabix/ram RSS | tabix/ram wall |
|---|---|---|---|---|---|---|---|---|---|---|---|
| epykit_lr | parquet_ooc | 1 | 1.88 | 435.9 | 416.8 | 100000 | 1.0× | 1.0× | 1.0× | — | — |
| methylkit_ram | ram | 1 | 78.15 | 1297.6 | 1279.0 | 99999 | 41.6× | 2.98× | 3.07× | — | — |
| methylkit_tabix | tabix_ooc | 1 | 79.18 | 1295.8 | 1277.1 | 99999 | 42.1× | 2.97× | 3.06× | 1.0× | 1.01× |
| methylkit_ram | ram | 8 | 21.18 | 8785.1 | 3165.3 | 99999 | 11.3× | 20.15× | 7.59× | — | — |
| methylkit_tabix | tabix_ooc | 8 | 25.14 | 10063.8 | 9075.6 | 99999 | 13.4× | 23.09× | 21.77× | 1.15× | 1.19× |
