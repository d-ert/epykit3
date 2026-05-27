"""One-shot finisher: load the existing chain_merge parquet and write the
remaining structured outputs (CSV, links, parameters.json, summary.md,
README.md). Use after a successful DMC+DMR run if write_outputs() failed
mid-way for any reason. Idempotent.
"""

from __future__ import annotations

import sys
from pathlib import Path

import polars as pl

# Reuse the helpers from the main script.
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from run_chain_merge_replication import (  # noqa: E402
    OUT_DIR, REFGENE, GENE_LINK_MAX_BP,
    parse_refgene_tss, build_100kb_gene_links,
    compare_to_paper, write_outputs, log,
)


def main() -> None:
    parquet = OUT_DIR / "dmr_chain_merge.parquet"
    log.info("Loading existing DMR parquet: %s", parquet)
    annotated = pl.read_parquet(str(parquet))
    log.info("    %d rows, %d cols", annotated.height, len(annotated.columns))

    log.info("[gene-link] Building 100 kb DMR-gene table from %s", REFGENE)
    gene_tss = parse_refgene_tss(REFGENE)
    log.info("    %d (chrom, gene) TSS entries", len(gene_tss))
    links = build_100kb_gene_links(annotated, gene_tss, GENE_LINK_MAX_BP)
    log.info("    %d (DMR, gene) links built; %d unique genes",
             len(links), links["gene"].nunique())

    paper_stats = compare_to_paper(annotated)
    if paper_stats.get("paper_table_available"):
        log.info(
            "[compare] paper=%d, ours=%d, recall=%.1f%%, precision=%.1f%%",
            paper_stats["paper_n_dmr"], paper_stats["our_n_dmr"],
            paper_stats["coord_recall_of_paper"] * 100,
            paper_stats["coord_precision"] * 100,
        )

    # Pipeline stats: we don't have DMC counts cached, but they're in the
    # log file. Pull them from the parquet for the relevant fields and
    # parse the log for the DMC numbers.
    n_dmr = annotated.height
    log_path = OUT_DIR / "run_log.txt"
    n_cpgs_tested = n_dmc_sig = 0
    dmc_seconds = dmr_seconds = 0.0
    if log_path.exists():
        for line in log_path.read_text(encoding="utf-8").splitlines():
            if "CpGs tested" in line:
                # "21,993,377 CpGs tested; 12,336 sig at q<0.05  (219.8s)"
                try:
                    parts = line.split("DMC: ")[1]
                    n_cpgs_tested = int(parts.split(" CpGs")[0].replace(",", ""))
                    n_dmc_sig = int(parts.split("; ")[1].split(" sig")[0].replace(",", ""))
                    dmc_seconds = float(parts.split("(")[1].split("s)")[0])
                except Exception:
                    pass
            if "candidates after BH" in line:
                try:
                    dmr_seconds = float(line.split("(")[1].split("s)")[0])
                except Exception:
                    pass
    pipeline_stats = dict(
        n_cpgs_tested=n_cpgs_tested,
        n_dmc_sig=n_dmc_sig,
        n_dmr=int(n_dmr),
        dmc_seconds=round(dmc_seconds, 1),
        dmr_seconds=round(dmr_seconds, 1),
    )
    log.info("Pipeline stats (parsed from log): %s", pipeline_stats)

    write_outputs(annotated, links, pipeline_stats, paper_stats)
    log.info("Finisher complete. Outputs in %s", OUT_DIR)


if __name__ == "__main__":
    main()
