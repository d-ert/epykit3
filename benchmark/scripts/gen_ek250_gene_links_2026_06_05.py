"""Generate the 100 kb DMR-gene link table for the post-rerun
ek-chain_merge dis.merge=250 DMR set (1,139 DMRs), matching the schema
of the committed dmr_gene_links_100kb.csv files.

Rule (paper Methods, identical to the other callers' link tables):
for each DMR midpoint, link every gene whose canonical TSS is within
+/- 100 kb of the midpoint. Canonical TSS per gene = most upstream
txStart (+ strand) / txEnd (- strand) across that gene's refGene rows.
"""

from __future__ import annotations
import gzip
import sys
import io
from pathlib import Path
import polars as pl

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

BENCH = Path(__file__).resolve().parents[1]
REFGENE = Path(r"D:/Coding/Projeler/methyl_lib/epykit2/GSE263850_RAW/refseq/refGene.txt.gz")
DMR_CSV = BENCH / "data" / "multi_thread_and_chain_sweep" / "chain_merge_dis_merge_sweep" / "dis_merge_250" / "dmr.csv"
OUT = DMR_CSV.parent / "dmr_gene_links_100kb.csv"
WINDOW = 100_000


def load_canonical_tss():
    """Return {chrom: [(gene, tss), ...]} with one canonical TSS per gene."""
    # gene -> (chrom, strand, most_upstream_tss)
    best = {}
    with gzip.open(REFGENE, "rt") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            # refGene: 0 bin,1 acc,2 chrom,3 strand,4 txStart,5 txEnd,...,12 name2
            chrom, strand = p[2], p[3]
            tx_start, tx_end = int(p[4]), int(p[5])
            gene = p[12]
            if "_" in chrom:  # skip alt/random contigs
                continue
            tss = tx_start if strand == "+" else tx_end
            key = gene
            if key not in best:
                best[key] = (chrom, strand, tss)
            else:
                _, _, prev = best[key]
                # most upstream: smaller tss on +, larger tss on -
                if strand == "+" and tss < prev:
                    best[key] = (chrom, strand, tss)
                elif strand == "-" and tss > prev:
                    best[key] = (chrom, strand, tss)
    by_chrom = {}
    for gene, (chrom, strand, tss) in best.items():
        by_chrom.setdefault(chrom, []).append((gene, tss))
    for chrom in by_chrom:
        by_chrom[chrom].sort(key=lambda t: t[1])
    return by_chrom


def main():
    tss_by_chrom = load_canonical_tss()
    dmr = pl.read_csv(DMR_CSV)
    rows = []
    for idx, r in enumerate(dmr.iter_rows(named=True)):
        chrom = r["chrom"]
        mid = (r["start"] + r["end"]) // 2
        for gene, tss in tss_by_chrom.get(chrom, []):
            d = tss - mid
            if abs(d) <= WINDOW:
                rows.append({
                    "dmr_index": idx,
                    "chrom": chrom,
                    "dmr_start": r["start"],
                    "dmr_end": r["end"],
                    "dmr_midpoint": mid,
                    "dmr_type": r["dmr_type"],
                    "mean_meth_diff": r["mean_meth_diff"],
                    "gene": gene,
                    "tss": tss,
                    "distance_tss_to_midpoint": d,
                    "abs_distance": abs(d),
                })
    out = pl.DataFrame(rows)
    out.write_csv(OUT)
    print(f"wrote {OUT}")
    print(f"  {out.height} (DMR, gene) pairs; {out['gene'].n_unique()} unique genes; "
          f"{out['dmr_index'].n_unique()} DMRs with >=1 link")


if __name__ == "__main__":
    main()
