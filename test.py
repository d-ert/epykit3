#Basic end-to-end WGBS methylation analysis workflow using epykit.

import epykit as ep

# 1. Ingestion: Bismark .cov.gz -> Partitioned Parquet store
md = ep.read_bismark(
    samplesheet="samplesheet.csv",
    treatment_group="Het_AKAP11_KO",
    control_group="WT",
    store_dir="results/methylstore",
)

# 2. Quality Control
ep.tl.qc(md)

# 3. Preprocessing (Coverage filtering & site intersection)
ep.pp.filter_coverage(md, lo_count=10, hi_perc=99.9)
ep.pp.set_unite_type(md, type="union")  # or "union" + min_samples_* guards

# 4. Differentially Methylated Cytosines (DMC)
# Default 'lr' (quasi-binomial LR); change to 'beta_binomial' to test the new engine
ep.tl.dmc(md, test="lr")

# 5. Differentially Methylated Regions (DMR)
ep.tl.dmr(
    md,
    method="regional_gls",
    window_bp=500,        # Smoothing window size in bp for region discovery (default: 500)
    min_cpgs=5,           # Minimum CpGs per DMR (default: 5)
    dis_merge_bp=500,     # Max gap to merge adjacent candidate regions (default: 500)
    min_abs_meth_diff=0.1,# Effect size floor for candidate seeding (default: 0.1)
    n_perm=20,            # Number of label permutations (default: 100; reduce for quick test)
    perm_seed=42,         # Random seed for reproducible permutations
    min_mean_qvalue=0.05, # FDR cutoff on combined_qvalue (or None to keep all)
)

# 6. Annotation (Optional: checks if GTF/CpG files exist)
# ep.tl.annotate(md, gtf="path/to/genes.gtf", cpg_islands="path/to/islands.bed")

# 7. Export result tables (TSV)
md.export_tables("results/tables", full=True)

# 8. Interactive HTML Report
ep.generate_report(md, output="results/report.html", title="GSE263850 Analysis Report")
