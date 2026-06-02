// epykit — thesis-progress presentation
// Audience: general bioinformatics (master's advisor)
// Outputs: benchmark/epykit_thesis_presentation.pptx

const path = require("path");
process.env.NODE_PATH = require("child_process")
  .execSync("npm root -g")
  .toString()
  .trim();
require("module").Module._initPaths();

const pptxgen = require("pptxgenjs");

const ROOT = path.resolve(__dirname, "..");
const FIG = (relpath) => path.join(ROOT, "figures", relpath);
const HEAD = (name) => FIG(path.join("summary", name));

const C = {
  navy: "0F2942",
  navyDeep: "081C2E",
  blue: "065A82",
  teal: "1C7293",
  gold: "F4B400",
  white: "FFFFFF",
  ink: "1E293B",
  muted: "64748B",
  rule: "CBD5E1",
  cardBg: "F8FAFC",
  cardBg2: "EFF6FB",
  good: "22C55E",
  warn: "F97316",
  bad:  "DC2626",
};

const F = {
  head: "Georgia",
  body: "Calibri",
  mono: "Consolas",
};

const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE";
pres.author = "Deniz Ertürk";
pres.title  = "epykit v1.0 — thesis progress";

const W = 13.3;
const H = 7.5;

function titleHeader(s, title, subtitle) {
  s.addText(title, {
    x: 0.6, y: 0.32, w: W - 1.2, h: 0.7,
    fontSize: 28, fontFace: F.head, bold: true, color: C.navy, margin: 0,
  });
  if (subtitle) {
    s.addText(subtitle, {
      x: 0.6, y: 0.92, w: W - 1.2, h: 0.45,
      fontSize: 13, fontFace: F.body, italic: true, color: C.muted, margin: 0,
    });
  }
}

function pageFooter(s, label, page) {
  s.addText(label, {
    x: 0.6, y: H - 0.42, w: 8, h: 0.3,
    fontSize: 9, fontFace: F.body, color: C.muted, margin: 0,
  });
  s.addText(String(page), {
    x: W - 1.2, y: H - 0.42, w: 0.6, h: 0.3,
    fontSize: 9, fontFace: F.body, color: C.muted, align: "right", margin: 0,
  });
}

function speakerNote(s, text) {
  s.addNotes(text);
}

// =========================================================================
// 1 — Title
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.navyDeep };

  s.addShape(pres.shapes.RECTANGLE, { x: 0.6, y: 6.4, w: 0.35, h: 0.35,
    fill: { color: C.teal }, line: { type: "none" } });
  s.addShape(pres.shapes.RECTANGLE, { x: 1.05, y: 6.4, w: 0.35, h: 0.35,
    fill: { color: C.blue }, line: { type: "none" } });
  s.addShape(pres.shapes.RECTANGLE, { x: 1.5, y: 6.4, w: 0.35, h: 0.35,
    fill: { color: C.gold }, line: { type: "none" } });

  s.addText("epykit", {
    x: 0.6, y: 1.6, w: 12, h: 1.3,
    fontSize: 96, fontFace: F.head, bold: true, color: C.white, margin: 0,
  });
  s.addText("A Python-native pipeline for differential DNA methylation analysis", {
    x: 0.6, y: 2.95, w: 12, h: 0.6,
    fontSize: 22, fontFace: F.body, color: "BFD7E5", margin: 0,
  });
  s.addText("Thesis progress  •  benchmarked on simulated and real WGBS data", {
    x: 0.6, y: 3.55, w: 12, h: 0.4,
    fontSize: 15, fontFace: F.body, italic: true, color: "8FB6CA", margin: 0,
  });

  const stats = [
    { num: "v1.0.0",   label: "stable API, MIT" },
    { num: "12.2×",    label: "faster than methylKit\non real WGBS" },
    { num: "3.8×",     label: "less peak RAM\non real WGBS" },
  ];
  stats.forEach((stat, i) => {
    const y = 4.6 + i * 0.78;
    s.addText(stat.num, {
      x: 8.0, y: y, w: 2.0, h: 0.6,
      fontSize: 32, fontFace: F.head, bold: true, color: C.gold,
      align: "right", margin: 0,
    });
    s.addText(stat.label, {
      x: 10.2, y: y, w: 2.6, h: 0.65,
      fontSize: 12, fontFace: F.body, color: "BFD7E5", margin: 0,
    });
  });

  s.addText("Deniz Ertürk  •  Master's thesis progress  •  June 2026", {
    x: 0.6, y: 6.95, w: 8, h: 0.3,
    fontSize: 11, fontFace: F.body, color: "8FB6CA", margin: 0,
  });

  speakerNote(s,
    "Opening: introduce yourself, the project name, and the headline. " +
    "epykit reached v1 last week. Three numbers worth remembering for the rest of the talk: " +
    "12.2× speedup on real data, 3.8× less RAM, and v1.0 with a stable API."
  );
}

// =========================================================================
// 2 — Biological motivation
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(s,
    "Why DNA methylation matters",
    "A reversible epigenetic mark with direct biomedical consequence"
  );

  s.addText([
    { text: "DNA methylation (5-methylcytosine, 5mC) ", options: { bold: true, color: C.ink } },
    { text: "is a chemical modification of CpG dinucleotides that silences gene transcription without changing the DNA sequence.", options: { color: C.ink, breakLine: true } },
    { text: "", options: { breakLine: true } },
    { text: "It is implicated in:", options: { bold: true, color: C.ink, breakLine: true } },
    { text: "Cancer — global hypomethylation + focal hypermethylation of tumour-suppressor promoters is a hallmark of malignant transformation.", options: { bullet: true, color: C.ink, breakLine: true } },
    { text: "Imprinting and X-chromosome inactivation — parent-of-origin-specific silencing.", options: { bullet: true, color: C.ink, breakLine: true } },
    { text: "Neuropsychiatric disease — schizophrenia, bipolar disorder, autism all show methylation-quantitative-trait-locus enrichment.", options: { bullet: true, color: C.ink, breakLine: true } },
    { text: "Aging — \"epigenetic clocks\" (Horvath, Hannum) estimate biological age from methylation at ~ 350 CpGs.", options: { bullet: true, color: C.ink } },
  ], { x: 0.6, y: 1.55, w: 7.4, h: 5.4, fontSize: 14, fontFace: F.body,
       paraSpaceAfter: 6, margin: 0 });

  s.addShape(pres.shapes.RECTANGLE, { x: 8.3, y: 1.7, w: 4.4, h: 5.0,
    fill: { color: C.cardBg }, line: { color: C.rule, width: 1 } });
  s.addText("Why measure it?", {
    x: 8.5, y: 1.85, w: 4, h: 0.4, fontSize: 14, fontFace: F.head, bold: true, color: C.navy, margin: 0,
  });
  s.addText([
    { text: "Most clinical biomarkers (eg blood-tumour signature) come from finding ", options: { color: C.ink } },
    { text: "differentially methylated CpGs (DMCs)", options: { bold: true, color: C.navy } },
    { text: " and the regions that contain them ", options: { color: C.ink } },
    { text: "(DMRs)", options: { bold: true, color: C.navy } },
    { text: " between healthy and disease groups.", options: { color: C.ink, breakLine: true } },
    { text: "", options: { breakLine: true } },
    { text: "The statistics is non-trivial: per-CpG read counts are over-dispersed binomial, sample sizes are small (often n ≤ 6), and the genome holds ≈ 28 M CpGs.", options: { italic: true, color: C.muted } },
  ], { x: 8.5, y: 2.3, w: 4, h: 4.3, fontSize: 11, fontFace: F.body, margin: 0 });

  pageFooter(s, "Framing — why this matters", 2);
  speakerNote(s,
    "Ground the work biologically. The advisor may not know that methylation is the basis " +
    "for several FDA-approved diagnostic tests (eg ColoGuard) and for epigenetic age. " +
    "Stress: the statistics is the hard part, which is what my package addresses."
  );
}

// =========================================================================
// 3 — WGBS pipeline primer
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(s,
    "From biological sample to differential calls — the WGBS pipeline",
    "Where epykit sits in the data-processing chain"
  );

  // Pipeline boxes
  const stages = [
    { label: "Bisulfite-treat\nDNA", desc: "Convert unmethylated C→T\nMethylated C unchanged",         color: "FED7AA" },
    { label: "Sequence",             desc: "Illumina WGBS\n~30× per base genome-wide",                 color: "FED7AA" },
    { label: "Align",                desc: "Bismark / BWA-meth\nto bisulfite reference",               color: "FDBA74" },
    { label: "Call methylation",     desc: "MethylDackel / Bismark\n.cov / .bedGraph",                 color: "FB923C" },
    { label: "DMC / DMR\nanalysis",  desc: "epykit • methylKit\nDSS • RADMeth",                        color: "F97316" },
    { label: "Interpret",            desc: "Gene annotation\nenrichment • report",                     color: "0EA5E9" },
  ];

  const x0 = 0.6;
  const y0 = 2.5;
  const boxW = 1.95;
  const boxH = 1.6;
  const gap = 0.10;

  stages.forEach((st, i) => {
    const x = x0 + i * (boxW + gap);
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, {
      x, y: y0, w: boxW, h: boxH,
      fill: { color: st.color }, line: { type: "none" },
      rectRadius: 0.10,
    });
    s.addText(st.label, {
      x, y: y0 + 0.08, w: boxW, h: 0.6,
      fontSize: 12, fontFace: F.body, bold: true, color: C.ink, align: "center", margin: 0,
    });
    s.addText(st.desc, {
      x, y: y0 + 0.7, w: boxW, h: 0.85,
      fontSize: 9, fontFace: F.body, color: C.ink, align: "center", margin: 0,
    });
    if (i < stages.length - 1) {
      const ax = x + boxW + 0.005;
      s.addText("▶", {
        x: ax, y: y0 + 0.6, w: 0.1, h: 0.4,
        fontSize: 12, fontFace: F.body, color: C.muted, margin: 0,
      });
    }
  });

  // Highlight DMC/DMR stage
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, {
    x: x0 + 4 * (boxW + gap) - 0.05, y: y0 - 0.05, w: boxW + 0.10, h: boxH + 0.10,
    fill: { type: "none" }, line: { color: C.gold, width: 2.5 },
    rectRadius: 0.12,
  });
  s.addText("epykit operates here", {
    x: x0 + 4 * (boxW + gap) - 0.4, y: y0 + boxH + 0.15, w: boxW + 0.8, h: 0.3,
    fontSize: 10, fontFace: F.body, italic: true, color: C.gold, align: "center", bold: true, margin: 0,
  });

  s.addText([
    { text: "Bisulfite treatment ", options: { bold: true, color: C.ink } },
    { text: "is the core trick: it deaminates ", options: { color: C.ink } },
    { text: "unmethylated cytosines into uracils ", options: { italic: true, color: C.navy } },
    { text: "(which read as T) while leaving ", options: { color: C.ink } },
    { text: "methylated cytosines untouched ", options: { italic: true, color: C.navy } },
    { text: "(still read as C). The C/T ratio at each genomic CpG position then encodes the methylation fraction (β) in that sample.", options: { color: C.ink } },
  ], { x: 0.6, y: 4.55, w: W - 1.2, h: 1.0, fontSize: 13, fontFace: F.body, margin: 0 });

  s.addText([
    { text: "What epykit takes as input: ", options: { bold: true, color: C.navy } },
    { text: "Bismark `.cov` or MethylDackel `.bedGraph` files (chrom, pos, β, n_methylated, n_total) for each sample, plus a samplesheet labelling case / control.", options: { color: C.ink, breakLine: true } },
    { text: "What it produces: ", options: { bold: true, color: C.navy } },
    { text: "per-CpG DMC table (pvalue, qvalue, Δβ), per-region DMR table, annotation (gene / feature / CpG island context), QC report, HTML summary.", options: { color: C.ink } },
  ], { x: 0.6, y: 5.7, w: W - 1.2, h: 1.3, fontSize: 12, fontFace: F.body, margin: 0 });

  pageFooter(s, "Framing — the WGBS pipeline", 3);
  speakerNote(s,
    "Six-stage pipeline diagram. epykit handles stages 5 and 6 — the analytics, not the alignment. " +
    "Stages 1-4 are someone else's tools (Trim Galore, Bismark, MethylDackel). " +
    "Mention that this is a deliberate scope choice: nf-core/methylseq already does the upstream well."
  );
}

// =========================================================================
// 4 — Δβ + coverage primer
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(s,
    "The two numbers that drive DMC calling: Δβ and coverage",
    "Effect size and read depth jointly determine the per-CpG signal-to-noise"
  );

  // Cartoon CpG: panel of small dots, methylated vs unmethylated
  // Left: case sample, Right: control sample

  s.addShape(pres.shapes.RECTANGLE, { x: 0.6, y: 1.7, w: 5.6, h: 5.0,
    fill: { color: C.cardBg }, line: { color: C.rule, width: 1 } });
  s.addText("Per-CpG β (methylation fraction)", {
    x: 0.8, y: 1.85, w: 5.2, h: 0.4,
    fontSize: 13, fontFace: F.head, bold: true, color: C.navy, margin: 0,
  });

  s.addText([
    { text: "β = ", options: { bold: true, color: C.ink, fontFace: F.mono } },
    { text: "n_methylated_reads / n_total_reads ", options: { color: C.ink, fontFace: F.mono, breakLine: true } },
    { text: "", options: { breakLine: true } },
    { text: "A site with coverage = 4 and 2 methylated reads has β = 0.50, but the sampling variance is enormous (Bernoulli with n=4). ", options: { color: C.ink, breakLine: true } },
    { text: "The same true β = 0.50 with coverage = 40 and 20 methylated reads is ", options: { color: C.ink } },
    { text: "much more reliable.", options: { italic: true, color: C.navy } },
  ], { x: 0.8, y: 2.4, w: 5.2, h: 3.0, fontSize: 12, fontFace: F.body, margin: 0 });

  s.addText([
    { text: "Δβ = β_case − β_control", options: { bold: true, color: C.navy, fontFace: F.mono, breakLine: true } },
    { text: "is the effect-size statistic. Reviewers typically demand |Δβ| ≥ 0.10 with q ≤ 0.05 for biological relevance.", options: { color: C.ink } },
  ], { x: 0.8, y: 5.3, w: 5.2, h: 1.2, fontSize: 12, fontFace: F.body, margin: 0 });

  s.addShape(pres.shapes.RECTANGLE, { x: 6.5, y: 1.7, w: 6.2, h: 5.0,
    fill: { color: C.cardBg2 }, line: { color: C.rule, width: 1 } });
  s.addText("Why coverage matters statistically", {
    x: 6.7, y: 1.85, w: 5.8, h: 0.4,
    fontSize: 13, fontFace: F.head, bold: true, color: C.navy, margin: 0,
  });

  const rows = [
    { cov: "5×",  pwr: "low — most methods false-negative",   col: C.bad },
    { cov: "10×", pwr: "marginal — engine choice matters",     col: C.warn },
    { cov: "15×", pwr: "comfortable for most engines",         col: C.good },
    { cov: "25×", pwr: "near-saturation — engines converge",   col: C.good },
  ];
  rows.forEach((r, i) => {
    const y = 2.5 + i * 0.6;
    s.addShape(pres.shapes.RECTANGLE, { x: 6.7, y, w: 0.85, h: 0.5,
      fill: { color: r.col }, line: { type: "none" } });
    s.addText(r.cov, {
      x: 6.7, y, w: 0.85, h: 0.5,
      fontSize: 13, fontFace: F.head, bold: true, color: C.white, align: "center", margin: 0,
    });
    s.addText(r.pwr, {
      x: 7.7, y: y + 0.07, w: 4.8, h: 0.4,
      fontSize: 12, fontFace: F.body, color: C.ink, margin: 0,
    });
  });

  s.addText([
    { text: "The benchmark sweeps coverage from 5× to 25× ", options: { bold: true, color: C.navy } },
    { text: "specifically because the low-coverage regime is where statistical engines diverge — and where epykit's `lr` engine outperforms.", options: { color: C.ink } },
  ], { x: 6.7, y: 5.2, w: 5.8, h: 1.3, fontSize: 12, fontFace: F.body, italic: true, margin: 0 });

  pageFooter(s, "Framing — the per-CpG signal", 4);
  speakerNote(s,
    "Δβ and coverage are the two numbers people fight about in methylation papers. " +
    "Δβ is the biological effect size; coverage is the sampling precision. " +
    "The benchmark we'll see varies both — at 5× coverage everyone struggles, at 25× everyone wins. " +
    "The interesting story is the middle."
  );
}

// =========================================================================
// 5 — Tooling landscape
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(s,
    "Why a new tool? The methylation ecosystem is R-monoculture",
    "Python-native scientific computing has no comparable mature option"
  );

  // Two-column layout: existing tools vs the gap
  const tools = [
    { name: "methylKit",  lang: "R / Bioconductor", pros: "logistic-regression + Fisher; widely cited; tile DMRs", cons: "in-memory R only; slow at WGBS scale" },
    { name: "DSS",        lang: "R / Bioconductor", pros: "smoothed beta-binomial; gold standard for low-rep designs", cons: "BSmooth coupling; serial processing" },
    { name: "BSmooth",    lang: "R / Bioconductor", pros: "local smoothing; permutation-based",                     cons: "memory blows up past 5M CpGs" },
    { name: "RADMeth",    lang: "C++ CLI",          pros: "GLM with covariates; competitive accuracy",              cons: "Linux-only; spartan I/O" },
    { name: "methylSig",  lang: "R / Bioconductor", pros: "beta-binomial with covariates",                          cons: "limited DMR support" },
    { name: "BiSeq",      lang: "R / Bioconductor", pros: "smooth + segmented",                                     cons: "RRBS-targeted; not WGBS-native" },
  ];

  const tableX = 0.6;
  const tableY = 1.55;
  const colW = [1.6, 1.7, 4.1, 3.3];
  const rowH = 0.5;

  // Header
  const headers = ["Tool", "Language", "Strength", "Limitation"];
  headers.forEach((h, c) => {
    let cx = tableX;
    for (let k = 0; k < c; k++) cx += colW[k];
    s.addShape(pres.shapes.RECTANGLE, { x: cx, y: tableY, w: colW[c], h: rowH,
      fill: { color: C.navy }, line: { type: "none" } });
    s.addText(h, {
      x: cx + 0.1, y: tableY, w: colW[c] - 0.2, h: rowH,
      fontSize: 11, fontFace: F.head, bold: true, color: C.white, valign: "middle", margin: 0,
    });
  });

  tools.forEach((t, i) => {
    const y = tableY + (i + 1) * rowH;
    const bg = i % 2 === 0 ? C.cardBg : C.white;
    let cx = tableX;
    [t.name, t.lang, t.pros, t.cons].forEach((cell, c) => {
      s.addShape(pres.shapes.RECTANGLE, { x: cx, y, w: colW[c], h: rowH,
        fill: { color: bg }, line: { color: C.rule, width: 0.5 } });
      s.addText(cell, {
        x: cx + 0.1, y, w: colW[c] - 0.2, h: rowH,
        fontSize: 10, fontFace: c === 0 ? F.head : F.body,
        bold: c === 0, color: C.ink, valign: "middle", margin: 0,
      });
      cx += colW[c];
    });
  });

  s.addShape(pres.shapes.RECTANGLE, {
    x: tableX, y: tableY + (tools.length + 1) * rowH + 0.2, w: W - 1.2, h: 1.55,
    fill: { color: C.navy }, line: { type: "none" },
  });
  s.addText("The gap: a Python-native, scanpy-style API for methylation", {
    x: tableX + 0.2, y: tableY + (tools.length + 1) * rowH + 0.3, w: 12, h: 0.4,
    fontSize: 14, fontFace: F.head, bold: true, color: C.gold, margin: 0,
  });
  s.addText([
    { text: "Single-cell methylation, multi-omics integration, and modern ML workflows all live in Python — `scanpy`, `anndata`, `mudata`. ", options: { color: "DCE7EF" } },
    { text: "Bridging into R for every methylation step kills reproducibility. ", options: { italic: true, color: "FFE893" } },
    { text: "epykit fills that gap.", options: { bold: true, color: C.gold } },
  ], { x: tableX + 0.2, y: tableY + (tools.length + 1) * rowH + 0.7, w: 12.6, h: 0.9,
       fontSize: 12, fontFace: F.body, margin: 0 });

  pageFooter(s, "Framing — the tool landscape", 5);
  speakerNote(s,
    "Don't bash R. R/Bioconductor methylation tools are excellent — methylKit is the gold standard. " +
    "The gap is one of language ecosystem, not algorithm quality. " +
    "Python single-cell tooling is now where the field is moving; methylation tooling should follow."
  );
}

// =========================================================================
// 6 — epykit positioning
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(s,
    "epykit's positioning",
    "What does it bring that the others don't?"
  );

  const features = [
    { icon: "1", title: "Parquet-backed storage", desc: "Per-chromosome, per-sample partitioned Parquet via Polars lazy scans. 22 M CpGs never load as a single in-memory frame." },
    { icon: "2", title: "Streaming DMC engines",  desc: "Per-chromosome processing with a manifest-based DMCStore; peak memory is O(largest chrom), not O(genome)." },
    { icon: "3", title: "scanpy-style API",        desc: "ep.pp.* for preprocessing, ep.tl.* for tools, ep.pl.* for plots. Mirrors the conventions Python bioinformaticians already know." },
    { icon: "4", title: "Five statistical engines", desc: "lr (quasi-binomial LR, default), lr+ (4-knob power stack), glm (IRLS with covariates), welch_t, fisher. Plus four DMR callers." },
    { icon: "5", title: "First-class interop",     desc: "AnnData / MuData export. MethylKit and MultiQC integration. nf-core/methylseq sample-sheet support out of the box." },
    { icon: "6", title: "Windows + Linux + macOS", desc: "Pure-Python core with optional C extensions; CI matrix covers Windows-py3.9 / Windows-py3.12 / Linux-py3.9 / Linux-py3.12." },
  ];

  const gridX0 = 0.6;
  const gridY0 = 1.55;
  const cardW = 4.0;
  const cardH = 2.55;
  const colGap = 0.10;
  const rowGap = 0.15;

  features.forEach((f, i) => {
    const col = i % 3;
    const row = Math.floor(i / 3);
    const x = gridX0 + col * (cardW + colGap);
    const y = gridY0 + row * (cardH + rowGap);

    s.addShape(pres.shapes.RECTANGLE, { x, y, w: cardW, h: cardH,
      fill: { color: C.cardBg }, line: { color: C.rule, width: 1 } });
    s.addShape(pres.shapes.RECTANGLE, { x, y, w: 0.10, h: cardH,
      fill: { color: C.gold }, line: { type: "none" } });
    s.addText(f.icon, {
      x: x + 0.25, y: y + 0.15, w: 0.5, h: 0.5,
      fontSize: 26, fontFace: F.head, bold: true, color: C.gold, margin: 0,
    });
    s.addText(f.title, {
      x: x + 0.85, y: y + 0.2, w: cardW - 1.0, h: 0.5,
      fontSize: 13, fontFace: F.head, bold: true, color: C.navy, margin: 0,
    });
    s.addText(f.desc, {
      x: x + 0.25, y: y + 0.85, w: cardW - 0.4, h: cardH - 1.0,
      fontSize: 10.5, fontFace: F.body, color: C.ink, margin: 0,
    });
  });

  pageFooter(s, "Framing — what epykit brings", 6);
  speakerNote(s,
    "Six bullets. The first two (Parquet storage + streaming engines) are the engineering. " +
    "The third is the API design choice. Four is the statistics surface. Five and six are user-facing concerns. " +
    "Skip details — say 'we'll see how each of these pays off in the benchmark'."
  );
}

// =========================================================================
// 7 — Thesis contribution
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.navy };
  titleHeader(s,
    "Thesis contribution in one sentence",
    null
  );
  // Override the title color since the bg is navy
  s.addText("Thesis contribution in one sentence", {
    x: 0.6, y: 0.32, w: W - 1.2, h: 0.7,
    fontSize: 28, fontFace: F.head, bold: true, color: C.white, margin: 0,
  });

  s.addText([
    { text: "epykit", options: { fontSize: 64, fontFace: F.head, bold: true, color: C.gold } },
    { text: " is the first Python-native WGBS differential-methylation pipeline whose accuracy matches the R/Bioconductor gold standards while running ", options: { fontSize: 30, fontFace: F.body, color: C.white } },
    { text: "12× faster", options: { fontSize: 36, fontFace: F.head, bold: true, color: C.gold } },
    { text: " on real data and integrating cleanly with the modern Python single-cell / multi-omics stack.", options: { fontSize: 30, fontFace: F.body, color: C.white } },
  ], { x: 1.5, y: 2.0, w: 10.3, h: 4.5, paraSpaceAfter: 0, valign: "top", margin: 0 });

  s.addText("Validated on three studies: 8-tool panel comparison on simulated data, " +
            "head-to-head against methylKit on the same simulator with fresh local runs, " +
            "and a real WGBS reproduction of Farhangdoost et al. 2024 (GSE263850).",
    { x: 1.5, y: 5.9, w: 10.3, h: 1.0,
      fontSize: 14, fontFace: F.body, italic: true, color: "BFD7E5", margin: 0 });

  pageFooter(s, "Framing — the one-line claim", 7);
  speakerNote(s,
    "Pause on this slide. This is the thesis claim — read it slowly. " +
    "Everything after is evidence for this sentence."
  );
}

// =========================================================================
// 8 — Architecture: methylstore
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(s,
    "Architecture: the methylstore is the source of truth",
    "Partitioned Parquet under .cache/, repointed after every pipeline step"
  );

  // Pipeline diagram showing data flow
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 0.6, y: 1.6, w: 2.3, h: 1.2,
    fill: { color: C.teal }, line: { type: "none" }, rectRadius: 0.10 });
  s.addText("Bismark .cov\nMethylDackel .bg\nnf-core/methylseq", {
    x: 0.6, y: 1.6, w: 2.3, h: 1.2,
    fontSize: 11, fontFace: F.body, color: C.white, align: "center", valign: "middle", bold: true, margin: 0,
  });

  s.addText("read_bismark()", {
    x: 3.0, y: 1.75, w: 1.4, h: 0.4,
    fontSize: 11, fontFace: F.mono, italic: true, color: C.navy, align: "center", margin: 0,
  });
  s.addText("→", {
    x: 3.0, y: 2.05, w: 1.4, h: 0.4,
    fontSize: 22, fontFace: F.head, color: C.navy, align: "center", margin: 0,
  });

  // Methylstore box
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 4.5, y: 1.4, w: 5.6, h: 1.6,
    fill: { color: C.navy }, line: { color: C.gold, width: 2 }, rectRadius: 0.10 });
  s.addText("MethylStore  (Parquet)", {
    x: 4.5, y: 1.45, w: 5.6, h: 0.4,
    fontSize: 14, fontFace: F.head, bold: true, color: C.gold, align: "center", margin: 0,
  });
  s.addText(
    "<store>/.cache/raw/sample=<id>/chrom=<chr>/part-0.parquet\n" +
    "Polars lazy scans  •  zero whole-genome materialisation",
    {
      x: 4.5, y: 1.95, w: 5.6, h: 1.0,
      fontSize: 10, fontFace: F.mono, color: "BFD7E5", align: "center", margin: 0,
    });

  // Step boxes - left-to-right pipeline
  s.addText("Each pp.* / tl.* step writes a new cache directory and repoints md.store", {
    x: 0.6, y: 3.3, w: 12, h: 0.4,
    fontSize: 12, fontFace: F.body, italic: true, color: C.muted, align: "center", margin: 0,
  });

  const steps = [
    { name: "filter",    desc: ".cache/filtered" },
    { name: "normalize", desc: ".cache/normalized" },
    { name: "unite",     desc: ".cache/united" },
    { name: "tl.dmc",    desc: ".cache/dmc" },
    { name: "tl.dmr",    desc: "md.uns['dmr']" },
    { name: "annotate",  desc: "md.varm[\"annot\"]" },
  ];
  const stepW = 1.95;
  const stepGap = 0.10;
  const stepX0 = 0.6;
  const stepY = 3.85;
  steps.forEach((st, i) => {
    const x = stepX0 + i * (stepW + stepGap);
    s.addShape(pres.shapes.RECTANGLE, { x, y: stepY, w: stepW, h: 1.5,
      fill: { color: i < 3 ? "DBEAFE" : "BFDBFE" }, line: { color: C.blue, width: 1 } });
    s.addText(st.name, {
      x, y: stepY + 0.15, w: stepW, h: 0.5,
      fontSize: 14, fontFace: F.mono, bold: true, color: C.navy, align: "center", margin: 0,
    });
    s.addText(st.desc, {
      x, y: stepY + 0.75, w: stepW, h: 0.6,
      fontSize: 9, fontFace: F.mono, color: C.ink, align: "center", margin: 0,
    });
  });

  // Why it matters
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.6, y: 5.7, w: W - 1.2, h: 1.4,
    fill: { color: C.cardBg }, line: { color: C.rule, width: 1 } });
  s.addText("Why this layout matters", {
    x: 0.8, y: 5.78, w: 12, h: 0.4,
    fontSize: 13, fontFace: F.head, bold: true, color: C.navy, margin: 0,
  });
  s.addText([
    { text: "•  Peak RAM is O(largest chromosome), not O(genome) — a 22 M-CpG dataset never lives in a single frame.", options: { color: C.ink, breakLine: true } },
    { text: "•  Every step is auditable: ", options: { color: C.ink } },
    { text: "uns[\"_store_history\"]", options: { fontFace: F.mono, color: C.navy } },
    { text: " is the source of truth for what ran; preprocessing flags are derived from it, not stored independently.", options: { color: C.ink, breakLine: true } },
    { text: "•  Resumability — re-running the same pipeline reuses cached parquets if inputs and parameters haven't changed.", options: { color: C.ink } },
  ], { x: 0.8, y: 6.15, w: 12.0, h: 1.0, fontSize: 11, fontFace: F.body, margin: 0 });

  pageFooter(s, "Architecture", 8);
  speakerNote(s,
    "The data layout is the central engineering decision. Each step writes a new cache and updates md.store. " +
    "This means we never materialise the genome in RAM, the steps are resumable, " +
    "and the _store_history field gives us a tamper-resistant audit trail."
  );
}

// =========================================================================
// 9 — Engines + lr+ stack
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(s,
    "Five statistical engines + four DMR callers + the lr+ power stack",
    "All implemented in vectorised Python (NumPy / SciPy / Numba)"
  );

  s.addShape(pres.shapes.RECTANGLE, { x: 0.6, y: 1.55, w: 6.05, h: 5.4,
    fill: { color: C.cardBg }, line: { color: C.rule, width: 1 } });
  s.addText("DMC engines (per-CpG)", {
    x: 0.8, y: 1.65, w: 5.8, h: 0.4,
    fontSize: 14, fontFace: F.head, bold: true, color: C.navy, margin: 0,
  });

  const engines = [
    { name: "lr",       desc: "Quasi-binomial likelihood-ratio (default, n ≥ 2)" },
    { name: "lr+",      desc: "lr + 4-knob power stack (opt-in)" },
    { name: "glm",      desc: "IRLS binomial GLM with covariates (~ group + batch)" },
    { name: "welch_t",  desc: "Welch t-test on raw β" },
    { name: "fisher",   desc: "Pooled Fisher exact (n = 1 fallback)" },
  ];
  engines.forEach((e, i) => {
    const y = 2.15 + i * 0.55;
    s.addShape(pres.shapes.RECTANGLE, { x: 0.8, y, w: 1.0, h: 0.4,
      fill: { color: C.navy }, line: { type: "none" } });
    s.addText(e.name, {
      x: 0.8, y, w: 1.0, h: 0.4,
      fontSize: 12, fontFace: F.mono, bold: true, color: C.gold, align: "center", margin: 0,
    });
    s.addText(e.desc, {
      x: 1.95, y, w: 4.6, h: 0.4,
      fontSize: 11, fontFace: F.body, color: C.ink, valign: "middle", margin: 0,
    });
  });

  s.addText("DMR callers (per-region)", {
    x: 0.8, y: 5.0, w: 5.8, h: 0.4,
    fontSize: 14, fontFace: F.head, bold: true, color: C.navy, margin: 0,
  });
  const callers = [
    { name: "chain_merge",     desc: "DSS-compatible (default), three presets" },
    { name: "tile",            desc: "Read-pooled fixed-window tiles" },
    { name: "sliding_window",  desc: "Signed Stouffer combiner across CpG neighbourhoods" },
    { name: "segment",         desc: "HMM-based segmentation" },
  ];
  callers.forEach((e, i) => {
    const y = 5.5 + i * 0.35;
    s.addShape(pres.shapes.RECTANGLE, { x: 0.8, y, w: 1.65, h: 0.3,
      fill: { color: C.teal }, line: { type: "none" } });
    s.addText(e.name, {
      x: 0.8, y, w: 1.65, h: 0.3,
      fontSize: 11, fontFace: F.mono, bold: true, color: C.white, align: "center", margin: 0,
    });
    s.addText(e.desc, {
      x: 2.6, y, w: 4.0, h: 0.3,
      fontSize: 10, fontFace: F.body, color: C.ink, valign: "middle", margin: 0,
    });
  });

  // lr+ power stack
  s.addShape(pres.shapes.RECTANGLE, { x: 6.85, y: 1.55, w: 5.85, h: 5.4,
    fill: { color: C.navy }, line: { color: C.gold, width: 2 } });
  s.addText("lr+ power stack — four opt-in knobs", {
    x: 7.05, y: 1.65, w: 5.5, h: 0.4,
    fontSize: 14, fontFace: F.head, bold: true, color: C.gold, margin: 0,
  });

  const knobs = [
    { num: "1", name: "Empirical-Bayes dispersion shrinkage", desc: "Inverse-Gamma posterior on per-CpG φ; lifts small-n TPR" },
    { num: "2", name: "Sign-aware Stouffer neighbour combiner", desc: "Combines same-direction neighbours ≤ 200 bp" },
    { num: "3", name: "Separation-aware Fisher fallback", desc: "Re-test sites with |Δβ| ≥ 0.9 via Fisher; take min(p)" },
    { num: "4", name: "Storey two-stage Benjamini-Hochberg", desc: "Data-driven π₀ estimate (vs BH's π₀ = 1)" },
  ];
  knobs.forEach((k, i) => {
    const y = 2.2 + i * 1.05;
    s.addText(k.num, {
      x: 7.05, y, w: 0.5, h: 0.7,
      fontSize: 28, fontFace: F.head, bold: true, color: C.gold, margin: 0,
    });
    s.addText(k.name, {
      x: 7.55, y, w: 5.1, h: 0.4,
      fontSize: 12, fontFace: F.body, bold: true, color: C.white, margin: 0,
    });
    s.addText(k.desc, {
      x: 7.55, y: y + 0.4, w: 5.1, h: 0.6,
      fontSize: 10, fontFace: F.body, italic: true, color: "BFD7E5", margin: 0,
    });
  });

  s.addText("Activated together via `power_stack=\"lr+\"`. Per-knob attribution is not isolatable from current eval data — they're presented as a bundle. Trade-off in F10.", {
    x: 7.05, y: 6.45, w: 5.5, h: 0.45,
    fontSize: 9, fontFace: F.body, italic: true, color: "8FB6CA", margin: 0,
  });

  pageFooter(s, "Architecture", 9);
  speakerNote(s,
    "Five engines + four DMR callers + the four-knob lr+ stack. " +
    "Don't try to memorise all of these — the headline is lr (default) and lr+ (opt-in extra power). " +
    "Mention that all four lr+ knobs activate together; we cannot isolate per-knob impact from the data."
  );
}

// =========================================================================
// 10 — API example
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(s,
    "The API surface, end to end",
    "Five lines from raw `.cov` files to annotated DMRs"
  );

  s.addShape(pres.shapes.RECTANGLE, { x: 0.6, y: 1.5, w: W - 1.2, h: 5.2,
    fill: { color: "0F172A" }, line: { type: "none" } });

  const code = [
    { line: "import epykit as ep", comment: "" },
    { line: "", comment: "" },
    { line: "md = ep.io.read_bismark(", comment: "  # 1. ingest .cov files into Parquet methylstore" },
    { line: "    samplesheet='samples.csv', store_dir='store/', assembly='hg38')", comment: "" },
    { line: "", comment: "" },
    { line: "ep.pp.filter(md, min_coverage=5, max_coverage=500)", comment: "  # 2. coverage filter" },
    { line: "ep.pp.unite(md, mode='intersect')", comment: "  # 3. intersect CpGs across samples" },
    { line: "", comment: "" },
    { line: "ep.tl.dmc(md, engine='lr', formula='~ group')", comment: "  # 4. per-CpG differential test" },
    { line: "ep.tl.dmr(md, engine='chain_merge', preset='default')", comment: "  # 5. region calling" },
    { line: "", comment: "" },
    { line: "ep.tl.annotate(md, gtf='refGene.gtf.gz')", comment: "  # gene + CpG-island context" },
    { line: "ep.report.html(md, out='report.html')", comment: "  # self-contained HTML report" },
  ];

  let y = 1.75;
  code.forEach((l) => {
    if (l.line.trim() === "" && l.comment === "") {
      y += 0.25;
      return;
    }
    s.addText([
      { text: l.line, options: { color: l.line.startsWith("#") ? "#7F8C8D" : "E2E8F0", fontFace: F.mono, fontSize: 14 } },
      { text: l.comment, options: { color: "#94A3B8", fontFace: F.mono, fontSize: 11, italic: true } },
    ], { x: 1.0, y, w: W - 1.5, h: 0.35, margin: 0 });
    y += 0.35;
  });

  pageFooter(s, "Architecture — API", 10);
  speakerNote(s,
    "Walk through line by line. Read the input files, filter to good-coverage sites, intersect, " +
    "run the DMC engine, call regions, annotate, write a report. " +
    "Mention this mirrors scanpy's pp/tl/pl convention — Python bioinformaticians can pattern-match instantly."
  );
}

// =========================================================================
// 11 — Study 1 design
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(s,
    "Study 1 — panel comparison on the Piao 2021 simulator",
    "Reproducing an established benchmark to position epykit among eight peer tools"
  );

  // Left: study design card
  s.addShape(pres.shapes.RECTANGLE, { x: 0.6, y: 1.55, w: 6.05, h: 5.4,
    fill: { color: C.cardBg }, line: { color: C.rule, width: 1 } });
  s.addText("Design", {
    x: 0.8, y: 1.65, w: 5.8, h: 0.4,
    fontSize: 14, fontFace: F.head, bold: true, color: C.navy, margin: 0,
  });
  s.addText([
    { text: "Simulator:", options: { bold: true, color: C.navy } },
    { text: " Piao et al. 2021 ", options: { color: C.ink } },
    { text: "(NAR Genom Bioinform; 100 k CpGs / scenario; 20 k spike-in true DMCs)", options: { color: C.muted, breakLine: true } },
    { text: "Coverage sweep:", options: { bold: true, color: C.navy } },
    { text: " 5, 10, 15, 20, 25 × at fixed n = 3 vs 3", options: { color: C.ink, breakLine: true } },
    { text: "Replicate sweep:", options: { bold: true, color: C.navy } },
    { text: " n_total = 2, 4, 6, 8, 10 at fixed cov = 10×", options: { color: C.ink, breakLine: true } },
    { text: "Effect-size bins:", options: { bold: true, color: C.navy } },
    { text: " 0.2–0.4 / 0.4–0.6 / 0.6–0.8 / 0.8–1.0", options: { color: C.ink, breakLine: true } },
    { text: "Metrics:", options: { bold: true, color: C.navy } },
    { text: " TPR, FPR, F1, AUROC at q ≤ 0.05", options: { color: C.ink, breakLine: true } },
    { text: "", options: { breakLine: true } },
    { text: "Tools:", options: { bold: true, color: C.navy, breakLine: true } },
    { text: "  • epykit lr / lr+ / welch_t / fisher", options: { fontFace: F.mono, color: C.ink, breakLine: true } },
    { text: "  • methylKit (plain + tuned)", options: { fontFace: F.mono, color: C.ink, breakLine: true } },
    { text: "  • DSS, RADMeth, BiSeq, methylSig, Fisher (Piao 2021 baselines)", options: { fontFace: F.mono, color: C.ink } },
  ], { x: 0.8, y: 2.1, w: 5.8, h: 5.0, fontSize: 12, fontFace: F.body, paraSpaceAfter: 4, margin: 0 });

  // Right: methodology card
  s.addShape(pres.shapes.RECTANGLE, { x: 6.85, y: 1.55, w: 5.85, h: 5.4,
    fill: { color: C.navy }, line: { type: "none" } });
  s.addText("Why this matters", {
    x: 7.05, y: 1.65, w: 5.5, h: 0.4,
    fontSize: 14, fontFace: F.head, bold: true, color: C.gold, margin: 0,
  });
  s.addText([
    { text: "•  Piao et al. 2021 already evaluated DSS, RADMeth, BiSeq, methylSig and Fisher on this simulator. We transcribe their tables (Supp Tables S1 & S2, ", options: { color: "DCE7EF" } },
    { text: "audited cell-by-cell", options: { italic: true, color: C.gold } },
    { text: ", 0 transcription errors).", options: { color: "DCE7EF", breakLine: true } },
    { text: "", options: { breakLine: true } },
    { text: "•  We add fresh measurements for the four epykit engines, plus a fresh local run of methylKit and a parameter-tuned methylKit (`methylkit_tuned`).", options: { color: "DCE7EF", breakLine: true } },
    { text: "", options: { breakLine: true } },
    { text: "•  This gives ", options: { color: "DCE7EF" } },
    { text: "apples-to-apples cross-tool comparisons ", options: { bold: true, color: C.white } },
    { text: "on identical simulated input.", options: { color: "DCE7EF" } },
  ], { x: 7.05, y: 2.1, w: 5.5, h: 4.5, fontSize: 12, fontFace: F.body, margin: 0 });

  pageFooter(s, "Study 1 — simulator panel", 11);
  speakerNote(s,
    "Study 1 anchors us to an established benchmark. " +
    "Piao 2021 published a panel comparison; we re-run it with epykit added, " +
    "and we cross-verify by also re-running methylKit on the same input."
  );
}

// =========================================================================
// 12 — F01: Cross-tool TPR/FPR
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(s,
    "Cross-tool TPR / FPR vs coverage — epykit `lr` matches the strongest baselines",
    "Smallest effect-size bin (0.2 – 0.4), the regime that discriminates between methods"
  );

  s.addImage({
    path: HEAD("F01_tool_panel_TPR_FPR.png"),
    x: 0.5, y: 1.45, w: 12.3, h: 5.2, sizing: { type: "contain", w: 12.3, h: 5.2 },
  });

  s.addText([
    { text: "Take-away:  ", options: { bold: true, color: C.navy } },
    { text: "epykit `lr` and `lr+` (blues) sit at the top of the TPR curve from 5× onward and stay below 0.005 FPR. methylKit plain (red, n=3) is a peer competitor. RADMeth and DSS catch up by 10× coverage.", options: { color: C.ink } },
  ], { x: 0.5, y: 6.65, w: 12.3, h: 0.5, fontSize: 11, fontFace: F.body, margin: 0 });

  pageFooter(s, "Study 1 — accuracy headline", 12);
  speakerNote(s,
    "Read off the figure: at coverage 5×, epykit lr is already at TPR ≈ 0.83. " +
    "methylKit (plain n=3) is at ≈ 0.27. DSS at ≈ 0.07. " +
    "By 10× everything has converged. The FPR row shows nobody is calling false positives at q=0.05. " +
    "F1 row only has epykit engines — baselines from Piao didn't publish F1."
  );
}

// =========================================================================
// 13 — F08: Pareto frontier
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(s,
    "Speed × accuracy Pareto frontier at the hardest scenario",
    "5× coverage, n = 3 vs 3 — where method choice matters most"
  );

  s.addImage({
    path: HEAD("F08_pareto_runtime_vs_f1.png"),
    x: 1.2, y: 1.45, w: 11.0, h: 5.2, sizing: { type: "contain", w: 11.0, h: 5.2 },
  });

  s.addText([
    { text: "Pareto frontier:  ", options: { bold: true, color: C.navy } },
    { text: "epykit `lr` is the only point in the fast-and-accurate corner (F1 > 0.85, runtime < 1 minute). methylKit_tuned reaches higher F1 but takes ≈ 2 min per scenario. RADMeth is accurate but slow. DSS is fast but its F1 lags.", options: { color: C.ink } },
  ], { x: 0.5, y: 6.7, w: 12.3, h: 0.45, fontSize: 11, fontFace: F.body, margin: 0 });

  pageFooter(s, "Study 1 — Pareto frontier", 13);
  speakerNote(s,
    "If there's one figure to remember, it's this. " +
    "Top-left = fast and accurate; bottom-right = slow and inaccurate. " +
    "epykit lr is alone in the top-left at 0.91 F1 / 0.01 minutes. " +
    "methylkit_tuned beats it on F1 but is 100× slower."
  );
}

// =========================================================================
// 14 — F09: Replicate scaling
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(s,
    "Replicate scaling — what happens at n = 2?",
    "The regime where overdispersion estimates break"
  );

  s.addImage({
    path: HEAD("F09_n2_power_gap.png"),
    x: 0.5, y: 1.45, w: 12.3, h: 5.0, sizing: { type: "contain", w: 12.3, h: 5.0 },
  });

  s.addText([
    { text: "Honest read-out:  ", options: { bold: true, color: C.navy } },
    { text: "at n = 2, only methylKit_tuned and epykit's Fisher fallback recover non-zero TPR. epykit's lr / lr+ engines require n ≥ 4 per group. All tools converge by n ≥ 6.", options: { color: C.ink } },
  ], { x: 0.5, y: 6.55, w: 12.3, h: 0.55, fontSize: 11, fontFace: F.body, margin: 0 });

  pageFooter(s, "Study 1 — replicate scaling", 14);
  speakerNote(s,
    "This is the honesty slide. n=2 is methylKit_tuned's domain (TPR = 0.96). " +
    "epykit at n=2 has to fall back to the Fisher exact engine (TPR = 0.29). " +
    "From n ≥ 4 onward, the tools are indistinguishable on accuracy — at which point speed and ergonomics decide. " +
    "Mention: lr / lr+ requires n ≥ 4 because the overdispersion estimator needs degrees of freedom."
  );
}

// =========================================================================
// 15 — Study 2 design
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(s,
    "Study 2 — fresh head-to-head against methylKit on the same simulator",
    "Removing the variable: identical input, identical machine, same wall clock"
  );

  s.addText([
    { text: "Motivation: ", options: { bold: true, color: C.navy } },
    { text: "the Piao 2021 baselines are accuracy-only — their tables don't include runtime or memory. To compare ", options: { color: C.ink } },
    { text: "speed", options: { italic: true, color: C.navy } },
    { text: " honestly, we need fresh local runs.", options: { color: C.ink, breakLine: true } },
    { text: "", options: { breakLine: true } },
    { text: "Design:", options: { bold: true, color: C.navy, breakLine: true } },
    { text: "•  Same simulator scenarios as Study 1 (coverage 5–25×, n_total 2–10).", options: { color: C.ink, breakLine: true } },
    { text: "•  Same OS-level resource tracker (psutil snapshots every second).", options: { color: C.ink, breakLine: true } },
    { text: "•  Same machine: 16-thread Windows desktop, 32 GB RAM.", options: { color: C.ink, breakLine: true } },
    { text: "•  Two methylKit configurations: ", options: { color: C.ink } },
    { text: "plain ", options: { italic: true, color: C.navy } },
    { text: "(default parameters) and ", options: { color: C.ink } },
    { text: "tuned ", options: { italic: true, color: C.navy } },
    { text: "(swept overdispersion + q-cutoff).", options: { color: C.ink } },
  ], { x: 0.6, y: 1.55, w: W - 1.2, h: 3.0, fontSize: 13, fontFace: F.body, paraSpaceAfter: 4, margin: 0 });

  s.addShape(pres.shapes.RECTANGLE, { x: 0.6, y: 4.6, w: W - 1.2, h: 2.3,
    fill: { color: C.cardBg2 }, line: { color: C.rule, width: 1 } });
  s.addText("Why this is not just \"more numbers\"", {
    x: 0.8, y: 4.7, w: 12, h: 0.4,
    fontSize: 13, fontFace: F.head, bold: true, color: C.navy, margin: 0,
  });
  s.addText([
    { text: "Most benchmark papers compare wall-clock numbers transcribed from each tool's own publication, where each was run on different hardware with different I/O assumptions. Those numbers are mostly noise. Study 2 strips that out — ", options: { color: C.ink } },
    { text: "every wall-clock number on the next slide comes from the same machine on the same day with the same input.", options: { bold: true, color: C.navy } },
  ], { x: 0.8, y: 5.1, w: 12.0, h: 1.7, fontSize: 12, fontFace: F.body, margin: 0 });

  pageFooter(s, "Study 2 — head-to-head", 15);
  speakerNote(s,
    "Study 2 exists because runtime claims from R-tool papers are often cherry-picked. " +
    "By running on the same machine we get clean comparisons. " +
    "And we run methylKit twice (plain and tuned) because methylkit-tuned beats plain on Piao's grid."
  );
}

// =========================================================================
// 16 — F08 reused for sim speed (or use Study 2 figure)
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(s,
    "Study 2 — head-to-head runtime on the simulator",
    "epykit engines vs methylKit, identical input, identical hardware"
  );

  // Reuse the existing study2 figure if it has the right info
  const study2Fig = path.join(ROOT, "figures", "study2_simulated_headToHead", "F6_runtime.png");
  s.addImage({
    path: study2Fig,
    x: 1.2, y: 1.45, w: 11.0, h: 5.2, sizing: { type: "contain", w: 11.0, h: 5.2 },
  });

  s.addText([
    { text: "On 100 k-CpG simulator scenarios:  ", options: { bold: true, color: C.navy } },
    { text: "epykit lr completes a scenario in ≈ 0.9 seconds; methylKit takes 60-120 seconds. lr+ is slower than lr (≈ 7 s) because of the four-knob overhead, but still well under methylKit.", options: { color: C.ink } },
  ], { x: 0.5, y: 6.7, w: 12.3, h: 0.5, fontSize: 11, fontFace: F.body, margin: 0 });

  pageFooter(s, "Study 2 — runtime", 16);
  speakerNote(s,
    "On the simulator, epykit is 1-2 orders of magnitude faster per scenario. " +
    "Note lr+ overhead — the four knobs add about 7× wall time vs bare lr but still beat methylKit handily."
  );
}

// =========================================================================
// 17 — Study 3 design
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(s,
    "Study 3 — real WGBS on GSE263850",
    "Reproducing Farhangdoost et al. 2024 (AKAP11 heterozygous LOF in iPSC-derived cortical neurons)"
  );

  s.addShape(pres.shapes.RECTANGLE, { x: 0.6, y: 1.55, w: 6.05, h: 5.4,
    fill: { color: C.cardBg }, line: { color: C.rule, width: 1 } });
  s.addText("Dataset", {
    x: 0.8, y: 1.65, w: 5.8, h: 0.4,
    fontSize: 14, fontFace: F.head, bold: true, color: C.navy, margin: 0,
  });
  s.addText([
    { text: "GEO:  ", options: { bold: true, color: C.navy, fontFace: F.mono } },
    { text: "GSE263850 (6 samples)", options: { color: C.ink, fontFace: F.mono, breakLine: true } },
    { text: "Design:  ", options: { bold: true, color: C.navy } },
    { text: "3 Het-AKAP11-KO × 3 WT iPSC-derived cortical neurons", options: { color: C.ink, breakLine: true } },
    { text: "Assembly:  ", options: { bold: true, color: C.navy } },
    { text: "hg38", options: { color: C.ink, fontFace: F.mono, breakLine: true } },
    { text: "CpGs after `unite`:  ", options: { bold: true, color: C.navy } },
    { text: "15,597,046", options: { color: C.ink, fontFace: F.mono, breakLine: true } },
    { text: "", options: { breakLine: true } },
    { text: "Original paper:  ", options: { bold: true, color: C.navy } },
    { text: "Farhangdoost et al. 2024", options: { color: C.ink, breakLine: true } },
    { text: "  •  Used DSS for DMC + callDMR ", options: { italic: true, color: C.muted, breakLine: true } },
    { text: "  •  Published 813 DMRs in Supp Table 5", options: { italic: true, color: C.muted } },
  ], { x: 0.8, y: 2.1, w: 5.8, h: 5.0, fontSize: 12, fontFace: F.body, paraSpaceAfter: 3, margin: 0 });

  s.addShape(pres.shapes.RECTANGLE, { x: 6.85, y: 1.55, w: 5.85, h: 5.4,
    fill: { color: C.navy }, line: { type: "none" } });
  s.addText("Three-way comparison", {
    x: 7.05, y: 1.65, w: 5.5, h: 0.4,
    fontSize: 14, fontFace: F.head, bold: true, color: C.gold, margin: 0,
  });
  s.addText([
    { text: "We run three tools end-to-end on the same input:", options: { color: "DCE7EF", breakLine: true } },
    { text: "", options: { breakLine: true } },
    { text: "1.  ", options: { bold: true, color: C.gold } },
    { text: "epykit 1.0", options: { bold: true, color: C.white } },
    { text: "  (lr + chain_merge)", options: { color: "DCE7EF", breakLine: true } },
    { text: "2.  ", options: { bold: true, color: C.gold } },
    { text: "methylKit 1.36.0", options: { bold: true, color: C.white } },
    { text: "  (default pipeline)", options: { color: "DCE7EF", breakLine: true } },
    { text: "3.  ", options: { bold: true, color: C.gold } },
    { text: "DSS (R/Bioconductor)", options: { bold: true, color: C.white } },
    { text: "  (smoothing + callDMR)", options: { color: "DCE7EF", breakLine: true } },
    { text: "", options: { breakLine: true } },
    { text: "Plus the paper's own 813 published DMRs as a ground-truth comparator.", options: { italic: true, color: "BFD7E5", breakLine: true } },
    { text: "", options: { breakLine: true } },
    { text: "Hardware caveat:", options: { bold: true, color: "FFE893" } },
    { text: " methylKit ran on a 24-CPU / 64 GB Linux workstation; epykit and DSS ran on a 16-thread Windows desktop.", options: { color: "DCE7EF" } },
  ], { x: 7.05, y: 2.1, w: 5.5, h: 4.7, fontSize: 12, fontFace: F.body, paraSpaceAfter: 4, margin: 0 });

  pageFooter(s, "Study 3 — real WGBS", 17);
  speakerNote(s,
    "Real data. AKAP11 heterozygous LOF in iPSC-cortical neurons, n=6. " +
    "Three-way: our tool, methylKit, DSS, plus the published call set. " +
    "Hardware caveat is important to flag — methylKit ran on a much beefier machine. " +
    "The speed gap is therefore conservative."
  );
}

// =========================================================================
// 18 — F02: Real-data speed/memory
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(s,
    "Real-data speed and memory — the headline number",
    "epykit completes the full pipeline in 18 min using 12.3 GB peak RAM"
  );

  s.addImage({
    path: HEAD("F02_speed_memory_real_data.png"),
    x: 0.5, y: 1.45, w: 12.3, h: 5.2, sizing: { type: "contain", w: 12.3, h: 5.2 },
  });

  s.addText([
    { text: "Headline:  ", options: { bold: true, color: C.navy } },
    { text: "epykit is 12.15× faster than methylKit and 3.82× lower peak RAM. DSS is fast (40 min) and slightly lower memory than epykit (9.1 GB) but produces only 922 DMRs vs epykit's 3,433 lenient / 257 strict.", options: { color: C.ink } },
  ], { x: 0.5, y: 6.7, w: 12.3, h: 0.5, fontSize: 11, fontFace: F.body, margin: 0 });

  pageFooter(s, "Study 3 — headline", 18);
  speakerNote(s,
    "The headline number. 12× faster wall-clock, 3.8× less memory on a 15.6 M-CpG real dataset. " +
    "DSS is the closest competitor on memory — slightly better than us. " +
    "Mention the hardware caveat: methylKit ran on a 24-CPU box, so the wall-clock comparison is conservative."
  );
}

// =========================================================================
// 19 — F04: per-CpG + DMR concordance
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(s,
    "Per-CpG and DMR-level agreement with established tools",
    "When epykit calls a site, do the other tools agree?"
  );

  s.addImage({
    path: HEAD("F04_per_cpg_and_dmr_concordance.png"),
    x: 0.5, y: 1.40, w: 12.3, h: 5.3, sizing: { type: "contain", w: 12.3, h: 5.3 },
  });

  s.addText([
    { text: "Verdict:  ", options: { bold: true, color: C.navy } },
    { text: "per-CpG concordance with methylKit is r = 1.000 on the 15,688-site intersection (top). DMR-level Δβ concordance is r = 0.887 vs methylKit, r = 0.972 vs DSS, both with 100% direction agreement on matched DMRs (bottom).", options: { color: C.ink } },
  ], { x: 0.5, y: 6.75, w: 12.3, h: 0.45, fontSize: 11, fontFace: F.body, margin: 0 });

  pageFooter(s, "Study 3 — agreement", 19);
  speakerNote(s,
    "Per-CpG: r = 1.000. Sites that pass q ≤ 0.05 in both tools agree perfectly on direction AND magnitude. " +
    "DMR level: r = 0.972 vs DSS, 0.887 vs methylKit. " +
    "Disagreements are about WHERE to draw region boundaries, not about the sign of the methylation change."
  );
}

// =========================================================================
// 20 — F05: DMR overlap with paper
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(s,
    "Recovery of the published DMR set",
    "How many of Farhangdoost et al.'s 813 DMRs does each tool independently rediscover?"
  );

  s.addImage({
    path: HEAD("F05_dmr_overlap_recall.png"),
    x: 0.5, y: 1.45, w: 12.3, h: 5.2, sizing: { type: "contain", w: 12.3, h: 5.2 },
  });

  s.addText([
    { text: "epykit chain_merge recovers 428 / 813 ≈ 52.6% ", options: { bold: true, color: C.navy } },
    { text: "of the paper's DMRs with 100% direction agreement on matched intervals. The remaining 47% are mostly short-length DMRs (paper median = 239 bp; ours = 123 bp) — a morphology, not a biology, gap.", options: { color: C.ink } },
  ], { x: 0.5, y: 6.7, w: 12.3, h: 0.55, fontSize: 11, fontFace: F.body, margin: 0 });

  pageFooter(s, "Study 3 — DMR recall", 20);
  speakerNote(s,
    "Half-recall sounds modest until you remember the paper used DSS with smoothing while we used chain_merge without smoothing. " +
    "The disagreement is mostly about where to draw region boundaries. " +
    "100% direction agreement means we never call a DMR with the wrong sign — when both call the same locus."
  );
}

// =========================================================================
// 21 — F07: Pipeline step breakdown
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(s,
    "Where does each tool spend its time?",
    "Per-step wall-clock decomposition on GSE263850"
  );

  s.addImage({
    path: HEAD("F07_pipeline_step_breakdown.png"),
    x: 0.5, y: 1.45, w: 12.3, h: 5.2, sizing: { type: "contain", w: 12.3, h: 5.2 },
  });

  s.addText([
    { text: "The bottleneck is the DMC test:  ", options: { bold: true, color: C.navy } },
    { text: "methylKit's `calculateDiffMeth` step alone takes 177 minutes on 14 M CpGs; epykit's equivalent takes 3.5 minutes. The other 11 steps (read, QC, filter, unite, tile, annotate) take comparable time across both tools.", options: { color: C.ink } },
  ], { x: 0.5, y: 6.7, w: 12.3, h: 0.55, fontSize: 11, fontFace: F.body, margin: 0 });

  pageFooter(s, "Study 3 — where time goes", 21);
  speakerNote(s,
    "This is where the speed gap actually comes from. " +
    "methylKit spends 90% of its time in the DMC step — calling per-CpG logistic regressions one CpG at a time in R. " +
    "epykit vectorises that same computation in NumPy. Same maths, 50× less wall-clock."
  );
}

// =========================================================================
// 22 — F03: FDR calibration
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(s,
    "Are the q-values honest? Permutation null calibration",
    "Shuffle labels, re-run, count false discoveries at q ≤ 0.05"
  );

  s.addImage({
    path: HEAD("F03_fdr_calibration_forest.png"),
    x: 0.5, y: 1.45, w: 12.3, h: 5.2, sizing: { type: "contain", w: 12.3, h: 5.2 },
  });

  s.addText([
    { text: "All engines stay well below the nominal 0.05 line ", options: { bold: true, color: C.navy } },
    { text: "across three datasets (real GSE263850, Piao 2021 simulator, and our internal simulator). Observed FDR ranges from 0 (welch_t, glm) to 1.5 × 10⁻⁵ (lr, lr+) — q-values can be trusted.", options: { color: C.ink } },
  ], { x: 0.5, y: 6.7, w: 12.3, h: 0.55, fontSize: 11, fontFace: F.body, margin: 0 });

  pageFooter(s, "Robustness — calibration", 22);
  speakerNote(s,
    "Standard sanity check: permute the case/control labels, re-run, count how many sites still come back significant. " +
    "All engines deliver observed FDR < 0.0001 on all three datasets. " +
    "q-values are honest. This isn't a flashy figure but it's load-bearing for the paper."
  );
}

// =========================================================================
// 23 — F06: lr vs lr+
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(s,
    "The lr+ power stack — does it actually help?",
    "Comparing bare `lr` to the four-knob `lr+` bundle"
  );

  s.addImage({
    path: HEAD("F06_lr_vs_lrplus_ablation.png"),
    x: 0.5, y: 1.45, w: 12.3, h: 5.2, sizing: { type: "contain", w: 12.3, h: 5.2 },
  });

  s.addText([
    { text: "lr+ buys ~15-20 pp TPR at coverage 5× ", options: { bold: true, color: C.navy } },
    { text: "(top-left), at the cost of a small FPR uptick in the low-coverage / low-replicate regime (bottom panels). At cov ≥ 15× or n ≥ 6, lr and lr+ are indistinguishable — only use lr+ when low-power scenarios bite.", options: { color: C.ink } },
  ], { x: 0.5, y: 6.7, w: 12.3, h: 0.55, fontSize: 11, fontFace: F.body, margin: 0 });

  pageFooter(s, "Robustness — lr+ benefit", 23);
  speakerNote(s,
    "The lr+ stack pays off at low coverage and low replicate counts. " +
    "At normal coverages (15× and up) it's not worth the overhead. " +
    "Mention that lr+ is opt-in by design — bare lr is the default, lr+ has to be asked for."
  );
}

// =========================================================================
// 24 — F10: lr+ tradeoff
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(s,
    "When the lr+ stack costs you — the honesty slide",
    "Neighbour-combine inflates FPR when DMCs cluster densely"
  );

  s.addImage({
    path: HEAD("F10_lrplus_tradeoff.png"),
    x: 0.5, y: 1.45, w: 12.3, h: 5.2, sizing: { type: "contain", w: 12.3, h: 5.2 },
  });

  s.addText([
    { text: "Known calibration issue:  ", options: { bold: true, color: C.navy } },
    { text: "lr+'s sign-aware Stouffer combiner amplifies clustered same-direction neighbours — when DMCs are densely packed (eg whole-DMR-embedded simulator scenarios), realised FDR can reach ~25%. The simulator-grid scenarios stay below 0.05 but the regime is documented in paper §4.", options: { color: C.ink } },
  ], { x: 0.5, y: 6.7, w: 12.3, h: 0.55, fontSize: 11, fontFace: F.body, margin: 0 });

  pageFooter(s, "Robustness — lr+ trade-off", 24);
  speakerNote(s,
    "Honesty slide. lr+ has a known failure mode when DMCs cluster densely — the Stouffer combiner over-counts evidence from neighbouring sites. " +
    "Documented in the paper. Mitigation is to use bare lr in those settings."
  );
}

// =========================================================================
// 25 — Limitations
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(s,
    "Limitations — when not to use epykit",
    "Honest scope boundaries"
  );

  const lims = [
    { title: "Not for n = 2",
      desc: "epykit's lr and lr+ engines require n ≥ 4. At n_total = 2, only the Fisher fallback works — and methylKit_tuned outperforms it. Use methylKit for paired n = 2 studies." },
    { title: "Not a smoothing tool",
      desc: "DSS and BSmooth do local-smoothing of β-values; epykit does not. For very low coverage (< 5×) where smoothing is the right call, use DSS." },
    { title: "Reduced-Representation BS not first-class",
      desc: "RRBS data works, but the unite step assumes WGBS-like density. For RRBS-specific designs, BiSeq remains the recommendation." },
    { title: "lr+ is not free",
      desc: "The four-knob stack has documented calibration issues on signal-dense simulators. Stay with bare lr for routine analysis; reach for lr+ only when low-power scenarios bite." },
    { title: "Linux extras",
      desc: "pyBigWig, methylKit interop, and BAM-input via pysam are Linux/macOS only. Windows installs use the pure-Python core." },
    { title: "GPU acceleration is preliminary",
      desc: "_glm_gpu (CuPy / JAX) exists but is experimental; CPU path is the supported route at v1.0." },
  ];

  const gridX0 = 0.6;
  const gridY0 = 1.55;
  const cardW = 4.0;
  const cardH = 2.55;
  const colGap = 0.10;
  const rowGap = 0.15;

  lims.forEach((f, i) => {
    const col = i % 3;
    const row = Math.floor(i / 3);
    const x = gridX0 + col * (cardW + colGap);
    const y = gridY0 + row * (cardH + rowGap);

    s.addShape(pres.shapes.RECTANGLE, { x, y, w: cardW, h: cardH,
      fill: { color: C.cardBg }, line: { color: C.rule, width: 1 } });
    s.addShape(pres.shapes.RECTANGLE, { x, y, w: 0.10, h: cardH,
      fill: { color: C.warn }, line: { type: "none" } });
    s.addText(f.title, {
      x: x + 0.25, y: y + 0.2, w: cardW - 0.4, h: 0.5,
      fontSize: 13, fontFace: F.head, bold: true, color: C.navy, margin: 0,
    });
    s.addText(f.desc, {
      x: x + 0.25, y: y + 0.85, w: cardW - 0.4, h: cardH - 1.0,
      fontSize: 10.5, fontFace: F.body, color: C.ink, margin: 0,
    });
  });

  pageFooter(s, "Wrap-up — limitations", 25);
  speakerNote(s,
    "Six honest limits. Lead with n=2 (the methylKit-tuned case from earlier). " +
    "Mention that 'when not to use' slides typically de-risk the work — advisors trust talks that acknowledge limits."
  );
}

// =========================================================================
// 26 — Future work + thanks
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.navyDeep };
  s.addText("Future work + acknowledgments", {
    x: 0.6, y: 0.32, w: W - 1.2, h: 0.7,
    fontSize: 28, fontFace: F.head, bold: true, color: C.white, margin: 0,
  });

  s.addShape(pres.shapes.RECTANGLE, { x: 0.6, y: 1.55, w: 6.05, h: 5.4,
    fill: { color: C.navy }, line: { color: C.gold, width: 1 } });
  s.addText("v1.1 roadmap", {
    x: 0.8, y: 1.65, w: 5.8, h: 0.4,
    fontSize: 16, fontFace: F.head, bold: true, color: C.gold, margin: 0,
  });
  s.addText([
    { text: "•  Single-cell methylation (sc-WGBS): AnnData + MuData first-class, scvi-tools integration.", options: { color: "DCE7EF", breakLine: true } },
    { text: "", options: { breakLine: true } },
    { text: "•  Allele-specific methylation engines (haplotype-aware DMR calling).", options: { color: "DCE7EF", breakLine: true } },
    { text: "", options: { breakLine: true } },
    { text: "•  GPU path for the IRLS GLM (currently experimental → supported).", options: { color: "DCE7EF", breakLine: true } },
    { text: "", options: { breakLine: true } },
    { text: "•  Per-knob ablation of the lr+ stack (decompose the four-knob bundle into single-knob attribution).", options: { color: "DCE7EF", breakLine: true } },
    { text: "", options: { breakLine: true } },
    { text: "•  CLI surface for lr+ knobs (currently API-only).", options: { color: "DCE7EF", breakLine: true } },
    { text: "", options: { breakLine: true } },
    { text: "•  Benchmark v2: head-to-head against MOABS, biscuit, metilene.", options: { color: "DCE7EF" } },
  ], { x: 0.8, y: 2.15, w: 5.7, h: 4.7, fontSize: 12, fontFace: F.body, margin: 0 });

  s.addShape(pres.shapes.RECTANGLE, { x: 6.85, y: 1.55, w: 5.85, h: 5.4,
    fill: { color: "0E3656" }, line: { color: C.gold, width: 1 } });
  s.addText("Acknowledgments", {
    x: 7.05, y: 1.65, w: 5.5, h: 0.4,
    fontSize: 16, fontFace: F.head, bold: true, color: C.gold, margin: 0,
  });
  s.addText([
    { text: "•  Master's advisor — guidance, scope, and reading drafts.", options: { color: "DCE7EF", breakLine: true } },
    { text: "", options: { breakLine: true } },
    { text: "•  Authors of methylKit, DSS, RADMeth, BSmooth — the tools that defined the field and against which epykit is benchmarked.", options: { color: "DCE7EF", breakLine: true } },
    { text: "", options: { breakLine: true } },
    { text: "•  Piao et al. 2021 for an extensively documented simulator that made apples-to-apples comparison possible.", options: { color: "DCE7EF", breakLine: true } },
    { text: "", options: { breakLine: true } },
    { text: "•  Farhangdoost et al. 2024 for the GSE263850 dataset with a complete reproducible analysis trail.", options: { color: "DCE7EF", breakLine: true } },
    { text: "", options: { breakLine: true } },
    { text: "•  The Polars, scanpy, and scientific-Python communities for the abstractions epykit stands on.", options: { color: "DCE7EF" } },
  ], { x: 7.05, y: 2.15, w: 5.55, h: 4.7, fontSize: 11.5, fontFace: F.body, margin: 0 });

  s.addText("Thank you  —  questions?", {
    x: 0.6, y: 7.0, w: W - 1.2, h: 0.4,
    fontSize: 18, fontFace: F.head, italic: true, color: C.gold, align: "center", margin: 0,
  });

  speakerNote(s,
    "Close with the roadmap and acknowledgments. The single-cell direction is the natural follow-up — " +
    "single-cell methylation tooling is even sparser than bulk."
  );
}

// =========================================================================
// Save
// =========================================================================
pres.writeFile({ fileName: path.join(ROOT, "epykit_thesis_presentation.pptx") })
  .then((p) => {
    console.log(`\nWrote: ${p}`);
  })
  .catch((err) => {
    console.error("ERROR:", err);
    process.exit(1);
  });
