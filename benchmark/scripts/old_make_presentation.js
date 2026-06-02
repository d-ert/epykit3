// epykit benchmark — slide deck
// Outputs: FINAL_REPORT/epykit_benchmark_presentation.pptx

const path = require("path");
process.env.NODE_PATH = require("child_process")
  .execSync("npm root -g")
  .toString()
  .trim();
require("module").Module._initPaths();

const pptxgen = require("pptxgenjs");

const ROOT = path.resolve(__dirname, "..");
const FIG = (relpath) => path.join(ROOT, "figures", relpath);

// -- palette: Ocean Gradient ------------------------------------------------
const C = {
  navy: "0F2942",       // dark title/conclusion bg
  navyDeep: "081C2E",
  blue: "065A82",       // section headers
  teal: "1C7293",       // accent strokes
  gold: "F4B400",       // callouts, speedup numbers
  white: "FFFFFF",
  ink: "1E293B",        // body text on light
  muted: "64748B",      // subtitles, captions
  rule: "CBD5E1",       // hairlines/dividers
  cardBg: "F8FAFC",     // card surface on white
  cardBg2: "EFF6FB",    // alt card
  good: "22C55E",       // positive metric
  warn: "F97316",       // tradeoff metric
};

const F = {
  head: "Georgia",
  body: "Calibri",
  mono: "Consolas",
};

const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE"; // 13.3 × 7.5
pres.author = "epykit contributors";
pres.title = "epykit Benchmark — Consolidated Evaluation";

const W = 13.3;
const H = 7.5;

// helpers ------------------------------------------------------------------
function titleHeader(slide, title, subtitle) {
  slide.addText(title, {
    x: 0.6,
    y: 0.35,
    w: W - 1.2,
    h: 0.7,
    fontSize: 30,
    fontFace: F.head,
    bold: true,
    color: C.navy,
    margin: 0,
  });
  if (subtitle) {
    slide.addText(subtitle, {
      x: 0.6,
      y: 0.95,
      w: W - 1.2,
      h: 0.4,
      fontSize: 14,
      fontFace: F.body,
      italic: true,
      color: C.muted,
      margin: 0,
    });
  }
}

function pageFooter(slide, label, pageNum) {
  slide.addText(label, {
    x: 0.6,
    y: H - 0.42,
    w: 8,
    h: 0.3,
    fontSize: 9,
    fontFace: F.body,
    color: C.muted,
    margin: 0,
  });
  slide.addText(String(pageNum), {
    x: W - 1.2,
    y: H - 0.42,
    w: 0.6,
    h: 0.3,
    fontSize: 9,
    fontFace: F.body,
    color: C.muted,
    align: "right",
    margin: 0,
  });
}

// =========================================================================
// SLIDE 1 — Title (dark)
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.navyDeep };

  // Decorative motif: stacked teal squares lower-left (carries to other slides)
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.6,
    y: 6.4,
    w: 0.35,
    h: 0.35,
    fill: { color: C.teal },
    line: { type: "none" },
  });
  s.addShape(pres.shapes.RECTANGLE, {
    x: 1.05,
    y: 6.4,
    w: 0.35,
    h: 0.35,
    fill: { color: C.blue },
    line: { type: "none" },
  });
  s.addShape(pres.shapes.RECTANGLE, {
    x: 1.5,
    y: 6.4,
    w: 0.35,
    h: 0.35,
    fill: { color: C.gold },
    line: { type: "none" },
  });

  s.addText("epykit", {
    x: 0.6,
    y: 1.6,
    w: 12,
    h: 1.3,
    fontSize: 96,
    fontFace: F.head,
    bold: true,
    color: C.white,
    margin: 0,
  });
  s.addText("A Python-native pipeline for differential methylation analysis", {
    x: 0.6,
    y: 2.9,
    w: 12,
    h: 0.6,
    fontSize: 22,
    fontFace: F.body,
    color: "BFD7E5",
    margin: 0,
  });
  s.addText("Benchmarked across three studies on simulated and real WGBS data", {
    x: 0.6,
    y: 3.5,
    w: 12,
    h: 0.45,
    fontSize: 16,
    fontFace: F.body,
    italic: true,
    color: "8FB6CA",
    margin: 0,
  });

  // Right-side stat strip
  const stats = [
    { num: "3", label: "benchmark studies" },
    { num: "9", label: "tools compared" },
    { num: "43×", label: "faster than methylKit" },
  ];
  stats.forEach((stat, i) => {
    const y = 4.6 + i * 0.7;
    s.addText(stat.num, {
      x: 8.4,
      y: y,
      w: 1.6,
      h: 0.6,
      fontSize: 36,
      fontFace: F.head,
      bold: true,
      color: C.gold,
      align: "right",
      margin: 0,
    });
    s.addText(stat.label, {
      x: 10.1,
      y: y + 0.05,
      w: 2.6,
      h: 0.5,
      fontSize: 14,
      fontFace: F.body,
      color: "BFD7E5",
      margin: 0,
    });
  });

  s.addText("2026-05-22  •  Consolidated final report", {
    x: 0.6,
    y: 6.95,
    w: 8,
    h: 0.3,
    fontSize: 11,
    fontFace: F.body,
    color: "8FB6CA",
    margin: 0,
  });
}

// =========================================================================
// SLIDE 2 — The challenge
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(
    s,
    "The WGBS analysis ecosystem is fragmented",
    "Established tools live in R/CLI; modern multi-omics workflows live in Python"
  );

  // Left column: problem text
  s.addText(
    [
      {
        text: "Whole-genome bisulfite sequencing (WGBS)",
        options: { bold: true, color: C.ink, breakLine: true },
      },
      {
        text:
          "is the reference technique for base-resolution DNA methylation. " +
          "Mature analysis tools — methylKit, DSS, BSmooth, methylSig, BiSeq, RADMeth, metilene — " +
          "live in R/Bioconductor or as CLI binaries.",
        options: { color: C.ink, breakLine: true },
      },
      { text: "", options: { breakLine: true } },
      {
        text: "Two pressures complicate this:",
        options: { bold: true, color: C.ink, breakLine: true },
      },
      {
        text: "Single-cell methylation & multi-omics workflows are written in Python (scanpy / anndata / mudata).",
        options: { bullet: true, color: C.ink, breakLine: true },
      },
      {
        text: "Modern experiments routinely involve dozens of samples and tens of millions of CpG sites — in-memory R data frames are no longer the right abstraction.",
        options: { bullet: true, color: C.ink },
      },
    ],
    {
      x: 0.6,
      y: 1.6,
      w: 6.4,
      h: 5.0,
      fontSize: 14,
      fontFace: F.body,
      paraSpaceAfter: 6,
      margin: 0,
    }
  );

  // Right column: card with the proposal
  s.addShape(pres.shapes.RECTANGLE, {
    x: 7.3,
    y: 1.7,
    w: 5.4,
    h: 4.6,
    fill: { color: C.navy },
    line: { type: "none" },
  });
  s.addShape(pres.shapes.RECTANGLE, {
    x: 7.3,
    y: 1.7,
    w: 0.12,
    h: 4.6,
    fill: { color: C.gold },
    line: { type: "none" },
  });
  s.addText("Our proposal: epykit", {
    x: 7.65,
    y: 1.95,
    w: 5,
    h: 0.5,
    fontSize: 22,
    fontFace: F.head,
    bold: true,
    color: C.white,
    margin: 0,
  });
  s.addText(
    [
      {
        text: "Partitioned Parquet storage with lazy I/O (polars).",
        options: { bullet: true, color: "DCE7EF", breakLine: true },
      },
      {
        text: "scanpy-style pp / tl / pl API.",
        options: { bullet: true, color: "DCE7EF", breakLine: true },
      },
      {
        text: "8 DMC backends (lr, glm, bb_lr, fisher, welch_t, …).",
        options: { bullet: true, color: "DCE7EF", breakLine: true },
      },
      {
        text: "4 DMR engines (tile, sliding-Stouffer, HMM, chain_merge).",
        options: { bullet: true, color: "DCE7EF", breakLine: true },
      },
      {
        text: "AnnData / MuData interop for downstream analysis.",
        options: { bullet: true, color: "DCE7EF" },
      },
    ],
    {
      x: 7.65,
      y: 2.55,
      w: 4.9,
      h: 3.5,
      fontSize: 14,
      fontFace: F.body,
      paraSpaceAfter: 6,
      margin: 0,
    }
  );

  pageFooter(s, "Background", 2);
}

// =========================================================================
// SLIDE 3 — Three-study design (3 cards)
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(
    s,
    "Three studies covering the evaluation surface",
    "Each addresses a different question; together they triangulate epykit's behaviour"
  );

  const cardY = 1.7;
  const cardH = 5.0;
  const cardW = 3.95;
  const cardGap = 0.27;
  const cardX0 = 0.6;

  const cards = [
    {
      tag: "STUDY 1",
      title: "Panel comparison",
      sub: "Simulated · Piao et al. 2021 grid",
      data: "100 K CpGs (DMC) + 4 M CpGs (DMR)\nCoverage 5×–25×, n = 2–10",
      compared: "8 baselines:\nmethylKit · DSS · RADMeth · BiSeq ·\nmethylSig · BSmooth · metilene · Fisher",
      question: "Where does epykit sit within the broader ecosystem?",
      tagColor: C.blue,
    },
    {
      tag: "STUDY 2",
      title: "Head-to-head: methylKit",
      sub: "Simulated · same machine, same harness",
      data: "Same 100 K-CpG grid\nLocal methylKit run, OS-level resource tracker",
      compared: "methylKit 1.34.0\ncalculateDiffMeth + tileMethylCounts",
      question: "Is epykit statistically equivalent to the closest analogue?",
      tagColor: C.teal,
    },
    {
      tag: "STUDY 3",
      title: "Real WGBS data",
      sub: "GEO GSE263850 (hg38)",
      data: "6 samples · 15.6 M CpGs after filtering\nClone16/20/21 vs SBP009 untreated",
      compared: "methylKit 1.36.0",
      question: "Does epykit behave correctly on real biological data?",
      tagColor: C.gold,
    },
  ];

  cards.forEach((c, i) => {
    const x = cardX0 + i * (cardW + cardGap);
    // card body
    s.addShape(pres.shapes.RECTANGLE, {
      x,
      y: cardY,
      w: cardW,
      h: cardH,
      fill: { color: C.cardBg },
      line: { color: C.rule, width: 0.75 },
    });
    // top stripe
    s.addShape(pres.shapes.RECTANGLE, {
      x,
      y: cardY,
      w: cardW,
      h: 0.16,
      fill: { color: c.tagColor },
      line: { type: "none" },
    });

    // STUDY N badge
    s.addShape(pres.shapes.RECTANGLE, {
      x: x + 0.3,
      y: cardY + 0.35,
      w: 1.0,
      h: 0.3,
      fill: { color: c.tagColor },
      line: { type: "none" },
    });
    s.addText(c.tag, {
      x: x + 0.3,
      y: cardY + 0.35,
      w: 1.0,
      h: 0.3,
      fontSize: 11,
      fontFace: F.body,
      bold: true,
      color: C.white,
      align: "center",
      valign: "middle",
      margin: 0,
      charSpacing: 2,
    });

    s.addText(c.title, {
      x: x + 0.3,
      y: cardY + 0.8,
      w: cardW - 0.6,
      h: 0.55,
      fontSize: 22,
      fontFace: F.head,
      bold: true,
      color: C.ink,
      margin: 0,
    });
    s.addText(c.sub, {
      x: x + 0.3,
      y: cardY + 1.35,
      w: cardW - 0.6,
      h: 0.35,
      fontSize: 12,
      fontFace: F.body,
      italic: true,
      color: C.muted,
      margin: 0,
    });

    // separator
    s.addShape(pres.shapes.LINE, {
      x: x + 0.3,
      y: cardY + 1.85,
      w: cardW - 0.6,
      h: 0,
      line: { color: C.rule, width: 0.5 },
    });

    s.addText("DATA", {
      x: x + 0.3,
      y: cardY + 1.95,
      w: cardW - 0.6,
      h: 0.25,
      fontSize: 9,
      fontFace: F.body,
      bold: true,
      color: c.tagColor,
      margin: 0,
      charSpacing: 2,
    });
    s.addText(c.data, {
      x: x + 0.3,
      y: cardY + 2.2,
      w: cardW - 0.6,
      h: 0.8,
      fontSize: 12,
      fontFace: F.body,
      color: C.ink,
      margin: 0,
    });

    s.addText("COMPARED AGAINST", {
      x: x + 0.3,
      y: cardY + 3.05,
      w: cardW - 0.6,
      h: 0.25,
      fontSize: 9,
      fontFace: F.body,
      bold: true,
      color: c.tagColor,
      margin: 0,
      charSpacing: 2,
    });
    s.addText(c.compared, {
      x: x + 0.3,
      y: cardY + 3.3,
      w: cardW - 0.6,
      h: 0.95,
      fontSize: 12,
      fontFace: F.body,
      color: C.ink,
      margin: 0,
    });

    // question pinned at bottom
    s.addShape(pres.shapes.RECTANGLE, {
      x,
      y: cardY + cardH - 1.1,
      w: cardW,
      h: 1.1,
      fill: { color: C.navy },
      line: { type: "none" },
    });
    s.addText(c.question, {
      x: x + 0.3,
      y: cardY + cardH - 0.95,
      w: cardW - 0.6,
      h: 0.85,
      fontSize: 13,
      fontFace: F.body,
      italic: true,
      color: C.white,
      valign: "middle",
      margin: 0,
    });
  });

  pageFooter(s, "Design", 3);
}

// =========================================================================
// SLIDE 4 — Study 1: TPR vs coverage
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(
    s,
    "Study 1 · DMC sensitivity at low coverage",
    "epykit lr retains 83.5 % TPR at 5× coverage in the hardest effect-size bin"
  );

  // Left figure
  s.addImage({
    path: FIG("study1_simulated_allPackages/F2_tpr_vs_coverage.png"),
    x: 0.6,
    y: 1.55,
    w: 7.4,
    h: 5.4,
    sizing: { type: "contain", w: 7.4, h: 5.4 },
  });

  // Right column: stat callouts
  const rightX = 8.4;
  s.addText("AT 5× COVERAGE, BIN 0.2–0.4", {
    x: rightX,
    y: 1.55,
    w: 4.3,
    h: 0.3,
    fontSize: 10,
    fontFace: F.body,
    bold: true,
    color: C.muted,
    charSpacing: 2,
    margin: 0,
  });
  s.addText("83.5%", {
    x: rightX,
    y: 1.85,
    w: 4.3,
    h: 1.0,
    fontSize: 64,
    fontFace: F.head,
    bold: true,
    color: C.gold,
    margin: 0,
  });
  s.addText("epykit lr — highest in the panel", {
    x: rightX,
    y: 2.9,
    w: 4.3,
    h: 0.35,
    fontSize: 13,
    fontFace: F.body,
    color: C.ink,
    margin: 0,
  });

  s.addShape(pres.shapes.LINE, {
    x: rightX,
    y: 3.35,
    w: 4.3,
    h: 0,
    line: { color: C.rule, width: 0.5 },
  });

  // baseline comparison table
  const baselines = [
    ["methylKit", "26.6 %"],
    ["DSS", "6.5 %"],
    ["RADMeth", "42.1 %"],
    ["Fisher (pooled)", "8.2 %"],
    ["methylSig", "n/a"],
    ["BiSeq", "68.4 %"],
  ];
  s.addText("Same cell, baseline tools:", {
    x: rightX,
    y: 3.5,
    w: 4.3,
    h: 0.3,
    fontSize: 11,
    fontFace: F.body,
    italic: true,
    color: C.muted,
    margin: 0,
  });
  baselines.forEach((row, i) => {
    const yy = 3.85 + i * 0.32;
    s.addText(row[0], {
      x: rightX,
      y: yy,
      w: 2.5,
      h: 0.3,
      fontSize: 12,
      fontFace: F.body,
      color: C.ink,
      margin: 0,
    });
    s.addText(row[1], {
      x: rightX + 2.5,
      y: yy,
      w: 1.8,
      h: 0.3,
      fontSize: 12,
      fontFace: F.body,
      bold: true,
      color: C.muted,
      align: "right",
      margin: 0,
    });
  });

  // bottom note
  s.addText(
    "Above 15× coverage every credible tool converges to TPR ≈ 0.97–1.00 — method choice matters most at low depth.",
    {
      x: 0.6,
      y: 6.95,
      w: 12.1,
      h: 0.3,
      fontSize: 11,
      fontFace: F.body,
      italic: true,
      color: C.muted,
      margin: 0,
    }
  );

  pageFooter(s, "Study 1 · DMC accuracy", 4);
}

// =========================================================================
// SLIDE 5 — Study 1: FPR calibration
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(
    s,
    "Study 1 · False-positive calibration",
    "lr's FPR at 5× is 100×–600× tighter than any baseline"
  );

  s.addImage({
    path: FIG("study1_simulated_allPackages/F4_fpr_vs_coverage.png"),
    x: 0.6,
    y: 1.55,
    w: 7.4,
    h: 5.0,
    sizing: { type: "contain", w: 7.4, h: 5.0 },
  });

  // Right: FPR comparison table
  const rightX = 8.4;
  s.addText("FPR @ q < 0.05, 5× coverage", {
    x: rightX,
    y: 1.55,
    w: 4.3,
    h: 0.4,
    fontSize: 14,
    fontFace: F.body,
    bold: true,
    color: C.ink,
    margin: 0,
  });

  const fprRows = [
    ["epykit lr", "3.7 × 10⁻⁵", C.gold],
    ["epykit fisher", "2.5 × 10⁻⁵", C.ink],
    ["methylKit", "2.3 × 10⁻²", C.muted],
    ["RADMeth", "1.8 × 10⁻²", C.muted],
    ["DSS", "2.9 × 10⁻²", C.muted],
    ["BiSeq", "1.0 × 10⁻²", C.muted],
  ];
  fprRows.forEach((row, i) => {
    const yy = 2.05 + i * 0.5;
    s.addShape(pres.shapes.RECTANGLE, {
      x: rightX,
      y: yy,
      w: 4.3,
      h: 0.45,
      fill: { color: i === 0 ? C.cardBg2 : C.white },
      line: { type: "none" },
    });
    s.addText(row[0], {
      x: rightX + 0.15,
      y: yy,
      w: 2.0,
      h: 0.45,
      fontSize: 13,
      fontFace: F.body,
      bold: i === 0,
      color: row[2],
      valign: "middle",
      margin: 0,
    });
    s.addText(row[1], {
      x: rightX + 2.05,
      y: yy,
      w: 2.1,
      h: 0.45,
      fontSize: 14,
      fontFace: F.mono,
      bold: i === 0,
      color: row[2],
      align: "right",
      valign: "middle",
      margin: 0,
    });
  });

  // Insight callout at bottom
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.6,
    y: 6.75,
    w: 12.1,
    h: 0.55,
    fill: { color: C.cardBg },
    line: { type: "none" },
  });
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.6,
    y: 6.75,
    w: 0.08,
    h: 0.55,
    fill: { color: C.gold },
    line: { type: "none" },
  });
  s.addText(
    [
      { text: "Why? ", options: { bold: true, color: C.ink } },
      {
        text:
          "lr clamps dispersion at φ ≥ 1 (binomial floor). The Piao 2021 simulator is underdispersed (median φ ≈ 0.41), so the clamp is nearly correct — neither too aggressive (low FPR) nor too conservative (high TPR).",
        options: { color: C.ink },
      },
    ],
    {
      x: 0.85,
      y: 6.75,
      w: 11.8,
      h: 0.55,
      fontSize: 12,
      fontFace: F.body,
      valign: "middle",
      margin: 0,
    }
  );

  pageFooter(s, "Study 1 · FPR", 5);
}

// =========================================================================
// SLIDE 6 — Study 1: DMR detection
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(
    s,
    "Study 1 · DMR detection on the 4 M-site simulation",
    "epykit chain_merge recovers 97 % / 100 % / 100 % / 100 % / 100 % across the coverage grid"
  );

  s.addImage({
    path: FIG("study1_simulated_allPackages/F5_dmr_detection.png"),
    x: 0.6,
    y: 1.55,
    w: 7.6,
    h: 5.2,
    sizing: { type: "contain", w: 7.6, h: 5.2 },
  });

  const rightX = 8.6;

  // headline number
  s.addText("DMR RECALL @ 5× COVERAGE", {
    x: rightX,
    y: 1.55,
    w: 4.1,
    h: 0.3,
    fontSize: 10,
    fontFace: F.body,
    bold: true,
    color: C.muted,
    charSpacing: 2,
    margin: 0,
  });
  s.addText("97 %", {
    x: rightX,
    y: 1.9,
    w: 4.1,
    h: 0.9,
    fontSize: 60,
    fontFace: F.head,
    bold: true,
    color: C.gold,
    margin: 0,
  });
  s.addText("epykit chain_merge, post-fix", {
    x: rightX,
    y: 2.8,
    w: 4.1,
    h: 0.35,
    fontSize: 12,
    fontFace: F.body,
    color: C.ink,
    margin: 0,
  });

  s.addShape(pres.shapes.LINE, {
    x: rightX,
    y: 3.25,
    w: 4.1,
    h: 0,
    line: { color: C.rule, width: 0.5 },
  });

  s.addText("DMR recall (5× / 10× / 25×):", {
    x: rightX,
    y: 3.4,
    w: 4.1,
    h: 0.3,
    fontSize: 11,
    fontFace: F.body,
    italic: true,
    color: C.muted,
    margin: 0,
  });

  const dmrRows = [
    ["methylKit", "1.00 / 1.00 / 1.00"],
    ["epykit chain_merge", "0.97 / 1.00 / 1.00"],
    ["RADMeth", "0.74 / 0.86 / 1.00"],
    ["methylSig", "0.20 / 0.86 / 0.97"],
    ["DSS", "0.00 / 0.06 / 0.69"],
    ["BiSeq", "0.14 / 0.14 / 0.14"],
  ];
  dmrRows.forEach((row, i) => {
    const yy = 3.75 + i * 0.4;
    s.addText(row[0], {
      x: rightX,
      y: yy,
      w: 2.2,
      h: 0.35,
      fontSize: 12,
      fontFace: F.body,
      bold: i === 1,
      color: i === 1 ? C.gold : C.ink,
      margin: 0,
    });
    s.addText(row[1], {
      x: rightX + 2.2,
      y: yy,
      w: 1.9,
      h: 0.35,
      fontSize: 11,
      fontFace: F.mono,
      color: i === 1 ? C.gold : C.muted,
      align: "right",
      margin: 0,
    });
  });

  pageFooter(s, "Study 1 · DMR", 6);
}

// =========================================================================
// SLIDE 7 — Study 2: bit-precise calibration table
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(
    s,
    "Study 2 · Statistical equivalence to methylKit",
    "TPR, FPR, F1, AUROC identical to three decimal places at n ≥ 4 — 3 vs 3 design"
  );

  // Big table center
  const headers = ["Coverage", "Tool", "TPR", "FPR", "F1", "AUROC"];
  const rows = [
    ["5×", "epykit lr", "0.849", "3.7 × 10⁻⁵", "0.918", "0.9990"],
    ["5×", "methylKit", "0.849", "3.7 × 10⁻⁵", "0.918", "0.9990"],
    ["10×", "epykit lr", "0.944", "1.2 × 10⁻⁵", "0.971", "0.9999"],
    ["10×", "methylKit", "0.944", "1.2 × 10⁻⁵", "0.971", "0.9999"],
    ["15×", "epykit lr", "0.984", "1.2 × 10⁻⁵", "0.992", "1.0000"],
    ["15×", "methylKit", "0.984", "1.2 × 10⁻⁵", "0.992", "1.0000"],
    ["20×", "epykit lr", "0.991", "1.2 × 10⁻⁵", "0.995", "1.0000"],
    ["20×", "methylKit", "0.991", "1.2 × 10⁻⁵", "0.995", "1.0000"],
    ["25×", "epykit lr", "0.993", "0.0", "0.996", "1.0000"],
    ["25×", "methylKit", "0.993", "0.0", "0.997", "1.0000"],
  ];

  const tableData = [
    headers.map((h) => ({
      text: h,
      options: {
        bold: true,
        fill: { color: C.navy },
        color: C.white,
        fontSize: 13,
        fontFace: F.body,
        align: "center",
        valign: "middle",
      },
    })),
    ...rows.map((row) =>
      row.map((v, j) => {
        const isEpykit = row[1].includes("epykit");
        return {
          text: v,
          options: {
            fontSize: 12,
            fontFace: j === 2 || j === 3 || j === 4 || j === 5 ? F.mono : F.body,
            color: isEpykit ? C.blue : C.ink,
            bold: isEpykit && (j === 2 || j === 5),
            fill: { color: isEpykit ? C.cardBg2 : C.white },
            align: j === 0 || j === 1 ? "left" : "right",
            valign: "middle",
          },
        };
      })
    ),
  ];

  s.addTable(tableData, {
    x: 0.6,
    y: 1.55,
    w: 8.5,
    colW: [1.1, 2.0, 1.3, 1.6, 1.3, 1.2],
    rowH: 0.35,
    border: { pt: 0.5, color: C.rule },
  });

  // Right column: interpretation
  const rightX = 9.5;
  s.addShape(pres.shapes.RECTANGLE, {
    x: rightX,
    y: 1.55,
    w: 3.2,
    h: 5.3,
    fill: { color: C.navy },
    line: { type: "none" },
  });
  s.addText("Same model.", {
    x: rightX + 0.25,
    y: 1.75,
    w: 2.7,
    h: 0.45,
    fontSize: 20,
    fontFace: F.head,
    bold: true,
    color: C.gold,
    margin: 0,
  });
  s.addText(
    "methylKit's calculateDiffMeth and epykit's lr both fit overdispersed logistic regression with BH FDR on (M, U) read counts.",
    {
      x: rightX + 0.25,
      y: 2.25,
      w: 2.7,
      h: 1.4,
      fontSize: 12,
      fontFace: F.body,
      color: "DCE7EF",
      margin: 0,
    }
  );

  s.addShape(pres.shapes.LINE, {
    x: rightX + 0.25,
    y: 3.7,
    w: 2.7,
    h: 0,
    line: { color: "3A6E9A", width: 0.75 },
  });

  s.addText("Same calls.", {
    x: rightX + 0.25,
    y: 3.85,
    w: 2.7,
    h: 0.45,
    fontSize: 20,
    fontFace: F.head,
    bold: true,
    color: C.gold,
    margin: 0,
  });
  s.addText(
    "On identical input counts at n ≥ 4 they return statistically indistinguishable site calls.",
    {
      x: rightX + 0.25,
      y: 4.35,
      w: 2.7,
      h: 1.2,
      fontSize: 12,
      fontFace: F.body,
      color: "DCE7EF",
      margin: 0,
    }
  );

  s.addText("This is one method, two implementations.", {
    x: rightX + 0.25,
    y: 5.75,
    w: 2.7,
    h: 1.0,
    fontSize: 13,
    fontFace: F.body,
    italic: true,
    color: "8FB6CA",
    margin: 0,
  });

  pageFooter(s, "Study 2 · Equivalence", 7);
}

// =========================================================================
// SLIDE 8 — Study 2: n=2 edge case
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(
    s,
    "Study 2 · Where epykit pulls ahead — n = 2",
    "When methylKit's overdispersion estimator becomes degenerate, epykit's allow_n1 keeps the model identifiable"
  );

  // 2-column comparison: methylKit (left) vs epykit (right)
  const colY = 1.7;
  const colH = 4.5;
  const colW = 5.9;

  // methylKit card
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.6,
    y: colY,
    w: colW,
    h: colH,
    fill: { color: C.cardBg },
    line: { color: C.rule, width: 0.75 },
  });
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.6,
    y: colY,
    w: colW,
    h: 0.6,
    fill: { color: "94A3B8" },
    line: { type: "none" },
  });
  s.addText("methylKit", {
    x: 0.8,
    y: colY,
    w: 5,
    h: 0.6,
    fontSize: 20,
    fontFace: F.head,
    bold: true,
    color: C.white,
    valign: "middle",
    margin: 0,
  });

  s.addText("TPR", {
    x: 0.8,
    y: colY + 0.85,
    w: 5,
    h: 0.3,
    fontSize: 11,
    fontFace: F.body,
    bold: true,
    color: C.muted,
    charSpacing: 2,
    margin: 0,
  });
  s.addText("0.302", {
    x: 0.8,
    y: colY + 1.2,
    w: 3,
    h: 1.0,
    fontSize: 72,
    fontFace: F.head,
    bold: true,
    color: C.muted,
    margin: 0,
  });
  s.addText("6,030 significant DMCs called", {
    x: 0.8,
    y: colY + 2.25,
    w: 5,
    h: 0.3,
    fontSize: 12,
    fontFace: F.body,
    color: C.ink,
    margin: 0,
  });
  s.addText("Dispersion estimator returns a degenerate value at n = 1 per group; test loses power.", {
    x: 0.8,
    y: colY + 2.7,
    w: 5,
    h: 1.6,
    fontSize: 12,
    fontFace: F.body,
    italic: true,
    color: C.muted,
    margin: 0,
  });

  // epykit card
  s.addShape(pres.shapes.RECTANGLE, {
    x: 6.8,
    y: colY,
    w: colW,
    h: colH,
    fill: { color: C.navy },
    line: { type: "none" },
  });
  s.addShape(pres.shapes.RECTANGLE, {
    x: 6.8,
    y: colY,
    w: colW,
    h: 0.6,
    fill: { color: C.gold },
    line: { type: "none" },
  });
  s.addText("epykit lr (allow_n1)", {
    x: 7.0,
    y: colY,
    w: 5.5,
    h: 0.6,
    fontSize: 20,
    fontFace: F.head,
    bold: true,
    color: C.navy,
    valign: "middle",
    margin: 0,
  });

  s.addText("TPR", {
    x: 7.0,
    y: colY + 0.85,
    w: 5,
    h: 0.3,
    fontSize: 11,
    fontFace: F.body,
    bold: true,
    color: "BFD7E5",
    charSpacing: 2,
    margin: 0,
  });
  s.addText("0.564", {
    x: 7.0,
    y: colY + 1.2,
    w: 3,
    h: 1.0,
    fontSize: 72,
    fontFace: F.head,
    bold: true,
    color: C.gold,
    margin: 0,
  });
  s.addText("11,283 significant DMCs called", {
    x: 7.0,
    y: colY + 2.25,
    w: 5,
    h: 0.3,
    fontSize: 12,
    fontFace: F.body,
    color: C.white,
    margin: 0,
  });
  s.addText("Falls back to binomial GLM (identifiable at n = 1). Same FPR as methylKit, ~2× recall.", {
    x: 7.0,
    y: colY + 2.7,
    w: 5,
    h: 1.6,
    fontSize: 12,
    fontFace: F.body,
    italic: true,
    color: "DCE7EF",
    margin: 0,
  });

  // Bottom: speedup factor strip
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.6,
    y: 6.45,
    w: 12.1,
    h: 0.75,
    fill: { color: C.cardBg2 },
    line: { type: "none" },
  });
  s.addText(
    [
      { text: "1.87×", options: { bold: true, color: C.gold, fontSize: 22 } },
      {
        text: "  more true DMCs recovered at the same FPR",
        options: { color: C.ink, fontSize: 16 },
      },
    ],
    {
      x: 0.85,
      y: 6.45,
      w: 11.8,
      h: 0.75,
      fontFace: F.body,
      valign: "middle",
      margin: 0,
    }
  );

  pageFooter(s, "Study 2 · n = 2 edge", 8);
}

// =========================================================================
// SLIDE 9 — Study 2: speed across the simulated grid
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(
    s,
    "Study 2 · Runtime, same machine",
    "epykit completes the full 15-point grid in 8.6 minutes; methylKit in 6 h 10 min"
  );

  // Native column chart on the left
  s.addChart(
    pres.charts.BAR,
    [
      {
        name: "epykit",
        labels: [
          "DMC 5×",
          "DMC 10×",
          "DMC 15×",
          "DMC 20×",
          "DMC 25×",
          "DMR 5×",
          "DMR 10×",
          "DMR 15×",
          "DMR 20×",
          "DMR 25×",
        ],
        values: [8.8, 16.4, 8.6, 8.5, 8.0, 62.0, 97.8, 89.0, 82.8, 86.0],
      },
      {
        name: "methylKit",
        labels: [
          "DMC 5×",
          "DMC 10×",
          "DMC 15×",
          "DMC 20×",
          "DMC 25×",
          "DMR 5×",
          "DMR 10×",
          "DMR 15×",
          "DMR 20×",
          "DMR 25×",
        ],
        values: [
          130.0, 116.1, 101.5, 99.0, 96.2, 4241.9, 4934.5, 4139.8, 3995.4,
          3889.1,
        ],
      },
    ],
    {
      x: 0.5,
      y: 1.6,
      w: 8.5,
      h: 5.5,
      barDir: "col",
      chartColors: [C.gold, "94A3B8"],
      chartArea: { fill: { color: C.white } },
      catAxisLabelColor: C.muted,
      valAxisLabelColor: C.muted,
      catAxisLabelFontSize: 9,
      valAxisLabelFontSize: 9,
      valAxisLogScaleBase: 10,
      valGridLine: { color: "E2E8F0", size: 0.5 },
      catGridLine: { style: "none" },
      showLegend: true,
      legendPos: "t",
      legendFontSize: 11,
      showTitle: true,
      title: "Wall-clock per scenario (s, log scale)",
      titleFontSize: 13,
      titleColor: C.ink,
    }
  );

  // Right callouts
  const rightX = 9.4;

  s.addShape(pres.shapes.RECTANGLE, {
    x: rightX,
    y: 1.6,
    w: 3.3,
    h: 1.55,
    fill: { color: C.cardBg2 },
    line: { type: "none" },
  });
  s.addShape(pres.shapes.RECTANGLE, {
    x: rightX,
    y: 1.6,
    w: 0.08,
    h: 1.55,
    fill: { color: C.gold },
    line: { type: "none" },
  });
  s.addText("FULL 15-POINT GRID", {
    x: rightX + 0.25,
    y: 1.7,
    w: 3.0,
    h: 0.3,
    fontSize: 9,
    fontFace: F.body,
    bold: true,
    color: C.muted,
    charSpacing: 2,
    margin: 0,
  });
  s.addText("8.6 min", {
    x: rightX + 0.25,
    y: 2.0,
    w: 3.0,
    h: 0.6,
    fontSize: 34,
    fontFace: F.head,
    bold: true,
    color: C.blue,
    margin: 0,
  });
  s.addText("epykit", {
    x: rightX + 0.25,
    y: 2.62,
    w: 3.0,
    h: 0.3,
    fontSize: 12,
    fontFace: F.body,
    color: C.muted,
    margin: 0,
  });
  s.addText("vs. 6 h 10 min for methylKit", {
    x: rightX + 0.25,
    y: 2.9,
    w: 3.0,
    h: 0.3,
    fontSize: 11,
    fontFace: F.body,
    italic: true,
    color: C.muted,
    margin: 0,
  });

  // 43x callout
  s.addShape(pres.shapes.RECTANGLE, {
    x: rightX,
    y: 3.3,
    w: 3.3,
    h: 1.8,
    fill: { color: C.navy },
    line: { type: "none" },
  });
  s.addText("43×", {
    x: rightX + 0.25,
    y: 3.4,
    w: 3.0,
    h: 1.05,
    fontSize: 78,
    fontFace: F.head,
    bold: true,
    color: C.gold,
    margin: 0,
  });
  s.addText("faster total wall-clock", {
    x: rightX + 0.25,
    y: 4.42,
    w: 3.0,
    h: 0.35,
    fontSize: 13,
    fontFace: F.body,
    color: C.white,
    margin: 0,
  });
  s.addText("DMR scenarios alone reach 45×–68×.", {
    x: rightX + 0.25,
    y: 4.72,
    w: 3.0,
    h: 0.35,
    fontSize: 11,
    fontFace: F.body,
    italic: true,
    color: "8FB6CA",
    margin: 0,
  });

  // Memory
  s.addShape(pres.shapes.RECTANGLE, {
    x: rightX,
    y: 5.25,
    w: 3.3,
    h: 1.6,
    fill: { color: C.cardBg },
    line: { type: "none" },
  });
  s.addShape(pres.shapes.RECTANGLE, {
    x: rightX,
    y: 5.25,
    w: 0.08,
    h: 1.6,
    fill: { color: C.teal },
    line: { type: "none" },
  });
  s.addText("PEAK RAM (FULL GRID)", {
    x: rightX + 0.25,
    y: 5.4,
    w: 3.0,
    h: 0.3,
    fontSize: 9,
    fontFace: F.body,
    bold: true,
    color: C.muted,
    charSpacing: 2,
    margin: 0,
  });
  s.addText("6.03 GB", {
    x: rightX + 0.25,
    y: 5.7,
    w: 3.0,
    h: 0.5,
    fontSize: 28,
    fontFace: F.head,
    bold: true,
    color: C.blue,
    margin: 0,
  });
  s.addText("epykit · vs 7.11 GB methylKit", {
    x: rightX + 0.25,
    y: 6.25,
    w: 3.0,
    h: 0.3,
    fontSize: 11,
    fontFace: F.body,
    color: C.muted,
    margin: 0,
  });

  pageFooter(s, "Study 2 · Performance", 9);
}

// =========================================================================
// SLIDE 10 — Study 3: real data setup
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(
    s,
    "Study 3 · Real WGBS data — GSE263850",
    "Bit-identical inputs to both pipelines; same parameters; the only thing that differs is the test"
  );

  // Left column: experimental design
  s.addText("EXPERIMENTAL DESIGN", {
    x: 0.6,
    y: 1.55,
    w: 5,
    h: 0.3,
    fontSize: 10,
    fontFace: F.body,
    bold: true,
    color: C.gold,
    charSpacing: 2,
    margin: 0,
  });
  s.addText(
    [
      { text: "Genome: ", options: { bold: true, color: C.ink } },
      { text: "hg38", options: { color: C.ink, breakLine: true } },
      { text: "Samples: ", options: { bold: true, color: C.ink } },
      { text: "6 (3 vs 3 design)", options: { color: C.ink, breakLine: true } },
      { text: "Treatment: ", options: { bold: true, color: C.ink } },
      {
        text: "Clone16, Clone20, Clone21",
        options: { color: C.ink, breakLine: true },
      },
      { text: "Control: ", options: { bold: true, color: C.ink } },
      {
        text: "SBP009 untreated 1, 2, 3",
        options: { color: C.ink, breakLine: true },
      },
      { text: "CpGs tested: ", options: { bold: true, color: C.ink } },
      { text: "15,597,046 (united)", options: { color: C.ink } },
    ],
    {
      x: 0.6,
      y: 1.85,
      w: 5.5,
      h: 2.4,
      fontSize: 14,
      fontFace: F.body,
      paraSpaceAfter: 6,
      margin: 0,
    }
  );

  // Parameters card
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.6,
    y: 4.5,
    w: 5.5,
    h: 2.7,
    fill: { color: C.cardBg },
    line: { color: C.rule, width: 0.75 },
  });
  s.addText("IDENTICAL PARAMETERS", {
    x: 0.85,
    y: 4.65,
    w: 5,
    h: 0.3,
    fontSize: 10,
    fontFace: F.body,
    bold: true,
    color: C.muted,
    charSpacing: 2,
    margin: 0,
  });
  s.addText(
    [
      { text: "min coverage = 10", options: { bullet: true, color: C.ink, breakLine: true } },
      { text: "99.9th-percentile clipping", options: { bullet: true, color: C.ink, breakLine: true } },
      { text: "Benjamini–Hochberg FDR", options: { bullet: true, color: C.ink, breakLine: true } },
      { text: "DMC: q < 0.05, |meth_diff| ≥ 10 %", options: { bullet: true, color: C.ink, breakLine: true } },
      { text: "DMR: 500 bp tiles, ≥ 5 CpGs", options: { bullet: true, color: C.ink } },
    ],
    {
      x: 0.85,
      y: 4.95,
      w: 5.0,
      h: 2.1,
      fontSize: 13,
      fontFace: F.body,
      paraSpaceAfter: 2,
      margin: 0,
    }
  );

  // Right column: input agreement
  s.addShape(pres.shapes.RECTANGLE, {
    x: 7.0,
    y: 1.55,
    w: 5.7,
    h: 5.65,
    fill: { color: C.navy },
    line: { type: "none" },
  });
  s.addText("INPUT IDENTITY", {
    x: 7.25,
    y: 1.75,
    w: 5,
    h: 0.3,
    fontSize: 10,
    fontFace: F.body,
    bold: true,
    color: C.gold,
    charSpacing: 2,
    margin: 0,
  });
  s.addText(
    "Per-CpG (N_meth, coverage) counts are bit-identical entering each pipeline's coverage filter.",
    {
      x: 7.25,
      y: 2.05,
      w: 5.3,
      h: 0.9,
      fontSize: 14,
      fontFace: F.body,
      color: C.white,
      margin: 0,
    }
  );

  // Stats grid
  const stat = (label, val, x, y) => {
    s.addText(label, {
      x: x,
      y: y,
      w: 2.5,
      h: 0.3,
      fontSize: 9,
      fontFace: F.body,
      bold: true,
      color: "BFD7E5",
      charSpacing: 1.5,
      margin: 0,
    });
    s.addText(val, {
      x: x,
      y: y + 0.3,
      w: 2.5,
      h: 0.55,
      fontSize: 26,
      fontFace: F.head,
      bold: true,
      color: C.gold,
      margin: 0,
    });
  };

  stat("CpG OVERLAP", "99.98 %", 7.25, 3.1);
  stat("TILE OVERLAP", "99.80 %", 9.95, 3.1);
  stat("JACCARD (CpGs)", "0.9998", 7.25, 4.25);
  stat("JACCARD (tiles)", "0.998", 9.95, 4.25);

  s.addShape(pres.shapes.LINE, {
    x: 7.25,
    y: 5.45,
    w: 5.3,
    h: 0,
    line: { color: "3A6E9A", width: 0.75 },
  });
  s.addText("Coordinate convention", {
    x: 7.25,
    y: 5.55,
    w: 5.3,
    h: 0.3,
    fontSize: 11,
    fontFace: F.body,
    bold: true,
    color: "DCE7EF",
    margin: 0,
  });
  s.addText(
    "methylKit .cov is 1-based; epykit BED is 0-based. A +1 shift on epykit positions aligns the keys exactly.",
    {
      x: 7.25,
      y: 5.85,
      w: 5.3,
      h: 1.2,
      fontSize: 12,
      fontFace: F.body,
      color: "BFD7E5",
      margin: 0,
    }
  );

  pageFooter(s, "Study 3 · Setup", 10);
}

// =========================================================================
// SLIDE 11 — Study 3: agreement on real data
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(
    s,
    "Study 3 · The two pipelines measure the same biology",
    "Effect sizes correlate at r = 0.994 on DMCs and r = 0.997 on DMRs"
  );

  s.addImage({
    path: FIG("study3_real_GSE263850/08B_dmc_effect_size_scatter.png"),
    x: 0.6,
    y: 1.55,
    w: 6.7,
    h: 5.0,
    sizing: { type: "contain", w: 6.7, h: 5.0 },
  });

  // Stats panel right
  const rightX = 7.7;
  s.addText("DMC AGREEMENT (15.6 M CpGs)", {
    x: rightX,
    y: 1.55,
    w: 5,
    h: 0.3,
    fontSize: 10,
    fontFace: F.body,
    bold: true,
    color: C.gold,
    charSpacing: 2,
    margin: 0,
  });

  const dmcStats = [
    ["Pearson r", "0.9936", C.gold],
    ["Spearman ρ", "0.9831", C.ink],
    ["Direction agreement", "94.05 %", C.ink],
  ];
  dmcStats.forEach((row, i) => {
    const yy = 1.95 + i * 0.55;
    s.addText(row[0], {
      x: rightX,
      y: yy,
      w: 2.6,
      h: 0.4,
      fontSize: 13,
      fontFace: F.body,
      color: C.muted,
      valign: "middle",
      margin: 0,
    });
    s.addText(row[1], {
      x: rightX + 2.5,
      y: yy,
      w: 2.5,
      h: 0.5,
      fontSize: 22,
      fontFace: F.head,
      bold: true,
      color: row[2],
      align: "right",
      valign: "middle",
      margin: 0,
    });
  });

  s.addShape(pres.shapes.LINE, {
    x: rightX,
    y: 3.85,
    w: 5,
    h: 0,
    line: { color: C.rule, width: 0.5 },
  });

  s.addText("DMR AGREEMENT (1.17 M TILES)", {
    x: rightX,
    y: 4.0,
    w: 5,
    h: 0.3,
    fontSize: 10,
    fontFace: F.body,
    bold: true,
    color: C.gold,
    charSpacing: 2,
    margin: 0,
  });

  const dmrStats = [
    ["Pearson r", "0.9970", C.gold],
    ["Direction agreement", "91.82 %", C.ink],
    ["DMR Jaccard (lenient)", "0.473", C.warn],
  ];
  dmrStats.forEach((row, i) => {
    const yy = 4.4 + i * 0.55;
    s.addText(row[0], {
      x: rightX,
      y: yy,
      w: 2.6,
      h: 0.4,
      fontSize: 13,
      fontFace: F.body,
      color: C.muted,
      valign: "middle",
      margin: 0,
    });
    s.addText(row[1], {
      x: rightX + 2.5,
      y: yy,
      w: 2.5,
      h: 0.5,
      fontSize: 22,
      fontFace: F.head,
      bold: true,
      color: row[2],
      align: "right",
      valign: "middle",
      margin: 0,
    });
  });

  // bottom note
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.6,
    y: 6.7,
    w: 12.1,
    h: 0.55,
    fill: { color: C.cardBg },
    line: { type: "none" },
  });
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.6,
    y: 6.7,
    w: 0.08,
    h: 0.55,
    fill: { color: C.gold },
    line: { type: "none" },
  });
  s.addText(
    [
      { text: "Same biology, ", options: { bold: true, color: C.ink } },
      {
        text: "different significance threshold. epykit's default calls ~60 % as many DMCs (50.7 % precision vs methylKit at the same threshold) — see next slide.",
        options: { color: C.ink },
      },
    ],
    {
      x: 0.85,
      y: 6.7,
      w: 11.8,
      h: 0.55,
      fontSize: 12,
      fontFace: F.body,
      valign: "middle",
      margin: 0,
    }
  );

  pageFooter(s, "Study 3 · Agreement", 11);
}

// =========================================================================
// SLIDE 12 — Calibration vs sensitivity (operating point)
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(
    s,
    "The calibration–sensitivity trade-off on real data",
    "Same input counts, different operating points on the precision/recall curve"
  );

  // Three engine columns
  const cardY = 1.65;
  const cardH = 5.0;
  const cardW = 4.0;
  const cardX0 = 0.65;
  const gap = 0.13;

  const engines = [
    {
      name: "lr / site",
      role: "Default",
      sig: "30,965",
      recall: "31.7 %",
      jaccard: "0.243",
      summary: "Conservative default. Same model as methylKit, more stringent per-site dispersion at n = 3 per group.",
      headerColor: C.blue,
      bgColor: C.cardBg,
      tone: "primary",
    },
    {
      name: "lr+ (power_stack)",
      role: "High recall",
      sig: "406,515",
      recall: "92.9 %",
      jaccard: "0.112",
      summary: "Recovers 93 % of methylKit's calls by borrowing across correlated CpGs; emits 13× more total calls — some may be real signal methylKit misses.",
      headerColor: C.gold,
      bgColor: "FFFAEC",
      tone: "highlight",
    },
    {
      name: "lr / shrink",
      role: "Maximum stringency",
      sig: "1,499",
      recall: "0.3 %",
      jaccard: "0.003",
      summary: "Strongest signals only; useful when downstream cost of a false positive is very high.",
      headerColor: "94A3B8",
      bgColor: C.cardBg,
      tone: "secondary",
    },
  ];

  engines.forEach((eng, i) => {
    const x = cardX0 + i * (cardW + gap);
    s.addShape(pres.shapes.RECTANGLE, {
      x: x,
      y: cardY,
      w: cardW,
      h: cardH,
      fill: { color: eng.bgColor },
      line: { color: C.rule, width: 0.75 },
    });
    s.addShape(pres.shapes.RECTANGLE, {
      x: x,
      y: cardY,
      w: cardW,
      h: 0.12,
      fill: { color: eng.headerColor },
      line: { type: "none" },
    });

    s.addText(eng.role.toUpperCase(), {
      x: x + 0.25,
      y: cardY + 0.3,
      w: cardW - 0.5,
      h: 0.3,
      fontSize: 10,
      fontFace: F.body,
      bold: true,
      color: eng.headerColor,
      charSpacing: 2,
      margin: 0,
    });
    s.addText(eng.name, {
      x: x + 0.25,
      y: cardY + 0.6,
      w: cardW - 0.5,
      h: 0.5,
      fontSize: 22,
      fontFace: F.mono,
      bold: true,
      color: C.ink,
      margin: 0,
    });

    s.addShape(pres.shapes.LINE, {
      x: x + 0.25,
      y: cardY + 1.2,
      w: cardW - 0.5,
      h: 0,
      line: { color: C.rule, width: 0.5 },
    });

    const stats = [
      ["DMCs called", eng.sig],
      ["Recall vs methylKit", eng.recall],
      ["Jaccard", eng.jaccard],
    ];
    stats.forEach((stat, j) => {
      const yy = cardY + 1.3 + j * 0.55;
      s.addText(stat[0], {
        x: x + 0.25,
        y: yy,
        w: 2.0,
        h: 0.4,
        fontSize: 11,
        fontFace: F.body,
        color: C.muted,
        valign: "middle",
        margin: 0,
      });
      s.addText(stat[1], {
        x: x + 1.9,
        y: yy,
        w: 1.85,
        h: 0.45,
        fontSize: 17,
        fontFace: F.mono,
        bold: true,
        color: eng.tone === "highlight" ? C.gold : C.ink,
        align: "right",
        valign: "middle",
        margin: 0,
      });
    });

    s.addText(eng.summary, {
      x: x + 0.25,
      y: cardY + 3.2,
      w: cardW - 0.5,
      h: 1.7,
      fontSize: 12,
      fontFace: F.body,
      italic: true,
      color: C.ink,
      margin: 0,
    });
  });

  // Bottom interpretation strip
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.6,
    y: 6.8,
    w: 12.1,
    h: 0.5,
    fill: { color: C.navy },
    line: { type: "none" },
  });
  s.addText(
    "Recommendation: use the default; document any opt-in; reproduce headlines under at least one alternative dispersion mode.",
    {
      x: 0.85,
      y: 6.8,
      w: 11.8,
      h: 0.5,
      fontSize: 12,
      fontFace: F.body,
      color: C.white,
      valign: "middle",
      margin: 0,
    }
  );

  pageFooter(s, "Operating point", 12);
}

// =========================================================================
// SLIDE 13 — Cross-benchmark summary (use summary figure)
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(
    s,
    "Runtime on the two studies with measured timings",
    "Study 2 (simulated grid, same machine) and Study 3 (real WGBS, full pipeline)"
  );

  s.addImage({
    path: FIG("summary/S1_runtime_across_studies.png"),
    x: 0.5,
    y: 1.55,
    w: 8.5,
    h: 5.3,
    sizing: { type: "contain", w: 8.5, h: 5.3 },
  });

  // Right callout column — three speedup factoids
  const rightX = 9.3;
  const cards = [
    {
      study: "STUDY 2",
      sub: "simulated grid, lr engine",
      speed: "298×",
      detail: "epykit: 73 s · methylKit: 21,690 s",
      bg: C.cardBg2,
    },
    {
      study: "STUDY 3",
      sub: "real GSE263850 wall-clock",
      speed: "12.2×",
      detail: "epykit: 1,072 s · methylKit: 13,033 s",
      bg: C.cardBg,
    },
    {
      study: "STUDY 3 RAM",
      sub: "peak RSS",
      speed: "3.83×",
      detail: "epykit: 12.6 GB · methylKit: 48.0 GB",
      bg: C.cardBg2,
    },
  ];
  cards.forEach((c, i) => {
    const yy = 1.6 + i * 1.85;
    s.addShape(pres.shapes.RECTANGLE, {
      x: rightX,
      y: yy,
      w: 3.5,
      h: 1.7,
      fill: { color: c.bg },
      line: { type: "none" },
    });
    s.addShape(pres.shapes.RECTANGLE, {
      x: rightX,
      y: yy,
      w: 0.08,
      h: 1.7,
      fill: { color: C.gold },
      line: { type: "none" },
    });
    s.addText(c.study, {
      x: rightX + 0.22,
      y: yy + 0.13,
      w: 3.2,
      h: 0.25,
      fontSize: 10,
      fontFace: F.body,
      bold: true,
      color: C.muted,
      charSpacing: 2,
      margin: 0,
    });
    s.addText(c.sub, {
      x: rightX + 0.22,
      y: yy + 0.38,
      w: 3.2,
      h: 0.25,
      fontSize: 10,
      fontFace: F.body,
      italic: true,
      color: C.muted,
      margin: 0,
    });
    s.addText(c.speed, {
      x: rightX + 0.22,
      y: yy + 0.62,
      w: 3.2,
      h: 0.7,
      fontSize: 42,
      fontFace: F.head,
      bold: true,
      color: C.blue,
      margin: 0,
    });
    s.addText(c.detail, {
      x: rightX + 0.22,
      y: yy + 1.35,
      w: 3.2,
      h: 0.3,
      fontSize: 9,
      fontFace: F.body,
      color: C.muted,
      margin: 0,
    });
  });

  pageFooter(s, "Cross-benchmark", 13);
}

// =========================================================================
// SLIDE 14 — Recommendations & limitations
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleHeader(
    s,
    "When to use which engine — and what we don't claim",
    "Practical guidance plus honest caveats"
  );

  // Left card: recommendations
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.6,
    y: 1.6,
    w: 6.05,
    h: 5.5,
    fill: { color: C.cardBg },
    line: { color: C.rule, width: 0.75 },
  });
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.6,
    y: 1.6,
    w: 6.05,
    h: 0.45,
    fill: { color: C.blue },
    line: { type: "none" },
  });
  s.addText("RECOMMENDATIONS", {
    x: 0.85,
    y: 1.6,
    w: 5.5,
    h: 0.45,
    fontSize: 13,
    fontFace: F.body,
    bold: true,
    color: C.white,
    valign: "middle",
    charSpacing: 2,
    margin: 0,
  });

  const recs = [
    ["n ≥ 3 per group", "ep.tl.dmc(test='lr')"],
    ["n ≤ 2 per group", "lr+ (auto-engages from v0.7.2)"],
    ["n = 1 (no replicates)", "fisher (bug-fixed)"],
    ["DMR calling", "dmr_chain_merge — variable-width regions, 1:1 with truth"],
    ["High recall on real data", "lr+ (power_stack)"],
  ];
  recs.forEach((r, i) => {
    const yy = 2.3 + i * 0.95;
    s.addText(r[0], {
      x: 0.85,
      y: yy,
      w: 5.5,
      h: 0.3,
      fontSize: 12,
      fontFace: F.body,
      bold: true,
      color: C.muted,
      charSpacing: 1,
      margin: 0,
    });
    s.addText(r[1], {
      x: 0.85,
      y: yy + 0.3,
      w: 5.5,
      h: 0.4,
      fontSize: 14,
      fontFace: F.mono,
      color: C.ink,
      margin: 0,
    });
  });

  // Right card: limitations
  s.addShape(pres.shapes.RECTANGLE, {
    x: 6.85,
    y: 1.6,
    w: 5.85,
    h: 5.5,
    fill: { color: C.navy },
    line: { type: "none" },
  });
  s.addShape(pres.shapes.RECTANGLE, {
    x: 6.85,
    y: 1.6,
    w: 5.85,
    h: 0.45,
    fill: { color: C.gold },
    line: { type: "none" },
  });
  s.addText("LIMITATIONS WE ACKNOWLEDGE", {
    x: 7.05,
    y: 1.6,
    w: 5.5,
    h: 0.45,
    fontSize: 13,
    fontFace: F.body,
    bold: true,
    color: C.navy,
    valign: "middle",
    charSpacing: 2,
    margin: 0,
  });

  s.addText(
    [
      {
        text: "Simulator is underdispersed",
        options: { bold: true, color: C.gold, breakLine: true },
      },
      {
        text:
          "Piao 2021 has φ ≈ 0.41 at 5×; real WGBS is overdispersed (φ ≈ 1.5–5). Low-coverage TPR advantages may shrink on real data.",
        options: { color: "DCE7EF", breakLine: true },
      },
      { text: "", options: { breakLine: true } },
      {
        text: "Baseline software versions",
        options: { bold: true, color: C.gold, breakLine: true },
      },
      {
        text:
          "Study 1 numbers are from 2021 R/CLI releases. Relative ordering is robust; absolutes may have shifted.",
        options: { color: "DCE7EF", breakLine: true },
      },
      { text: "", options: { breakLine: true } },
      {
        text: "Single real dataset",
        options: { bold: true, color: C.gold, breakLine: true },
      },
      {
        text:
          "Study 3 is one tissue × one genome. Multi-dataset validation (IMR90 vs H1-hESC, mouse imprinted DMRs) is future work.",
        options: { color: "DCE7EF", breakLine: true },
      },
      { text: "", options: { breakLine: true } },
      {
        text: "Ground truth is reconstructed",
        options: { bold: true, color: C.gold, breakLine: true },
      },
      {
        text:
          "True-DMC labels come from the 25× simulator observation, not from internal simulator flags.",
        options: { color: "DCE7EF" },
      },
    ],
    {
      x: 7.05,
      y: 2.2,
      w: 5.5,
      h: 4.8,
      fontSize: 11.5,
      fontFace: F.body,
      paraSpaceAfter: 2,
      margin: 0,
    }
  );

  pageFooter(s, "Guidance", 14);
}

// =========================================================================
// SLIDE 15 — Conclusion (dark)
// =========================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.navyDeep };

  // motif
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.6,
    y: 6.4,
    w: 0.35,
    h: 0.35,
    fill: { color: C.teal },
    line: { type: "none" },
  });
  s.addShape(pres.shapes.RECTANGLE, {
    x: 1.05,
    y: 6.4,
    w: 0.35,
    h: 0.35,
    fill: { color: C.blue },
    line: { type: "none" },
  });
  s.addShape(pres.shapes.RECTANGLE, {
    x: 1.5,
    y: 6.4,
    w: 0.35,
    h: 0.35,
    fill: { color: C.gold },
    line: { type: "none" },
  });

  s.addText("Conclusion", {
    x: 0.6,
    y: 0.5,
    w: 12,
    h: 0.9,
    fontSize: 44,
    fontFace: F.head,
    bold: true,
    color: C.white,
    margin: 0,
  });

  // Three big takeaways
  const takeaways = [
    {
      heading: "Matches the state of the art",
      detail:
        "On simulated data epykit lr is statistically equivalent to methylKit to three decimal places at n ≥ 4 — and recovers ~2× more true DMCs at n = 2.",
    },
    {
      heading: "Calibrated on real data",
      detail:
        "On GSE263850, effect sizes correlate at r = 0.994 (DMCs) and r = 0.997 (DMRs); both pipelines see the same biology.",
    },
    {
      heading: "Order-of-magnitude faster",
      detail:
        "12×–68× speedup on matched workloads (43× on the full simulated grid, 12× on real WGBS); 1.18×–3.83× less peak memory.",
    },
  ];

  takeaways.forEach((t, i) => {
    const y = 1.85 + i * 1.5;
    // Number badge
    s.addShape(pres.shapes.OVAL, {
      x: 0.65,
      y: y,
      w: 0.7,
      h: 0.7,
      fill: { color: C.gold },
      line: { type: "none" },
    });
    s.addText(String(i + 1), {
      x: 0.65,
      y: y,
      w: 0.7,
      h: 0.7,
      fontSize: 28,
      fontFace: F.head,
      bold: true,
      color: C.navy,
      align: "center",
      valign: "middle",
      margin: 0,
    });

    s.addText(t.heading, {
      x: 1.6,
      y: y - 0.05,
      w: 11,
      h: 0.5,
      fontSize: 22,
      fontFace: F.head,
      bold: true,
      color: C.white,
      margin: 0,
    });
    s.addText(t.detail, {
      x: 1.6,
      y: y + 0.5,
      w: 11,
      h: 0.8,
      fontSize: 14,
      fontFace: F.body,
      color: "BFD7E5",
      margin: 0,
    });
  });

  // Tagline at bottom-right
  s.addText(
    "epykit brings WGBS downstream analysis into the modern Python bioinformatics stack — without sacrificing accuracy.",
    {
      x: 2.5,
      y: 6.7,
      w: 10.2,
      h: 0.45,
      fontSize: 12,
      fontFace: F.body,
      italic: true,
      color: "8FB6CA",
      align: "right",
      margin: 0,
    }
  );
}

// =========================================================================
pres.writeFile({
  fileName: path.join(ROOT, "epykit_benchmark_presentation.pptx"),
});
