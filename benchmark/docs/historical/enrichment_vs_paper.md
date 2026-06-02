# Pathway enrichment — epykit DMRs vs the AKAP11 paper

Submitted gene sets to **Enrichr** (Maayan lab) and pulled the top 10 enriched terms from **Reactome 2022** (vs paper panel D) and **GO Molecular Function 2023** (vs paper panel E). Paper-reported headline terms are flagged with a ★ when a related term appears in our top 10.

## Paper's reported pathways (target)

**Panel D (Reactome):** Class A/1 Rhodopsin-like receptors · Peptide ligand-binding receptors · GPCR ligand binding · GPCR downstream signalling · Signalling by GPCR · G alpha i signalling events

**Panel E (GO MF):** Sequence-specific DNA binding · DNA-binding TF activity · RNA Pol II cis-regulatory region DNA binding · TF binding · DNA-binding TF activator activity


## Input gene set: `epykit_top500_DMRs` (612 genes)

### Panel D (Reactome) — top 15 (★ = matches paper)

| Rank | Term | adj p | n / overlap | Sample overlap genes |
|---:|---|---:|---:|---|
| 1 | MECP2 Regulates Transcription Factors R-HSA-9022707 | 1.00e+00 | 2 | MEF2C, RBFOX1 |
| 2 | Reelin Signaling Pathway R-HSA-8866376 | 1.00e+00 | 2 | RELN, FYN |
| 3 | L1CAM Interactions R-HSA-373760 | 1.00e+00 | 8 | NRP2, RPS6KA2, CNTN1, LAMB1, SCN5A, AP2B1, SPTBN1, ITGA9 |
| 4 | HS-GAG Biosynthesis R-HSA-2022928 | 1.00e+00 | 4 | NDST4, HS3ST5, HS3ST4, HS3ST1 |
| 5 | FGFR3b Ligand Binding And Activation R-HSA-190371 | 1.00e+00 | 2 | FGF18, FGF1 |
| 6 | Transport Of RCbl Within Body R-HSA-9758890 | 1.00e+00 | 2 | LMBRD1, LRP2 |
| 7 | Cobalamin (Cbl, Vitamin B12) Transport And Metabolism R-HSA-196741 | 1.00e+00 | 3 | CUBN, LMBRD1, LRP2 |
| 8 | Nephrin Family Interactions R-HSA-373753 | 1.00e+00 | 3 | NPHS1, FYN, SPTBN1 |
| 9 | Activation Of Ca-permeable Kainate Receptor R-HSA-451308 | 1.00e+00 | 2 | GRIK3, GRIK4 |
| 10 | Cell-Cell Communication R-HSA-1500931 | 1.00e+00 | 8 | SDK1, NPHS1, PXDN, CTNNB1, CLDN17, FYN, CDH8, SPTBN1 |
| 11 | Constitutive Signaling By Aberrant PI3K In Cancer R-HSA-2219530 | 1.00e+00 | 6 | PDGFRA, FGF18, IRS2, FYN, FGF1, ICOS |
| 12 | Sodium/Calcium Exchangers R-HSA-425561 | 1.00e+00 | 2 | SLC24A4, SLC8A1 |
| 13 | Uptake Of Dietary Cobalamins Into Enterocytes R-HSA-9758881 | 1.00e+00 | 2 | CUBN, LMBRD1 |
| 14 | PI3K/AKT Signaling In Cancer R-HSA-2219528 | 1.00e+00 | 7 | PDGFRA, MAPKAP1, FGF18, IRS2, FYN, FGF1, ICOS |
| 15 | Platelet Calcium Homeostasis R-HSA-418360 | 1.00e+00 | 3 | TRPC7, P2RX4, SLC8A1 |

### Panel E (GO MF) — top 15 (★ = matches paper)

| Rank | Term | adj p | n / overlap | Sample overlap genes |
|---:|---|---:|---:|---|
| 1 | Calcium:Sodium Antiporter Activity (GO:0005432) | 7.01e-01 | 3 | SLC24A2, SLC24A4, SLC8A1 |
| 2 | Heparan Sulfate Sulfotransferase Activity (GO:0034483) | 7.01e-01 | 3 | NDST4, HS3ST4, HS3ST1 |
| 3 | CXCR Chemokine Receptor Binding (GO:0045236) | 7.01e-01 | 3 | CXCL11, CXCL13, CXCL2 |
| 4 | Vascular Endothelial Growth Factor Receptor Activity (GO:0005021) | 7.01e-01 | 2 | NELL1, PDGFRA |
| 5 | Beta-Galactoside (CMP) Alpha-2,3-Sialyltransferase Activity (GO:0003836) | 7.01e-01 | 2 | ST3GAL6, ST3GAL1 |
| 6 | Voltage-Gated Calcium Channel Activity (GO:0005245) | 7.01e-01 | 5 | CACNA2D2, CACNA1D, SCN5A, CACNA2D4, CACNG3 |
| 7 | Ankyrin Binding (GO:0030506) | 7.01e-01 | 3 | SCN5A, SLC8A1, SPTBN1 |
| 8 | Calcium Channel Activity (GO:0005262) | 7.01e-01 | 7 | SLC24A2, SLC24A4, P2RX4, CACNA2D2, CACNA1D, CACNA2D4, CACNG3 |
| 9 | Sialyltransferase Activity (GO:0008373) | 7.01e-01 | 3 | ST6GAL2, ST3GAL6, ST3GAL1 |
| 10 | Ligand-Gated Monoatomic Cation Channel Activity (GO:0099094) | 7.01e-01 | 6 | KCNJ6, P2RX4, GRIK3, GRIK4, HTR3A, KCNN2 |
| 11 | Thioesterase Binding (GO:0031996) | 7.01e-01 | 2 | TRAF3, TRAF2 |
| 12 | Tumor Necrosis Factor Receptor Binding (GO:0005164) | 7.01e-01 | 2 | TRAF3, TRAF2 |
| 13 | Histone H3K36 Demethylase Activity (GO:0051864) | 7.01e-01 | 2 | KDM2B, RIOX2 |
| 14 | Alpha-Actinin Binding (GO:0051393) | 7.01e-01 | 3 | PALLD, CACNA1D, KCNN2 |
| 15 | Transmitter-Gated Monoatomic Ion Channel Activity (GO:0022824) | 7.01e-01 | 4 | GRIK3, GRIK4, HTR3A, GABRG3 |


## Input gene set: `epykit_top813_DMRs` (767 genes)

### Panel D (Reactome) — top 15 (★ = matches paper)

| Rank | Term | adj p | n / overlap | Sample overlap genes |
|---:|---|---:|---:|---|
| 1 | L1CAM Interactions R-HSA-373760 | 3.91e-01 | 12 | NRP1, NRP2, ALCAM, RPS6KA2, SCN8A, CNTN1, LAMB1, SCN5A (+4) |
| 2 | MECP2 Regulates Transcription Factors R-HSA-9022707 | 1.00e+00 | 2 | MEF2C, RBFOX1 |
| 3 | Reelin Signaling Pathway R-HSA-8866376 | 1.00e+00 | 2 | RELN, FYN |
| 4 | Regulation Of IFNG Signaling R-HSA-877312 | 1.00e+00 | 3 | SOCS1, IFNGR1, PTPN2 |
| 5 | SEMA3A-Plexin Repulsion Signaling By Inhibiting Integrin Adhesion R-HSA-399955 | 1.00e+00 | 3 | NRP1, PLXNA2, FYN |
| 6 | CRMPs In Sema3A Signaling R-HSA-399956 | 1.00e+00 | 3 | NRP1, PLXNA2, FYN |
| 7 | Gap Junction Assembly R-HSA-190861 | 1.00e+00 | 3 | GJD2, GJD4, GJB6 |
| 8 | Sema3A PAK Dependent Axon Repulsion R-HSA-399954 | 1.00e+00 | 3 | NRP1, PLXNA2, FYN |
| 9 | Interaction Between L1 And Ankyrins R-HSA-445095 | 1.00e+00 | 4 | SCN8A, SCN5A, ANK1, SPTBN1 |
| 10 | HS-GAG Biosynthesis R-HSA-2022928 | 1.00e+00 | 4 | NDST4, HS3ST5, HS3ST4, HS3ST1 |
| 11 | FGFR3b Ligand Binding And Activation R-HSA-190371 | 1.00e+00 | 2 | FGF18, FGF1 |
| 12 | Inactivation Of CDC42 And RAC1 R-HSA-428543 | 1.00e+00 | 2 | SLIT2, SRGAP1 |
| 13 | Neurofascin Interactions R-HSA-447043 | 1.00e+00 | 2 | CNTN1, ANK1 |
| 14 | NrCAM Interactions R-HSA-447038 | 1.00e+00 | 2 | NRP2, ANK1 |
| 15 | Transport Of RCbl Within Body R-HSA-9758890 | 1.00e+00 | 2 | LMBRD1, LRP2 |

### Panel E (GO MF) — top 15 (★ = matches paper)

| Rank | Term | adj p | n / overlap | Sample overlap genes |
|---:|---|---:|---:|---|
| 1 | Calcium:Sodium Antiporter Activity (GO:0005432) | 8.15e-01 | 3 | SLC24A2, SLC24A4, SLC8A1 |
| 2 | Sialyltransferase Activity (GO:0008373) | 8.15e-01 | 4 | ST6GAL2, ST8SIA3, ST3GAL6, ST3GAL1 |
| 3 | Voltage-Gated Calcium Channel Activity (GO:0005245) | 8.15e-01 | 6 | SCN8A, CACNA2D2, CACNA1D, SCN5A, CACNA2D4, CACNG3 |
| 4 | Heparan Sulfate Sulfotransferase Activity (GO:0034483) | 8.15e-01 | 3 | NDST4, HS3ST4, HS3ST1 |
| 5 | Sodium Channel Activity (GO:0005272) | 8.15e-01 | 5 | SCN8A, GRIK3, GRIK4, SCN5A, ASIC2 |
| 6 | Carboxypeptidase Activity (GO:0004180) | 8.15e-01 | 5 | CPB2, CPB1, CPXM2, BLMH, CPZ |
| 7 | Vascular Endothelial Growth Factor Receptor Activity (GO:0005021) | 8.15e-01 | 2 | NELL1, PDGFRA |
| 8 | Aminophospholipid Flippase Activity (GO:0015247) | 8.15e-01 | 2 | ATP8A2, ATP8A1 |
| 9 | Beta-Galactoside (CMP) Alpha-2,3-Sialyltransferase Activity (GO:0003836) | 8.15e-01 | 2 | ST3GAL6, ST3GAL1 |
| 10 | CXCR Chemokine Receptor Binding (GO:0045236) | 8.15e-01 | 3 | CXCL11, CXCL13, CXCL2 |
| 11 | Metallocarboxypeptidase Activity (GO:0004181) | 8.15e-01 | 4 | CPB2, CPB1, CPXM2, CPZ |
| 12 | Ligand-Gated Monoatomic Cation Channel Activity (GO:0099094) | 8.15e-01 | 7 | KCNJ6, P2RX4, GRIK3, GRIK4, HTR3A, KCNN2, ASIC2 |
| 13 | Ankyrin Binding (GO:0030506) | 8.15e-01 | 3 | SCN5A, SLC8A1, SPTBN1 |
| 14 | Ligand-Gated Calcium Channel Activity (GO:0099604) | 8.15e-01 | 3 | P2RX4, RASA3, RYR3 |
| 15 | Gap Junction Channel Activity (GO:0005243) | 8.15e-01 | 3 | GJD2, GJD4, GJB6 |


## Input gene set: `epykit_top1000_DMCs_genes` (919 genes)

### Panel D (Reactome) — top 15 (★ = matches paper)

| Rank | Term | adj p | n / overlap | Sample overlap genes |
|---:|---|---:|---:|---|
| 1 | Netrin Mediated Repulsion Signals R-HSA-418886 | 1.00e+00 | 3 | DCC, UNC5C, UNC5D |
| 2 | PI Metabolism R-HSA-1483255 | 1.00e+00 | 9 | TPTE2, INPP4A, VAC14, TNFAIP8, PNPLA7, SACM1L, CLSTN1, TNFAIP8L3 (+1) |
| 3 | Activated NTRK2 Signals Thru CDK5 R-HSA-9032845 | 1.00e+00 | 2 | NTRK2, TIAM1 |
| 4 | Activated NTRK2 Signals Thru FYN R-HSA-9032500 | 1.00e+00 | 2 | NTRK2, GRIN2B |
| 5 | WNT Mediated Activation Of DVL R-HSA-201688 | 1.00e+00 | 2 | CSNK2A2, PIP5K1B |
| 6 | Gap Junction Assembly R-HSA-190861 | 1.00e+00 | 3 | GJD2, GJD4, GJB6 |
| 7 | Synthesis Of PIPs At Golgi Membrane R-HSA-1660514 | 1.00e+00 | 3 | TPTE2, VAC14, SACM1L |
| 8 | Activation Of Kainate Receptors Upon Glutamate Binding R-HSA-451326 | 1.00e+00 | 4 | GNG10, DLG1, GRIK4, GNG11 |
| 9 | RORA Activates Gene Expression R-HSA-1368082 | 1.00e+00 | 3 | CREBBP, CARM1, RORA |
| 10 | Ligand-receptor Interactions R-HSA-5632681 | 1.00e+00 | 2 | SHH, CDON |
| 11 | Phospholipid Metabolism R-HSA-1483257 | 1.00e+00 | 15 | TPTE2, VAC14, TNFAIP8, CHKA, PNPLA7, OSBPL5, SACM1L, CSNK2A2 (+7) |
| 12 | Activation Of Ca-permeable Kainate Receptor R-HSA-451308 | 1.00e+00 | 2 | DLG1, GRIK4 |
| 13 | BH3-only Proteins Associate With And Inactivate Anti-Apoptotic BCL-2 Members R-HSA-111453 | 1.00e+00 | 2 | BCL2L11, BCL2 |
| 14 | Beta-oxidation Of pristanoyl-CoA R-HSA-389887 | 1.00e+00 | 2 | ACOXL, ACOX3 |
| 15 | Regulation Of Commissural Axon Pathfinding By SLIT And ROBO R-HSA-428542 | 1.00e+00 | 2 | ROBO2, DCC |

### Panel E (GO MF) — top 15 (★ = matches paper)

| Rank | Term | adj p | n / overlap | Sample overlap genes |
|---:|---|---:|---:|---|
| 1 | Oxysterol Binding (GO:0008142) | 9.18e-01 | 3 | OSBPL5, INSIG2, RORA |
| 2 | Cyclic-Nucleotide Phosphodiesterase Activity (GO:0004112) | 9.18e-01 | 3 | CNP, PDE3A, PDE9A |
| 3 | 3',5'-Cyclic-Nucleotide Phosphodiesterase Activity (GO:0004114) | 9.18e-01 | 4 | PDE3A, PDE4B, PDE7B, PDE9A |
| 4 | Ligand-Gated Monoatomic Anion Channel Activity (GO:0099095) | 9.29e-01 | 4 | GABRB3, GLRA3, GLRB, GABRA5 |
| 5 | Ligand-Gated Monoatomic Ion Channel Activity Involved In Regulation Of Presynaptic Membrane Potential (GO:0099507) | 9.29e-01 | 3 | GABRA5, GRIK4, GRIN2B |
| 6 | 3',5'-cyclic-AMP Phosphodiesterase Activity (GO:0004115) | 9.29e-01 | 3 | PDE3A, PDE4B, PDE7B |
| 7 | Inhibitory Extracellular Ligand-Gated Monoatomic Ion Channel Activity (GO:0005237) | 9.29e-01 | 3 | GLRA3, GLRB, GABRA5 |
| 8 | Protein Kinase Activator Activity (GO:0030295) | 9.29e-01 | 9 | NEK9, GPRC5B, ALS2, PRKAG2, SPRY2, TAB2, CALM2, EGFR (+1) |
| 9 | Chloride Channel Activity (GO:0005254) | 9.29e-01 | 7 | GABRB3, GLRA3, TTYH2, GLRB, GABRA5, ANO2, ANO3 |
| 10 | acyl-CoA Oxidase Activity (GO:0003997) | 9.29e-01 | 2 | ACOXL, ACOX3 |
| 11 | Aspartic-Type Endopeptidase Inhibitor Activity (GO:0019828) | 9.29e-01 | 2 | CRB2, ROCK1 |
| 12 | Pyrimidine Nucleotide-Sugar Transmembrane Transporter Activity (GO:0015165) | 9.29e-01 | 2 | SLC35B4, SLC35D2 |
| 13 | Tau Protein Binding (GO:0048156) | 9.29e-01 | 5 | ROCK1, MAP2, DYRK1A, LGMN, PICALM |
| 14 | Transmitter-Gated Monoatomic Ion Channel Activity (GO:0022824) | 9.29e-01 | 5 | GABRB3, GLRA3, GABRA5, GRIK4, GRIN2B |
| 15 | Delayed Rectifier Potassium Channel Activity (GO:0005251) | 9.29e-01 | 4 | KCNG1, KCNB2, KCNQ1, KCNA4 |

## Summary table

| Input gene set | Paper-matching terms in Reactome top 15 | Paper-matching in GO MF top 15 |
|---|---:|---:|
| `epykit_top500_DMRs` | 0 | 0 |
| `epykit_top813_DMRs` | 0 | 0 |
| `epykit_top1000_DMCs_genes` | 0 | 0 |