# 🧬 Reference-Based RNA-Seq Data Analysis (Galaxy)

> A complete step-by-step guide to detecting differentially expressed genes from raw RNA-seq reads — using the Galaxy web platform, no programming required.

🔗 **Tutorial Link:** https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/ref-based/tutorial.html

**Authors:** Bérénice Batut, Mallory Freeberg, Mo Heydarian, Anika Erxleben, Pavankumar Videm, Clemens Blank, Maria Doyle, Nicola Soranzo, Peter van Heusden, Lucille Delisle
**Last Updated:** December 2025 | **License:** MIT | **Platform:** Galaxy (usegalaxy.org / Galaxy Europe)

---

## 📚 Table of Contents

- [What Is Reference-Based RNA-Seq?](#what-is-reference-based-rna-seq)
- [Biological Background — The Pasilla Experiment](#biological-background--the-pasilla-experiment)
- [Dataset Overview](#dataset-overview)
- [Full Pipeline at a Glance](#full-pipeline-at-a-glance)
- [Step-by-Step Workflow](#step-by-step-workflow)
  - [Phase 1: Setup — Data Upload & Collection Organization](#phase-1-setup--data-upload--collection-organization)
  - [Phase 2: Quality Control — Falco + MultiQC](#phase-2-quality-control--falco--multiqc)
  - [Phase 3: Trimming — Cutadapt](#phase-3-trimming--cutadapt)
  - [Phase 4: Mapping — RNA STAR](#phase-4-mapping--rna-star)
  - [Phase 5: Post-Mapping QC — IGV, MarkDuplicates, RSeQC](#phase-5-post-mapping-qc--igv-markduplicates-rseqc)
  - [Phase 6: Read Counting — featureCounts](#phase-6-read-counting--featurecounts)
  - [Phase 7: Differential Expression — DESeq2](#phase-7-differential-expression--deseq2)
  - [Phase 8: Visualization — Heatmap2 & Volcano Plot](#phase-8-visualization--heatmap2--volcano-plot)
  - [Phase 9: Functional Enrichment — goseq](#phase-9-functional-enrichment--goseq)
- [Tool Summary Table](#tool-summary-table)
- [DESeq2 Output Explained](#deseq2-output-explained)
- [References](#references)

---

## What Is Reference-Based RNA-Seq?

RNA sequencing (RNA-seq) measures the expression level of every gene in a sample at the same time — a technology that captures a snapshot of the entire transcriptome (the set of all RNA molecules in a cell). One of the most common goals is **Differential Expression (DE) analysis**: comparing gene expression between two conditions (e.g., treated vs. untreated) to find which genes are turned up or down.

In **reference-based RNA-seq**, reads are mapped to a known reference genome rather than assembled from scratch. This is appropriate when a good reference genome exists for your organism (human, mouse, *Drosophila*, etc.).

**Why is this more complex than DNA sequencing?**
Because mRNA in eukaryotes has been spliced — introns have been removed. A read that spans a splice junction cannot be mapped continuously to genomic DNA. You need a **splice-aware aligner** like STAR (used in this tutorial) that knows about splice sites and can map reads that cross exon-exon boundaries.

The full pipeline goes from **raw FASTQ reads → quality control → adapter trimming → genome alignment → read counting → statistical DE analysis → visualization → pathway enrichment**.

---

## Biological Background — The Pasilla Experiment

This tutorial is based on the study by **Brooks et al. (2011)**, which investigated the role of the **Pasilla (PS) gene** in *Drosophila melanogaster*.

- **What is Pasilla?** The Drosophila homologue of mammalian splicing regulators Nova-1 and Nova-2 — proteins that regulate alternative splicing of pre-mRNA
- **What did the study do?** Knocked down the PS gene using **RNA interference (RNAi)** in *Drosophila* cells, then sequenced total RNA from both treated (PS-depleted) and untreated cells
- **Goal:** Identify which genes change in expression when PS is removed — revealing which genes PS normally regulates

**This tutorial uses 7 samples from the original study:**

| Sample ID | Condition | Sequencing Type |
|-----------|-----------|-----------------|
| GSM461176 | Untreated | Paired-end |
| GSM461177 | Untreated | Paired-end |
| GSM461178 | Untreated | Paired-end |
| GSM461182 | Untreated | Single-end |
| GSM461179 | Treated (PS depleted) | Paired-end |
| GSM461180 | Treated (PS depleted) | Paired-end |
| GSM461181 | Treated (PS depleted) | Single-end |

> ⚠️ The tutorial includes both paired-end AND single-end samples in the same experiment. This is handled separately in some steps (particularly for STAR alignment parameters) and as a combined factor in DESeq2.

---

## Dataset Overview

| Property | Value |
|----------|-------|
| Organism | *Drosophila melanogaster* |
| Reference genome | dm6 (BDGP6.32) |
| Annotation file | `Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz` |
| Read type | Paired-end (4 samples) + Single-end (3 samples) |
| Total samples | 7 (4 untreated, 3 treated) |
| Biological replicates | ≥ 3 per condition |
| Data source | Zenodo (doi: 10.5281/zenodo.1185122) |

---

## Full Pipeline at a Glance

```
 Raw FASTQ Files (7 samples: paired-end + single-end)
         │
         ▼
 ┌───────────────────┐
 │  Phase 1: Setup   │  ← Upload files, organize into Galaxy collections
 └───────────────────┘
         │
         ▼
 ┌──────────────────────────────┐
 │  Phase 2: Quality Control    │  ← Falco / FastQC + MultiQC
 └──────────────────────────────┘
         │
         ▼
 ┌─────────────────────────┐
 │  Phase 3: Trimming      │  ← Cutadapt (adapter + quality trimming)
 └─────────────────────────┘
         │
         ▼
 ┌─────────────────────────────────────┐
 │  Phase 4: Splice-aware Mapping      │  ← RNA STAR → BAM files
 └─────────────────────────────────────┘
         │
         ▼
 ┌────────────────────────────────────────────────────┐
 │  Phase 5: Post-Mapping QC                          │
 │  IGV (visual check) + MarkDuplicates + RSeQC       │
 └────────────────────────────────────────────────────┘
         │
         ▼
 ┌──────────────────────────────────┐
 │  Phase 6: Read Counting          │  ← featureCounts → count matrix
 └──────────────────────────────────┘
         │
         ▼
 ┌─────────────────────────────────────────┐
 │  Phase 7: Differential Expression       │  ← DESeq2
 └─────────────────────────────────────────┘
         │
         ▼
 ┌─────────────────────────────────────────────────┐
 │  Phase 8: Visualization                         │
 │  Heatmap2 (sample clustering) + Volcano Plot    │
 └─────────────────────────────────────────────────┘
         │
         ▼
 ┌───────────────────────────────────┐
 │  Phase 9: Functional Enrichment   │  ← goseq → GO term analysis
 └───────────────────────────────────┘
```

---

## Step-by-Step Workflow

---

### Phase 1: Setup — Data Upload & Collection Organization

**Goal:** Get the raw FASTQ files and annotation files into Galaxy and organize them so multi-sample tools can process them efficiently.

#### Step 1.1 — Create a New Galaxy History

Log in at usegalaxy.org (or Galaxy Europe). Click the **+** icon in the History panel and name your history (e.g., "Reference-based RNA-Seq — Pasilla").

#### Step 1.2 — Import FASTQ Files from Zenodo

```
Data source: https://doi.org/10.5281/zenodo.1185122
```

Use **Upload Data → Paste/Fetch Data** and paste the Zenodo URLs. Make sure each file is set to format `fastqsanger` (not `fastq` — Galaxy treats these differently for quality score encoding).

#### Step 1.3 — Import the GTF Annotation File

Import `Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz` from Zenodo or Galaxy's Shared Data Library.

**Why do you need a GTF file?**
The GTF (Gene Transfer Format) file is a text annotation file that tells the aligner and counting tools where every gene, exon, and transcript is located in the reference genome. It is used twice: first by STAR during alignment (to identify splice sites), and again by featureCounts during read counting (to assign reads to genes). It is critical that the GTF file matches the same genome build (dm6) used for alignment — coordinates in different builds are different and mixing them gives wrong results.

#### Step 1.4 — Organize Files into Paired Collections

For the paired-end samples, Galaxy needs the R1 (forward) and R2 (reverse) reads linked together in a **paired collection**. Use **Build List of Dataset Pairs** to match each `_R1.fastq` with its corresponding `_R2.fastq` by sample name.

**Why use collections?**
Collections allow Galaxy to run a tool on all samples at once in one operation, rather than manually submitting each sample separately. This is essential for a 7-sample experiment and becomes critical for larger studies.

---

### Phase 2: Quality Control — Falco + MultiQC

**Goal:** Assess the quality of the raw reads to detect problems before any processing. Problems at this stage — if missed — corrupt all downstream results.

#### Step 2.1 — Run Falco on All Samples

| Field | Value |
|-------|-------|
| **Tool** | `Falco` (efficiency-optimized rewrite of FastQC) |
| **Input** | Raw FASTQ files (collection) |
| **Output** | One HTML quality report per sample |
| **Galaxy version** | Current |

**What Falco checks and why each matters:**

| QC Module | What It Checks | What to Look For |
|-----------|---------------|-----------------|
| Per-base sequence quality | Phred quality score at each position in the read | Should be >20 (ideally >28) across all positions; drop at 3' end is normal |
| Per-sequence quality scores | Distribution of mean quality per read | Most reads should have mean quality >30 |
| Per-base sequence content | A/T/C/G proportion at each position | Should be roughly equal; bias at first few positions is normal for RNA-seq |
| GC content | Overall GC% distribution | Should follow a bell curve; sharp peaks can indicate contamination |
| Sequence duplication levels | % reads that are exact duplicates | High duplication in RNA-seq is normal (highly expressed genes); flag only if extreme |
| Adapter content | % reads containing adapter sequences | Any adapter > 5% means trimming is needed |
| Overrepresented sequences | Sequences appearing far more than expected | Could be adapters, rRNA, or contamination |

#### Step 2.2 — Aggregate with MultiQC

| Field | Value |
|-------|-------|
| **Tool** | `MultiQC` |
| **Input** | Collection of all Falco reports |
| **Output** | Single unified HTML report for all samples side by side |

**Why MultiQC is essential:**
Reviewing 7 individual Falco reports would be slow and it's easy to miss patterns. MultiQC combines them into one interactive report where you can immediately see if one sample has a quality problem compared to others. It also works with many other tools (STAR, featureCounts, Cutadapt) to aggregate their outputs in later steps.

---

### Phase 3: Trimming — Cutadapt

**Goal:** Remove adapter sequences and low-quality bases from the ends of reads before alignment. Adapters are synthetic sequences from the library preparation protocol — if not removed, they interfere with mapping.

| Field | Value |
|-------|-------|
| **Tool** | `Cutadapt` |
| **Input** | Raw FASTQ collection (or paired collection for PE data) |
| **Key Parameters** | Adapter sequence (Illumina TruSeq: `AGATCGGAAGAGC`); quality cutoff 20; minimum length 20 |
| **Output** | Trimmed FASTQ collection |

**What Cutadapt does step by step:**
1. **Adapter trimming:** Scans for the adapter sequence at the 3' end of each read and removes it along with everything downstream. Without this, the aligner would try to map the adapter sequence to the genome, fail, and discard the read.
2. **Quality trimming:** Removes low-quality bases from the 3' end using a sliding window approach. Poor-quality bases at the end of reads reduce alignment accuracy.
3. **Length filtering:** Discards reads that become too short after trimming (< 20 bp by default). Very short reads map to multiple places in the genome (multi-mapping) and are unreliable.

**After trimming:**
Re-run Falco + MultiQC on the trimmed reads to confirm that adapter content is gone and quality has improved.

| Before Cutadapt | After Cutadapt |
|-----------------|---------------|
| Adapter sequences visible in Falco | No adapter content |
| Low-quality bases at read 3' end | Quality remains high to end |
| All reads same length | Reads have variable lengths (different amounts trimmed) |

> 💡 **Alternative tools:** Trim Galore! or Trimmomatic can be used instead of Cutadapt for the same purpose.

---

### Phase 4: Mapping — RNA STAR

**Goal:** Align the trimmed reads to the *Drosophila* reference genome (dm6) using a splice-aware aligner, producing BAM files showing exactly where each read maps.

**Why a splice-aware aligner?**
Most reads in RNA-seq come from processed mRNA — meaning introns have been removed. A read at an exon-exon junction spans a region that does not exist as a continuous sequence in the genome. Standard DNA aligners cannot handle this. STAR uses the GTF annotation to know where splice sites are and maps junction-spanning reads correctly.

| Field | Value |
|-------|-------|
| **Tool** | `RNA STAR` (Galaxy version 2.7.11b+galaxy0) |
| **Input** | Trimmed FASTQ collection · Reference genome: dm6 · GTF annotation file |
| **Key Parameters** | Junction overhang = read length − 1 (36 for 37bp reads) · Output: Per gene read counts (GeneCounts) · Single-end or Paired-end as appropriate |
| **Output** | BAM file (aligned reads) · Alignment log · Splice junction file · Gene count table (from STAR directly) |

**Key STAR parameters explained:**

| Parameter | Why It Matters |
|-----------|----------------|
| `--sjdbGTFfile` | GTF file for splice junction database; STAR uses this to know valid splice sites |
| `--sjdbOverhang` | Should be read_length − 1; defines how far a read can extend past a junction |
| `--outSAMtype BAM SortedByCoordinate` | Produces a coordinate-sorted BAM, required by downstream tools |
| Per gene read counts (GeneCounts) | Optionally counts reads per gene during mapping; equivalent to HTSeq-count output |

**What to check in the STAR alignment log:**
- **% Uniquely mapped reads:** Should be >80%. Lower rates suggest poor quality data, wrong reference genome, or significant contamination.
- **% Multi-mapped reads:** Reads mapping to >1 location. A small % is expected (repetitive elements); very high % is a problem.
- **% Reads unmapped:** Reads that could not be placed. Investigate if >10%.

---

### Phase 5: Post-Mapping QC — IGV, MarkDuplicates, RSeQC

**Goal:** Verify the quality of the alignment more deeply — checking that reads are going to the right places, that duplicate rates are acceptable, and that the read distribution across genomic features makes sense for RNA-seq data.

#### Step 5.1 — Visual Inspection with IGV (Integrative Genomics Viewer)

| Field | Value |
|-------|-------|
| **Tool** | IGV (external desktop application or via Galaxy) |
| **Input** | BAM file from STAR + reference genome (dm6) |
| **Output** | Visual alignment browser |

**What you're checking:**
Load the BAM file in IGV and navigate to a known well-expressed gene. You should see reads piling up on the exons — not uniformly across the gene. The view reveals **sashimi plots** (arc plots showing reads that span splice junctions), which confirm that splicing is being captured correctly. If reads pile up on introns or the distributions look wrong, this could indicate a wrong strand setting or annotation mismatch.

#### Step 5.2 — Duplicate Read Marking with MarkDuplicates

| Field | Value |
|-------|-------|
| **Tool** | `MarkDuplicates` (Picard suite) |
| **Input** | Sorted BAM files (collection from STAR) |
| **Output** | BAM with duplicate reads flagged · Duplicate metrics log |

**What are duplicate reads?**
Duplicate reads are multiple copies of the same original DNA fragment that were amplified during PCR library preparation and then sequenced multiple times. They are not independent measurements of gene expression.

**Why keep them in RNA-seq (usually)?**
In RNA-seq, duplicates mostly come from highly expressed genes — a gene expressed at very high levels will naturally produce many identical fragments even without PCR over-amplification. Removing them would unfairly penalize high-expression genes. For this reason, duplicates are typically **marked but not removed** in RNA-seq differential expression analysis. They are only cause for concern if the overall duplication rate is very high (>60–70%), which could indicate over-amplification of a low-complexity library.

#### Step 5.3 — Read Distribution with RSeQC

| Field | Value |
|-------|-------|
| **Tool** | `RSeQC Read Distribution` |
| **Input** | BAM files (collection) · BED12 reference gene model (converted from GTF) |
| **Output** | Text report showing how reads are distributed across genomic features |

**Expected distribution for successful RNA-seq:**

| Feature | Expected % | Meaning if Different |
|---------|-----------|---------------------|
| CDS Exons | ~60–80% | Good: reads are in coding exons |
| 5' UTR / 3' UTR | ~5–15% | Normal; mRNA includes UTRs |
| Introns | ~2–10% | Low is good; high may suggest genomic DNA contamination or unspliced pre-mRNA |
| Intergenic regions | ~2–10% | Very high would suggest non-specific sequencing or wrong annotation |

> 🔑 **Expected result:** >80% of reads mapped to exons confirms the data is genuine RNA-seq and the alignment was successful.

Run **MultiQC** on the Read Distribution outputs to compare all samples side by side.

---

### Phase 6: Read Counting — featureCounts

**Goal:** Produce a count matrix — a table showing how many reads mapped to each gene in each sample. This is the direct numerical input to differential expression analysis.

| Field | Value |
|-------|-------|
| **Tool** | `featureCounts` (Subread package) |
| **Input** | BAM files (collection from STAR) · GTF annotation file |
| **Key Parameters** | Feature type: `exon` · Attribute: `gene_id` · Strand specificity: matching your library prep · Paired-end/Single-end appropriate for each sample group |
| **Output** | Count table per sample (genes as rows, read count as value) |

**What featureCounts does:**
For each gene in the GTF file, featureCounts counts the number of reads from the BAM file that overlap with the gene's annotated exons. It reports one count per gene — the total reads across all exons.

**Important nuances:**

- **Spliced reads at junctions:** A read spanning two exons is counted once for the gene, not twice. However, the count reflects fragments (in paired-end data), not individual reads.
- **Strand specificity:** If the library was prepared with a strand-specific protocol (most modern kits are), you must tell featureCounts which strand the reads should come from. Getting this wrong can halve your counts or assign reads to the wrong genes. Use RSeQC's `infer_experiment` to determine strandedness from your data.
- **GTF version must match:** The GTF must be the same version as the reference genome used for alignment. Different builds have different chromosome coordinates — mixing them gives wrong or zero counts.

**Output count table (example):**

| GeneID | GSM461176 | GSM461177 | GSM461178 | GSM461179 | GSM461180 |
|--------|-----------|-----------|-----------|-----------|-----------|
| FBgn0000003 | 0 | 0 | 0 | 0 | 0 |
| FBgn0000008 | 92 | 79 | 85 | 12 | 8 |
| FBgn0000014 | 4521 | 4831 | 4320 | 3892 | 4102 |

Run **MultiQC** on featureCounts outputs to check the assignment rates across samples.

> 💡 **Alternative tools:** HTSeq-count produces equivalent results; STAR can also count directly during alignment (identical to HTSeq-count). featureCounts is preferred for speed.

---

### Phase 7: Differential Expression — DESeq2

**Goal:** Using the count tables from all 7 samples, statistically identify genes that are significantly more or less expressed in the treated (PS-depleted) condition compared to untreated.

| Field | Value |
|-------|-------|
| **Tool** | `DESeq2` |
| **Input** | 7 count tables (one per sample) · Factor definitions (treatment status, sequencing type) |
| **Output** | Normalized count table · Summary results table (log2 fold change, p-values) · QC plots (PCA, sample distance heatmap, MA plot) |

#### Why DESeq2 (and not just comparing averages)?

Raw read counts cannot be directly compared between samples because:
1. Samples are sequenced to different depths — a sample with twice the total reads will have roughly double the counts for every gene, regardless of true expression
2. Gene expression data is highly variable and overdispersed — the variance of counts is much greater than the mean, violating assumptions of simple statistical tests
3. With only 3–4 replicates per condition, statistical power is limited and needs to be handled carefully

DESeq2 addresses all of these with a principled statistical framework based on the **negative binomial distribution**.

#### How DESeq2 Works

**Step 1 — Normalization (size factors):**
DESeq2 calculates a **size factor** for each sample that corrects for sequencing depth. It does this by computing the geometric mean of each gene across all samples, then dividing each sample's counts by this reference. The median of these ratios becomes the size factor. After dividing counts by their sample's size factor, all samples are on a comparable scale.

**Step 2 — Dispersion estimation:**
For each gene, DESeq2 estimates how variable its counts are across replicates. Genes with very low counts tend to be noisy, and DESeq2 "shrinks" their dispersion estimates toward a global trend to borrow statistical strength across all genes — this is especially important with few replicates.

**Step 3 — Design matrix (accounting for multiple factors):**
In this experiment there are two factors:
- **Treatment:** treated (PS depleted) vs. untreated — *the factor of interest*
- **Sequencing type:** paired-end vs. single-end — *a confounding variable to correct for*

The DESeq2 design `~ sequencing_type + treatment` tells the model to account for any systematic differences between PE and SE samples before testing for treatment effects. Without this correction, differences between PE and SE libraries would appear as false differential expression.

**Step 4 — Statistical testing (Wald test):**
For each gene, DESeq2 fits a generalized linear model and performs a Wald test for the treatment factor. This produces a p-value for the null hypothesis "this gene is not differentially expressed". P-values are adjusted for multiple testing using the Benjamini-Hochberg method (FDR control) to produce **adjusted p-values (padj)**.

#### Interpreting DESeq2 Results Table

| Column | Description |
|--------|-------------|
| `Gene ID` | Drosophila gene identifier (e.g., FBgn0000014) |
| `baseMean` | Average normalized count across all samples — indicates how much the gene is expressed overall |
| `log2FoldChange` | How many times more (positive) or less (negative) the gene is expressed in treated vs. untreated, on log2 scale. log2FC = 1 means 2× higher; log2FC = −2 means 4× lower |
| `lfcSE` | Standard error of the fold change estimate — smaller is more precise |
| `stat` | Wald statistic (log2FC / lfcSE) |
| `pvalue` | Raw p-value from the Wald test |
| `padj` | Benjamini-Hochberg adjusted p-value (FDR-corrected) |

**Standard cutoffs for calling a gene differentially expressed:**
- `padj < 0.05` — statistically significant (at most 5% false discovery rate)
- `|log2FoldChange| > 1` — at least 2-fold change in expression (biologically meaningful)

#### DESeq2 QC Plots

| Plot | What It Shows | What Good Looks Like |
|------|--------------|---------------------|
| **PCA plot** | Samples projected onto 2 principal components based on normalized expression | PC1 should separate treated from untreated; PC2 separates PE from SE. No outlier samples. |
| **Sample distance heatmap** | Hierarchical clustering of sample-to-sample expression similarity | Treated samples cluster together; untreated samples cluster together |
| **MA plot** | log2FC on y-axis vs. mean expression (baseMean) on x-axis for every gene | Most genes cluster around log2FC = 0 (unchanged); significant genes shown in red |
| **Dispersion plot** | Gene-wise dispersion vs. mean expression | Shows DESeq2's shrinkage of dispersions toward the fitted trend |

> 💡 **Alternative tools:** edgeR and limma-voom are valid alternatives to DESeq2 and often give similar results. limma-voom performs well with more replicates.

---

### Phase 8: Visualization — Heatmap2 & Volcano Plot

**Goal:** Produce publication-quality visualizations of the differentially expressed genes to understand patterns and communicate results.

#### Step 8.1 — Volcano Plot

| Field | Value |
|-------|-------|
| **Tool** | `Volcano Plot` |
| **Input** | DESeq2 results table |
| **Output** | Scatter plot of −log10(padj) on y-axis vs. log2FoldChange on x-axis |

**How to read a Volcano plot:**
- **X-axis:** Effect size — how strongly the gene changes. Positive = upregulated in treated; negative = downregulated.
- **Y-axis:** Statistical significance — higher = more significant. Points near the bottom are not significant.
- **Top-left:** Strongly downregulated, highly significant genes
- **Top-right:** Strongly upregulated, highly significant genes
- **Bottom centre:** Unchanged genes (most genes)

A well-conducted experiment shows a clear symmetric spread with a meaningful number of significant genes in the top corners.

#### Step 8.2 — Heatmap2

| Field | Value |
|-------|-------|
| **Tool** | `Heatmap2` |
| **Input** | Normalized counts for significant DE genes (filtered from DESeq2 output) |
| **Output** | Clustered heatmap showing expression patterns across all samples |

**What the heatmap shows:**
Rows are genes, columns are samples. Color intensity shows normalized expression level (typically red = high, blue = low, relative to row mean). Hierarchical clustering groups genes with similar expression patterns together and groups samples by similarity. The heatmap allows you to visually confirm that:
- Treated samples form one cluster and untreated form another
- Upregulated genes are uniformly high in treated, low in untreated
- Downregulated genes show the opposite pattern

---

### Phase 9: Functional Enrichment — goseq

**Goal:** Go beyond individual genes and ask: are there biological pathways, processes, or functions that are systematically enriched among the DE genes? This converts a gene list into biological insight.

| Field | Value |
|-------|-------|
| **Tool** | `goseq` |
| **Input** | List of significantly DE gene IDs · Background gene list (all expressed genes tested) · Organism: *Drosophila melanogaster* |
| **Output** | Table of enriched Gene Ontology (GO) terms with p-values |

**What is Gene Ontology (GO)?**
The Gene Ontology is a standardized vocabulary that describes biological attributes of genes across species. Every gene can be annotated with one or more GO terms in three categories:
- **Biological Process (BP):** What biological process the gene is involved in (e.g., "alternative mRNA splicing", "cell cycle")
- **Molecular Function (MF):** What the gene product does at the molecular level (e.g., "RNA binding", "kinase activity")
- **Cellular Component (CC):** Where the gene product is located in the cell (e.g., "nucleus", "mitochondria")

**What goseq does and why it's needed:**
A naive approach would be to test whether each GO term's genes are over-represented among your DE gene list using a Fisher's test. However, in RNA-seq data, **longer genes tend to get more reads** and are therefore more likely to be detected as differentially expressed — not because they are biologically more important. goseq corrects for this **gene-length bias** before performing enrichment testing, making the results much more reliable than simple Fisher's tests.

**Interpreting results:**
Each row is a GO term. A small adjusted p-value means the genes annotated with that term appear more often in your DE list than expected by chance. For a Pasilla knockdown (a splicing regulator), you would expect to see enrichment of terms like **"mRNA splicing"**, **"RNA processing"**, and **"alternative splicing"** — which confirms the biological validity of the analysis.

---

## Tool Summary Table

| Phase | Tool | Galaxy Version | Input | Output | Purpose |
|-------|------|---------------|-------|--------|---------|
| QC | Falco | Current | Raw FASTQ | Per-sample HTML QC report | Assess raw read quality |
| QC Aggregation | MultiQC | 1.27+galaxy4 | Falco reports | Combined HTML report | Compare all samples in one view |
| Trimming | Cutadapt | Current | Raw FASTQ | Trimmed FASTQ | Remove adapters and low-quality bases |
| Alignment | RNA STAR | 2.7.11b+galaxy0 | Trimmed FASTQ + GTF + dm6 | BAM + log + junctions | Splice-aware genome alignment |
| Alignment QC | IGV | External | BAM | Visual alignment browser | Visual check of mapping quality |
| Duplicate Marking | MarkDuplicates (Picard) | Current | BAM | Flagged BAM + metrics | Identify PCR duplicates |
| Read Distribution | RSeQC | Current | BAM + BED12 | Distribution report | Confirm reads map to exons |
| Read Counting | featureCounts | Current | BAM + GTF | Count matrix | Count reads per gene |
| Differential Expression | DESeq2 | Current | Count tables | Results table + QC plots | Statistical DE analysis |
| Visualization | Volcano Plot | Current | DESeq2 results | Volcano plot | Visualize significance vs. fold change |
| Visualization | Heatmap2 | Current | Normalized counts | Clustered heatmap | Show expression patterns |
| Enrichment | goseq | Current | DE gene list | GO enrichment table | Identify enriched pathways |

---

## DESeq2 Output Explained

DESeq2 produces three main outputs:

**Output 1 — Normalized counts table:**
A matrix of all genes (rows) × all samples (columns), where each value is the raw count divided by the sample's size factor. Use this for visualization (heatmaps) and for checking individual gene expression levels.

**Output 2 — Graphical summary:**
Contains the PCA plot, sample-distance heatmap, MA plot, and dispersion plot — all described in Phase 7 above. Review these before trusting the results table.

**Output 3 — Results summary table:**
The main output: one row per gene with log2FoldChange, p-value, and padj. Apply these filters to get your final DE gene list:
- `padj < 0.05` (statistically significant after multiple testing correction)
- `|log2FoldChange| ≥ 1` (at least 2× change — biologically meaningful)

**Number of replicates matters:**
This tutorial uses 3–4 biological replicates per condition. DESeq2 and similar tools require **at least 3 biological replicates** per condition to estimate dispersion reliably. With fewer replicates, statistical power drops sharply and the number of false negatives increases significantly.

---

## References

- Brooks AN et al. (2011). *Conservation of an RNA regulatory map between Drosophila and mammals*. Genome Research 21:193–202.
- Dobin A et al. (2013). *STAR: ultrafast universal RNA-seq aligner*. Bioinformatics 29:15–21.
- Liao Y, Smyth GK, Shi W (2014). *featureCounts: an efficient general purpose program for assigning sequence reads to genomic features*. Bioinformatics 30:923–930.
- Love MI, Huber W, Anders S (2014). *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2*. Genome Biology 15:550.
- Martin M (2011). *Cutadapt removes adapter sequences from high-throughput sequencing reads*. EMBnet.journal 17:10–12.
- Ewels P et al. (2016). *MultiQC: summarize analysis results for multiple tools and samples in a single report*. Bioinformatics 32:3047–3048.
- Batut B et al. (2018). *Community-Driven Data Analysis Training for Biology*. Cell Systems. doi: 10.1016/j.cels.2018.05.012.

---

*README compiled from the Galaxy Training Network tutorial. For the most current version of this tutorial, visit the [official tutorial page](https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/ref-based/tutorial.html).*
