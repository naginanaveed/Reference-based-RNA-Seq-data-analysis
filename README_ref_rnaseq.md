#  Reference-Based RNA-Seq Data Analysis (Galaxy)

> A complete step-by-step guide to detecting differentially expressed genes from raw RNA-seq reads — using the Galaxy web platform.

🔗 **Tutorial Link:** https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/ref-based/tutorial.html

---

## Table of Contents

- [Introduction](#Introduction)
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
- [DESeq2 Output](#deseq2-output)
- [References](#references)

---

## Introduction

RNA sequencing (RNA-seq) measures the expression level of every gene in a sample at the same time — a technology that captures a snapshot of the entire transcriptome (the set of all RNA molecules in a cell). One of the most common goals is **Differential Expression (DE) analysis**: comparing gene expression between two conditions (e.g., treated vs. untreated) to find which genes are turned up or down.

- In **reference-based RNA-seq**, reads are mapped to a known reference genome rather than assembled from scratch. This is appropriate when a good reference genome exists for your organism (human, mouse, *Drosophila*, etc.).

- The full pipeline goes from **raw FASTQ reads → quality control → adapter trimming → genome alignment → read counting → statistical DE analysis → visualization → pathway enrichment**.

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

>  The tutorial includes both paired-end AND single-end samples in the same experiment. This is handled separately in some steps (particularly for STAR alignment parameters) and as a combined factor in DESeq2.

---

## Dataset Overview

### Study Design
 
**Original Study:** Brooks et al. 2011  
**Organism:** *Drosophila melanogaster* (Fruit fly)  
**Condition:** Pasilla gene depletion via RNAi  
 
### Sample Information
 
| Condition | Samples | Type | Replicates |
|-----------|---------|------|-----------|
| **Untreated (Wild-type)** | GSM461176, GSM461177, GSM461178, GSM461182 | 2 Single-end, 2 Paired-end | 4 |
| **Treated (Pasilla-depleted)** | GSM461179, GSM461180, GSM461181 | 1 Single-end, 2 Paired-end | 3 |
| **Total** | 7 samples | Mixed sequencing types | 7 biological replicates |
 
### Data Formats
 
- **Input Format:** FASTQ (gzipped)
- **Output Format:** BAM (binary alignment), counts matrix (tabular), visualizations (HTML)
- **Reference Genome:** *Drosophila melanogaster* dm6 (Release 6)
- **Annotation File:** Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz
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

 
| Parameter | Value | Reason |
|-----------|-------|--------|
| **Format** | FASTQ (gzip) | Standard sequencing output format |
| **Import Method** | URL paste or upload | Galaxy auto-detects format |
| **Organization** | Create collections | Batch processing for multiple samples |
 
**Why This Step?**
- Ensures all raw sequencing data is available in Galaxy
- Establishes baseline dataset for tracking
- Creates organized history for reproducibility
**Expected Output:**
- FASTQ datasets (one per sample) in Galaxy history
- File size: typically 100 MB - 1 GB per sample
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


| Parameter | Value | Reason |
|-----------|-------|--------|
| **Input** | Falco/FastQC reports (collection) | Combines all sample QC data |
| **Report Type** | Interactive HTML | Easy visualization and filtering |
 
**Why This Tool?**
- **Aggregate Data:** Combines quality reports from all samples into one interactive report
- **Compare Samples:** Quickly identify outliers or problematic samples
- **Export Statistics:** Generate publication-quality plots
**Key Metrics Evaluated:**
 
| Metric | Ideal Value | Acceptable Range | Action if Poor |
|--------|-------------|------------------|----------------|
| **Per-base Q-score** | >30 (Q30+) | >28 | Trim low-quality bases |
| **GC Content** | Stable across reads | ±5% of mean | Indicates contamination if skewed |
| **Duplicate Rate** | 25-50% | <70% | High = possible PCR bias |
| **Read Length** | Uniform | ±2 bp variation | Low = adapter contamination |
 
**Output:** MultiQC HTML report with interactive visualizations
 
---

### Phase 3: Trimming — Cutadapt

| Parameter | Value | Reason |
|-----------|-------|--------|
| **Tool** | Cutadapt | Industry standard for adapter trimming |
| **Adapter Sequence** | Illumina TruSeq | Matches library preparation kit |
| **Quality Threshold** | PHRED ≥20 | Removes low-quality bases |
| **Minimum Length** | 20 bp | Discards too-short fragments |
| **Input Type** | Paired-end reads | Preserve read pairing information |
 
**Why This Step?**
 
1. **Remove Adapters:** Sequencing adapters don't map to genome; they create false alignments
2. **Quality Trimming:** Low-quality bases at read ends cause mapping errors
3. **Preserve Pairing:** Maintains forward/reverse read pairs for sensitive alignment
**Input:**
- R1 (forward) FASTQ
- R2 (reverse) FASTQ
**Output:**
- Trimmed R1 FASTQ
- Trimmed R2 FASTQ
- Single-end reads (if one read discarded during trimming)
- Trimming report (HTML)
**Quality Control Check:**
 
| Check | Expected Outcome |
|-------|-----------------|
| **Adapter Removal** | <1% adapter sequence remaining |
| **Length Distribution** | Peak at original length - a few bp |
| **Base Loss** | 5-15% of bases removed (typical) |
| **Pair Retention** | >90% of read pairs retained |
 
**Example Report Output:**
```
=== Cutadapt Summary ===
Reads processed: 12,345,678
Reads with adapters: 2,341,567 (19%)
Reads trimmed for quality: 3,456,789 (28%)
Reads discarded: 234,567 (1.9%)
Final read pairs: 12,110,111 (98%)
```
 
---
 

### Phase 4: Mapping — RNA STAR

**Goal:** Align the trimmed reads to the *Drosophila* reference genome (dm6) using a splice-aware aligner, producing BAM files showing exactly where each read maps.

| Parameter | Value | Reason |
|-----------|-------|--------|
| **Tool** | RNA STAR v2.7.11b+ | Ultra-fast splice-aware aligner |
| **Input Format** | Paired-end collections | Batch map all samples simultaneously |
| **Reference Genome** | *D. melanogaster* dm6 | Matches annotation file version |
| **Annotation File** | .GTF (gzip) | Enables junction detection |
| **Junction Window** | 36 bp | Length_of_reads - 1 |
| **Output Mode** | Per-gene read counts | Directly generate feature counts |
| **Threading** | 8+ cores | Galaxy provides HPC resources |
 
**Why STAR?**
 
1. **Handles Splice Junctions:** Most RNA reads span exon-exon boundaries (splice junctions)
2. **Speed:** Aligns millions of reads in minutes (vs. hours for other tools)
3. **Accuracy:** Superior sensitivity for novel splicing events
4. **Gene Quantification:** Outputs read counts per gene directly
   
**Input:**
- Trimmed R1 FASTQ (forward reads)
- Trimmed R2 FASTQ (reverse reads)
- Reference genome index (built automatically by Galaxy)
- Gene annotation GTF file
**Output:**
- **Aligned_sortedByCoord.bam:** Binary alignment file (sortable, indexed)
- **Log.final.out:** Mapping statistics (alignment rate, splice junctions, etc.)
- **ReadsPerGene.out.tab:** Per-gene raw read counts (unstranded)
- **SJ.out.tab:** Splice junction information
**Mapping Statistics Explanation:**
 
| Statistic | Ideal | Acceptable | Problem |
|-----------|-------|-----------|---------|
| **Uniquely Mapped** | >80% | >70% | <70% = contamination or wrong reference |
| **Multimapped Reads** | <10% | <20% | Repetitive sequences or low complexity |
| **Unmapped Reads** | <15% | <25% | Low quality sequencing |
| **Spliced Reads** | 40-60% | 30-70% | Organism-dependent, indicates splice detection |
 
**Example Output:**
```
STAR Mapping Results:
Total reads: 12,110,111
Uniquely mapped: 10,455,687 (86.4%)
Multimapped: 892,340 (7.4%)
Unmapped - too short: 634,291 (5.2%)
Unmapped - other: 127,793 (1.1%)
Spliced reads: 5,234,123 (43%)
Splice junctions: 156,234 novel
```
 
---

### Phase 5: Post-Mapping QC — IGV, MarkDuplicates, RSeQC

---

#### Step 5a — Visual Inspection with IGV (Integrative Genomics Viewer)

| Field | Value |
|-------|-------|
| **Tool** | IGV (external desktop application or via Galaxy) |
| **Input** | BAM file from STAR + reference genome (dm6) |
| **Output** | Visual alignment browser |

**What you're checking:**
Load the BAM file in IGV and navigate to a known well-expressed gene. You should see reads piling up on the exons — not uniformly across the gene. The view reveals **sashimi plots** (arc plots showing reads that span splice junctions), which confirm that splicing is being captured correctly. If reads pile up on introns or the distributions look wrong, this could indicate a wrong strand setting or annotation mismatch.

 #### **Step 5b: MarkDuplicates - Identify PCR Duplicates**
 
| Parameter | Value | Reason |
|-----------|-------|--------|
| **Tool** | Picard MarkDuplicates | Industry standard for duplicate detection |
| **Input** | BAM files (from STAR) | Requires aligned reads |
| **Assumptions** | Reads at same position = duplicates | Conservative approach |
| **Marking Mode** | Mark (not remove) | Preserves all reads for later filtering |
 
**Why This Step?**
 
In RNA-Seq, duplicate reads can arise from:
 
1. **Legitimate:** Highly expressed genes naturally produce many identical reads
2. **Artifacts:** PCR amplification during library prep (bad for differential expression)
The tool marks duplicates so downstream analysis can decide whether to filter them.
 
**Input:** BAM file from STAR  
**Output:** 
- BAM file with duplicate flags
- Duplicate metrics report (tabular)
**Expected Metrics:**
 
| Sample | Duplicate % | Interpretation |
|--------|-------------|----------------|
| **5-20%** | Normal | Low complexity library, normal expression variation |
| **20-50%** | Acceptable | Common in RNA-Seq, especially if lowly-sequenced |
| **>50%** | Concern | Possible PCR over-amplification or very biased sample |
 
** Report:**
```
Sample: GSM461177_paired
Total reads: 10,455,687
Duplicate reads: 2,708,000 (25.9%)
Duplicate rate indicates: Normal (within expected range)
```
 
---
#### **Step 5c: RSeQC - Quality Control Suite**
 
**i: Gene Body Coverage**
 
| Parameter | Value | Reason |
|-----------|-------|--------|
| **Tool** | RSeQC Gene Body Coverage | Detects 3' bias in cDNA prep |
| **Input** | BAM + GTF annotation | Requires gene coordinates |
| **Bins** | 100 (5% intervals) | Smooth coverage profile |
 
**Why This Tool?**
- **Detects Bias:** Some library prep methods show preferential 5' or 3' bias
- **RNA Degradation:** Non-uniform coverage suggests degraded RNA
- **Quality Indicator:** Uniform coverage = high-quality RNA
**Expected Output:**
- Coverage plot showing 5' → 3' distribution
- Ideally: Relatively flat line (uniform coverage)
- Problem: Steep 5' or 3' drop-off indicates RNA degradation
**Interpretation:**
```
 Expected: Relatively even coverage 5' to 3' (flat line)
 Concern: Sharp drop at 3' end → possible RNA degradation
 Concern: Drop at 5' end → unusual, may indicate technical issue
```
 
---
 
**ii: Infer Experiment - Strand Detection**
 
| Parameter | Value | Reason |
|-----------|-------|--------|
| **Tool** | RSeQC Infer Experiment | Determines library strandedness |
| **Input** | BAM file | Read strand information |
| **Sample Size** | 200,000 reads | Sufficient for accurate inference |
 
**Why This Step?**
 
Different library prep methods preserve strand information differently:
 
| Strandedness Type | Library Kit | Read Mapping Pattern |
|------------------|------------|-------------------|
| **Unstranded (IU)** | Standard mRNA-Seq | Reads from both strands equally |
| **Forward-stranded (SF)** | dUTP or TrueSeq | Read 1 = reverse strand, Read 2 = forward |
| **Reverse-stranded (SR)** | Older protocols | Read 1 = forward strand, Read 2 = reverse |
 
**Input:** BAM file  
**Output:** Strand fraction report
 
**Output:**
```
Fraction of reads explained by "1++,1--,2+-,2-+": 4.9%
Fraction of reads explained by "1+-,1-+,2++,2--": 95.1%
Fraction of reads explained by other combinations: 0.0%
 
Conclusion: This is a reverse-stranded library (SR type)
```
 
**Interpretation Decision:**
 
If >90% of reads fit one pattern → **Stranded library** (use in featureCounts/DESeq2)  
If ~50/50 split → **Unstranded library** (ignore strand in counting)
 
---
 
**iii: Read Distribution**
 
| Parameter | Value | Reason |
|-----------|-------|--------|
| **Tool** | RSeQC Read Distribution | Categorizes read mapping locations |
| **Input** | BAM + GTF file | Requires gene boundaries |
 
**Why This Step?**
 
Expected read distribution in clean RNA-Seq data:
 
| Feature | Expected % | Indicates |
|---------|-----------|----------|
| **Exons** | 70-90% | Correct target (mRNA) |
| **Introns** | 5-15% | Immature transcripts (normal small amount) |
| **Intergenic** | 5-15% | Contamination if high (genomic DNA) |
| **5' UTR** | 2-5% | Promoter artifacts |
| **3' UTR** | 5-10% | 3' biased amplification |
 
**Example Output:**
```
Exonic reads: 85.2% ✓
Intronic reads: 8.7% ✓
Intergenic reads: 6.1% ✓
```
 
#### **Step 5d: Aggregate with MultiQC**
 
Combines all post-mapping QC results into one interactive report for easy visualization and comparison across samples.
 
---

<img width="800" height="500" alt="mapping result igv" src="https://github.com/user-attachments/assets/89d1740c-2aec-41f9-ae8e-1d60e8d009ed" />

### Phase 6: Read Counting — featureCounts

**Goal:** Produce a count matrix — a table showing how many reads mapped to each gene in each sample. This is the direct numerical input to differential expression analysis.

| Field | Value |
|-------|-------|
| **Tool** | `featureCounts` (Subread package) |
| **Input** | BAM files (collection from STAR) · GTF annotation file |
| **Key Parameters** | Feature type: `exon` · Attribute: `gene_id` · Strand specificity: matching your library prep · Paired-end/Single-end appropriate for each sample group |
| **Output** | Count table per sample (genes as rows, read count as value) |

**Important nuances:**

- **Spliced reads at junctions:** A read spanning two exons is counted once for the gene, not twice. However, the count reflects fragments (in paired-end data), not individual reads.
- **Strand specificity:** If the library was prepared with a strand-specific protocol (most modern kits are), you must tell featureCounts which strand the reads should come from. Getting this wrong can halve your counts or assign reads to the wrong genes. Use RSeQC's `infer_experiment` to determine strandedness from your data.
- **GTF version must match:** The GTF must be the same version as the reference genome used for alignment. Different builds have different chromosome coordinates — mixing them gives wrong or zero counts.

 
**Output Matrix Format:**
 
```
Geneid          Length  GSM461176  GSM461177  GSM461179  GSM461180
FBgn0000003     1000    234        567        123        456
FBgn0000008     2500    890        1234       567        789
FBgn0000014     3200    1200       1567       234        567
...
```

### Phase 7: Differential Expression — DESeq2


| Parameter | Value | Reason |
|-----------|-------|--------|
| **Tool** | DESeq2 | Robust, widely-used R package; accounts for sequencing depth and variance |
| **Count Input** | featureCounts output | Raw counts (not normalized) |
| **Experimental Design** | 2-factor: Sample + Condition | Treated vs. Untreated |
| **Statistical Model** | Negative Binomial | Appropriate for count data |
| **Primary Factor** | Treatment (Treated/Untreated) | Effect of interest |
| **Normalization** | Media of Ratios (default) | Corrects for sequencing depth bias |
| **Dispersion Estimation** | Automatic | Galaxy estimates from data |
| **False Discovery Rate** | <0.05 (adjusted p-value) | Stringent multiple testing correction |
| **Log2 Fold Change** | ±2.0 (optional threshold) | Biologically meaningful cutoff |

**Why DESeq2?**

1. **Accounts for Variation:** Models the variance inherent in RNA-Seq counts
2. **Handles Replicates Well:** Uses all replicates to estimate gene-wise dispersion
3. **Robust:** Widely validated and cited (>50,000 citations)
4. **Multiple Corrections:** Properly adjusts p-values for multiple comparisons
5. **Normalization:** Factors out technical bias from sequencing depth differences

#### How DESeq2 Works

**Step 1 — Normalization (size factors):**
- DESeq2 calculates a **size factor** for each sample that corrects for sequencing depth. It does this by computing the geometric mean of each gene across all samples, then dividing each sample's counts by this reference. The median of these ratios becomes the size factor. After dividing counts by their sample's size factor, all samples are on a comparable scale.

**Step 2 — Dispersion estimation:**
- For each gene, DESeq2 estimates how variable its counts are across replicates. Genes with very low counts tend to be noisy, and DESeq2 "shrinks" their dispersion estimates toward a global trend to borrow statistical strength across all genes — this is especially important with few replicates.

**Step 3 — Design matrix (accounting for multiple factors):**
In this experiment there are two factors:
- **Treatment:** treated (PS depleted) vs. untreated — *the factor of interest*
- **Sequencing type:** paired-end vs. single-end — *a confounding variable to correct for*

The DESeq2 design `~ sequencing_type + treatment` tells the model to account for any systematic differences between PE and SE samples before testing for treatment effects. Without this correction, differences between PE and SE libraries would appear as false differential expression.

**Step 4 — Statistical testing (Wald test):**
For each gene, DESeq2 fits a generalized linear model and performs a Wald test for the treatment factor. This produces a p-value for the null hypothesis "this gene is not differentially expressed". P-values are adjusted for multiple testing using the Benjamini-Hochberg method (FDR control) to produce **adjusted p-values (padj)**.


**Input:**
- Count matrix (from featureCounts)
- Experimental design: Sample IDs + Condition (Treated/Untreated)
- Metadata: Which samples are biological replicates

**Output:**

1. **Statistical Results Table:**
   - Gene ID
   - Base Mean (average expression)
   - Log2 Fold Change (effect size)
   - Standard Error
   - Wald Test Statistic
   - P-value (raw)
   - Adjusted p-value (Benjamini-Hochberg)

2. **Visualizations:**
   - MA-plot (M = log2FC, A = mean log expression)
   - Volcano plot (log2FC vs. -log10 p-value)
   - Heatmap of top DE genes
   - PCA plot (sample relationships)

**Interpretation of Results:**

| Gene | log2FC | padj | Interpretation |
|------|--------|------|----------------|
| Gene A | +3.5 | 0.001 | **Significantly UP-regulated** (8.9x higher in treated) |
| Gene B | -2.1 | 0.01 | **Significantly DOWN-regulated** (4.3x lower in treated) |
| Gene C | +0.5 | 0.5 | Not significant (noisy, inconsistent) |
| Gene D | +5.0 | 0.2 | Large effect but p > 0.05 (few replicates?) |

**Statistical Significance Criteria:**

For a gene to be called **differentially expressed (DE):**
- Adjusted p-value < 0.05 (account for multiple testing)
- AND log2 Fold Change > 2 or < -0.5 (optional, depends on biological context)

**Example Results Summary:**
```
Total genes tested: 13,602
Genes up-regulated: 1,234 (padj < 0.05)
Genes down-regulated: 987 (padj < 0.05)
Total DE genes: 2,221 (16.3%)
```

<img width="500" height="500" alt="image" src="https://github.com/user-attachments/assets/4ce2b861-d523-4ae1-b1dd-0795cf8cb7c9" />



---

<img width="500" height="500" alt="image" src="https://github.com/user-attachments/assets/ad9dc2d7-ded5-48a7-9af1-a7d4f33fd900" />


---
### Phase 8: Visualization — Heatmap2 & Volcano Plot

**Goal:** Produce publication-quality visualizations of the differentially expressed genes to understand patterns and communicate results.

#### Step 8.1 — Volcano Plot

| Field | Value |
|-------|-------|
| **Tool** | `Volcano Plot` |
| **Input** | DESeq2 results table |
| **Output** | Scatter plot of −log10(padj) on y-axis vs. logFoldChange on x-axis |

**How to read a Volcano plot:**
- **X-axis:** Effect size — how strongly the gene changes. Positive = upregulated in treated; negative = downregulated.
- **Y-axis:** Statistical significance — higher = more significant. Points near the bottom are not significant.
- **Top-left:** Strongly downregulated, highly significant genes
- **Top-right:** Strongly upregulated, highly significant genes
- **Bottom centre:** Unchanged genes (most genes)


<img width="500" height="500" alt="image" src="https://github.com/user-attachments/assets/4b0c3d21-7155-4c96-a0e0-6efceff6ddec" />

---

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
<img width="500" height="500" alt="image" src="https://github.com/user-attachments/assets/71c3e622-d2fe-4ec0-8e73-2154008a0f53" />


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


**Interpreting results:**
Each row is a GO term. A small adjusted p-value means the genes annotated with that term appear more often in your DE list than expected by chance. For a Pasilla knockdown (a splicing regulator), you would expect to see enrichment of terms like **"mRNA splicing"**, **"RNA processing"**, and **"alternative splicing"** — which confirms the biological validity of the analysis.

---

<img width="500" height="500" alt="image" src="https://github.com/user-attachments/assets/65f2bd9d-963e-4921-9cdb-df8b0a3ce02a" />

---

## Tool Summary 

| Tool | Function | Input | Output | use |
|------|----------|-------|--------|----------|
| **Falco/FastQC** | Quality assessment | FASTQ | HTML report | Detect sequencing problems early |
| **MultiQC** | Aggregate reports | Multiple reports | Interactive HTML | Compare samples, spot outliers |
| **Cutadapt** | Trim adapters & low quality | FASTQ | Trimmed FASTQ | Remove non-genomic sequences |
| **STAR** | Map to genome (splice-aware) | FASTQ + reference | BAM, counts | Ultra-fast, handles introns |
| **Picard MarkDuplicates** | Flag PCR duplicates | BAM | Flagged BAM | Identify amplification bias |
| **RSeQC** | QC suite (coverage, strand, dist) | BAM + GTF | Multiple reports | Comprehensive quality evaluation |
| **featureCounts** | Count reads per gene | BAM + GTF | Count matrix | Fast, accurate quantification |
| **DESeq2** | Differential expression testing | Count matrix + metadata | DE table, plots | Statistical significance testing |
| **g:Profiler** | Gene ontology enrichment | Gene list | Enrichment results | Interpret biological function |

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

