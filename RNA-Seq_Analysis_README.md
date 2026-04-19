# Reference-Based RNA-Seq Data Analysis in Galaxy

> A comprehensive guide to performing RNA-Seq analysis using the Galaxy platform — from raw sequencing reads to differential expression analysis and functional enrichment

![Galaxy](https://img.shields.io/badge/Platform-Galaxy-blue)
![Status](https://img.shields.io/badge/Status-Complete-green)
![License](https://img.shields.io/badge/License-CC%20BY%204.0-yellow)

---

## 📋 Table of Contents

- [Overview](#overview)
- [Dataset Description](#dataset-description)
- [Workflow Architecture](#workflow-architecture)
- [Detailed Step-by-Step Analysis](#detailed-step-by-step-analysis)
  - [1. Data Import & QC](#1-data-import--quality-control)
  - [2. Read Trimming](#2-read-trimming)
  - [3. Read Mapping](#3-read-mapping)
  - [4. Post-Mapping Quality Checks](#4-post-mapping-quality-checks)
  - [5. Feature Counting](#5-feature-counting)
  - [6. Differential Expression Analysis](#6-differential-expression-analysis)
  - [7. Visualization & Enrichment](#7-visualization--enrichment-analysis)
- [Tool Overview](#tool-overview)
- [Quality Control Checkpoints](#quality-control-checkpoints)
- [Parameters & Best Practices](#parameters--best-practices)
- [Expected Results](#expected-results)
- [Troubleshooting](#troubleshooting)
- [References](#references)

---

## 🎯 Overview

This tutorial demonstrates a complete **reference-based RNA-Seq analysis workflow** using the Galaxy platform. The analysis:

- **Analyzes** gene expression changes in *Drosophila melanogaster* cells following *Pasilla* gene depletion
- **Compares** treated (PS-depleted via RNAi) vs. untreated samples
- **Identifies** differentially expressed genes and pathways
- **Provides** reproducible, web-based bioinformatics analysis

### Key Advantages of Galaxy

| Feature | Benefit |
|---------|---------|
| **No Installation Required** | All tools pre-installed and configured |
| **High Performance Computing** | Access to HPC clusters without local hardware |
| **Reproducibility** | Complete workflow history and dataset tracking |
| **Ease of Use** | Web interface eliminates command-line barriers |
| **Collaborative** | Easy sharing and publication of workflows |

---

## 📊 Dataset Description

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

## 🔄 Workflow Architecture

```
Raw FASTQ Data
       ↓
┌─────────────────────────────────────────┐
│  STEP 1: Quality Control Assessment     │
│  • Falco/FastQC → Quality reports       │
│  • MultiQC → Aggregated QC stats        │
└─────────────────────────────────────────┘
       ↓
┌─────────────────────────────────────────┐
│  STEP 2: Adapter & Quality Trimming     │
│  • Cutadapt → Remove adapters           │
│  • Trim low-quality bases               │
└─────────────────────────────────────────┘
       ↓
┌─────────────────────────────────────────┐
│  STEP 3: Read Mapping to Genome         │
│  • STAR/HISAT2 → BAM files              │
│  • Splice-aware alignment               │
└─────────────────────────────────────────┘
       ↓
┌─────────────────────────────────────────┐
│  STEP 4: Post-Mapping Quality Checks    │
│  • MarkDuplicates → Identify PCR dups   │
│  • RSeQC → Coverage & strandness        │
│  • InferExperiment → Strand detection   │
└─────────────────────────────────────────┘
       ↓
┌─────────────────────────────────────────┐
│  STEP 5: Feature Counting               │
│  • featureCounts → Gene-level counts    │
│  • Generate count matrix                │
└─────────────────────────────────────────┘
       ↓
┌─────────────────────────────────────────┐
│  STEP 6: Differential Expression        │
│  • DESeq2 → Statistical testing         │
│  • Identify DE genes                    │
└─────────────────────────────────────────┘
       ↓
Downstream Visualization & Enrichment
```

---

## 📖 Detailed Step-by-Step Analysis

### 1. Data Import & Quality Control

#### **Step 1a: Import FASTQ Files**

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

#### **Step 1b: Quality Assessment with Falco/FastQC**

| Parameter | Value | Reason |
|-----------|-------|--------|
| **Tool** | Falco (FastQC replacement) | Faster, lighter, identical output |
| **Read Length** | 36-100 bp | Standard Illumina read lengths |
| **Adapter Content** | Check enabled | Detect sequencing adapters |
| **Duplicate Content** | Check enabled | Expected in RNA-Seq (highly expressed genes) |

**Why This Tool?**
- **Assess Data Quality:** Identifies sequencing errors, adapter contamination, and quality issues
- **Detect Problems Early:** Catches issues before expensive computational analysis
- **Quality Metrics Evaluated:**
  - Per-base sequence quality (should be ≥30 Phred score)
  - GC content distribution
  - Adapter contamination levels
  - Duplicate reads percentage

**Input:** FASTQ files  
**Output:** FastQC report (HTML) with quality plots

**Interpretation Example:**

```
✓ PASS: Per base quality score >30 (high quality)
✓ PASS: GC content consistent with organism (~47% for Drosophila)
⚠ WARNING: High duplicate content (25-30% normal for RNA-Seq)
✓ PASS: No adapter contamination detected
```

---

#### **Step 1c: Aggregate Reports with MultiQC**

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

### 2. Read Trimming

#### **Step 2: Cutadapt - Trim Adapters & Low-Quality Bases**

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

### 3. Read Mapping

#### **Step 3: STAR - Spliced Alignment to Reference Genome**

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

STAR (Spliced Transcripts Alignment to a Reference) is the **gold standard** for RNA-Seq mapping because:

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

<img width="1000" height="500" alt="mapping result igv" src="https://github.com/user-attachments/assets/eb1f2c15-f35f-4084-ba91-c354b1b18c82" />


---
### 4. Post-Mapping Quality Checks

#### **Step 4a: MarkDuplicates - Identify PCR Duplicates**

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

**Example Report:**
```
Sample: GSM461177_paired
Total reads: 10,455,687
Duplicate reads: 2,708,000 (25.9%)
Duplicate rate indicates: Normal (within expected range)
```

---

#### **Step 4b: RSeQC - Quality Control Suite**

**Step 4b-i: Gene Body Coverage**

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
✓ Expected: Relatively even coverage 5' to 3' (flat line)
✗ Concern: Sharp drop at 3' end → possible RNA degradation
✗ Concern: Drop at 5' end → unusual, may indicate technical issue
```

---

**Step 4b-ii: Infer Experiment - Strand Detection**

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

**Example Output:**
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

**Step 4b-iii: Read Distribution**

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

**Interpretation:**
- ✓ If >70% exonic: High-quality cDNA
- ⚠ If >20% intergenic: Possible genomic DNA contamination
- ⚠ If >25% intronic: May indicate pre-mRNA contamination

---

#### **Step 4c: Aggregate with MultiQC**

Combines all post-mapping QC results into one interactive report for easy visualization and comparison across samples.

---

### 5. Feature Counting

#### **Step 5: featureCounts - Quantify Gene Expression**

| Parameter | Value | Reasoning |
|-----------|-------|-----------|
| **Tool** | featureCounts (Subread) | Fast, accurate, industry-standard |
| **Input** | BAM files (with duplicates flagged) | Aligned reads |
| **Annotation** | GTF file | Gene/exon coordinates |
| **Count Method** | Per-gene | Aggregate reads across all exons |
| **Multi-mapped Reads** | Count once per location | Conservative quantification |
| **Duplicate Handling** | Skip duplicates | Remove PCR bias (adjustable parameter) |
| **Strandedness** | Based on RSeQC result | Use detected library type |
| **Overlap Ambiguity** | NOT overlap (IntersectionStrict) | Read must overlap exon only, not adjacent feature |
| **Min Overlap** | 1 bp | Any overlap counts |
| **Pair-end Mode** | Yes | Count fragments, not individual reads |

**Why featureCounts?**

1. **Speed:** 10x faster than HTSeq while maintaining accuracy
2. **Accuracy:** Handles multimapped reads, overlapping features intelligently
3. **Summarization:** Automatically aggregates exons into gene-level counts
4. **Standard Output:** Compatible with all downstream analysis tools

**Input:**
- BAM file (aligned, sorted, indexed)
- GTF annotation file
- Library strandedness information

**Output:**

| File Type | Content | Use Case |
|-----------|---------|----------|
| **counts.txt** | Gene-level read counts | Differential expression analysis |
| **counts.txt.summary** | Classification summary | QC and troubleshooting |
| **counts.txt.jcounts** | Junction-based counts (optional) | Isoform-level analysis |

**Output Matrix Format:**

```
Geneid          Length  GSM461176  GSM461177  GSM461179  GSM461180
FBgn0000003     1000    234        567        123        456
FBgn0000008     2500    890        1234       567        789
FBgn0000014     3200    1200       1567       234        567
...
```

**Understanding the Summary Report:**

```
Status                  GSM461176  GSM461177  GSM461179  GSM461180
Assigned               8,234,567  9,123,456  7,234,567  8,456,789
Unassigned_Ambiguous     234,567    345,678    123,456    234,567
Unassigned_NoFeatures    567,890    234,567    678,901    345,678
Unassigned_Unmapped            0          0          0          0
Unassigned_MappingQuality 123,456    123,456    123,456    123,456
```

| Category | Expected % | Interpretation |
|----------|-----------|----------------|
| **Assigned** | 70-85% | Reads mapping to known genes |
| **Unassigned_Ambiguous** | 2-5% | Reads overlapping multiple genes (filter) |
| **Unassigned_NoFeatures** | 10-20% | Unmapped or intergenic reads (quality issue if >25%) |

**Expected Read Count Statistics:**

| Sample | Total Assigned | Median Count | Max Count |
|--------|---|---|---|
| High quality | 8-12M | 5-20 | 50,000-200,000 |
| Medium quality | 5-8M | 2-10 | 20,000-100,000 |
| Low quality | <5M | <2 | <20,000 |

---

### 6. Differential Expression Analysis

#### **Step 6: DESeq2 - Statistical Differential Expression Testing**

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

**Quality Plots Interpretation:**

**MA-Plot:** (M = log2FC, A = mean log2 expression)
```
- Red dots = significant DE genes
- Expected: Points mostly centered on y=0 (no change)
- Problem: Red dots heavily skewed up/down = possible batch effect
```

**Volcano Plot:** (x-axis = log2FC, y-axis = -log10 p-value)
```
- Red dots = significant genes (high |FC| AND low p-value)
- Expected: Two "wings" of red dots (up-regulated on right, down-regulated on left)
- Problem: All points on vertical line = no fold-change effect (only noise)
```

---

### 7. Visualization & Enrichment Analysis

#### **Step 7a: Extract and Visualize Top DE Genes**

| Tool | Purpose | Output |
|------|---------|--------|
| **Filter Results** | Select genes by p-value/FC threshold | Reduced table of DE genes |
| **Sort & Rank** | Order by log2FC or adjusted p-value | Ranked gene list |
| **Heatmap** | Visualize expression patterns | Clustered expression visualization |
| **Export** | Save for external tools (Excel, GraphPad) | Tabular format |

**Heatmap Interpretation:**

```
       Untreated         Treated
Gene1  ████████ (high)   ██ (low)       → UP in treated
Gene2  ██ (low)          ████████ (high) → DOWN in treated
Gene3  ████ (medium)     ████ (medium)  → No change
```

---

#### **Step 7b: Functional Enrichment Analysis (Gene Ontology)**

| Tool | Purpose | Input |
|------|---------|-------|
| **g:Profiler / Ensembl** | Identify enriched biological processes, cellular components, molecular functions | List of DE gene IDs |

**Gene Ontology Categories:**

| Ontology | Examples | Interpretation |
|----------|----------|-----------------|
| **Biological Process (BP)** | "DNA replication", "immune response", "apoptosis" | What biological functions are affected? |
| **Cellular Component (CC)** | "nucleus", "mitochondria", "cell membrane" | Where in the cell are changes happening? |
| **Molecular Function (MF)** | "kinase activity", "ATP binding", "transcription factor" | What molecular activities changed? |

**Example Enrichment Results:**

```
Top enriched GO terms in UP-regulated genes:
1. "RNA splicing" (p = 1.2e-8)
   - 45 genes affected
   - Fold enrichment = 2.3x
   
2. "Nuclear transport" (p = 3.4e-5)
   - 23 genes affected
   - Fold enrichment = 1.8x
```

**Interpretation:**
- Pasilla gene (a splicing regulator) is depleted
- Enrichment of splicing-related genes makes biological sense (validation!)
- Shows DESeq2 results are biologically meaningful

---

## 🛠️ Tool Overview

| Tool | Function | Input | Output | Why Use? |
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

## ✅ Quality Control Checkpoints

### Checkpoint 1: Raw Data Quality

| Check | Target | Pass/Fail Criteria |
|-------|--------|-------------------|
| **Per-base Q-score** | ≥30 | Q-score at 5' and 3' both >28 |
| **GC Content** | ±5% of mean | Within normal range for organism |
| **Adapter Content** | <5% | Few reads contain adapters |
| **Duplicate Rate** | 25-50% | >50% warrants investigation |

### Checkpoint 2: After Trimming

| Check | Target | Pass/Fail Criteria |
|-------|--------|-------------------|
| **Length Distribution** | >20 bp | <5% discarded as too short |
| **Adapter Removal** | >99% | Minimal residual adapters |
| **Pair Retention** | >90% | Most read pairs retained |

### Checkpoint 3: After Mapping

| Check | Target | Pass/Fail Criteria |
|-------|--------|-------------------|
| **Mapping Rate** | >85% | ≥75% acceptable, <70% is concerning |
| **Uniquely Mapped** | >80% | High specificity |
| **Splice Detection** | 40-60% | Organism-dependent |
| **Unmapped <2%** | <5% | Low error rate |

### Checkpoint 4: Post-Mapping QC

| Check | Target | Pass/Fail Criteria |
|-------|--------|-------------------|
| **Duplicate Rate** | <50% | Indicates PCR bias if higher |
| **Gene Body Coverage** | Uniform | Flat 5'→3' distribution |
| **Strandedness** | >90% one pattern | Confirms library type |
| **Exonic Reads** | >70% | <70% indicates potential contamination |

### Checkpoint 5: Feature Counting

| Check | Target | Pass/Fail Criteria |
|-------|--------|-------------------|
| **Assigned Reads** | 70-85% | ≥60% required |
| **Ambiguous Reads** | <5% | Remove in downstream analysis if high |
| **Count Distribution** | Median >5 | Too many zeros (low expression) |

### Checkpoint 6: Differential Expression

| Check | Target | Pass/Fail Criteria |
|-------|--------|-------------------|
| **DE Gene % | 5-30% | <1% = no real effect; >50% = possible batch effect |
| **Volcano Plot** | Symmetric wings | Asymmetric = possible batch/bias |
| **MA-plot** | Centered at M=0 | Heavy skew = technical artifact |
| **GO Enrichment** | Biologically relevant | Enrichment should match biological hypothesis |

---

## 🎯 Parameters & Best Practices

### Key Parameter Decisions

#### **1. Junction Window Size (STAR)**
```
Recommended: Read_Length - 1
Example: For 36 bp reads → use 35
Reason: Enables detection of splice junctions at maximum sensitivity
```

#### **2. Duplicate Handling**
```
Approach 1: Mark duplicates, filter later
- Preserves reads for re-analysis with different parameters
- Recommended for initial exploratory analysis

Approach 2: Remove duplicates before counting
- Reduces PCR bias in differential expression
- Use only if duplicate rate > 50%
```

#### **3. Strandedness Parameter in featureCounts**
```
If unstranded library:
  -s 0 (count reads from both strands)
  
If reverse-stranded (common):
  -s 2 (specific strand counting)
  
If forward-stranded:
  -s 1 (specific strand counting)
```

#### **4. DESeq2 Design Formula**
```
Simple comparison: ~ condition
  Design = formula specifying Treated vs. Untreated

Batch effect correction: ~ batch + condition
  If samples from different sequencing runs, add as variable
```

#### **5. Fold Change Threshold**
```
Biological Context     Recommended Log2FC
─────────────────────────────────────────
General discovery      ≥ 1.0 (2-fold change)
Conservative (strict)  ≥ 2.0 (4-fold change)
Exploratory           ≥ 0.5 (1.4-fold change)
```

### Best Practices

1. **Always Use Replicates**
   - Minimum 3 per condition (ideally 5+)
   - Biological replicates > technical replicates

2. **Quality First**
   - Do not skip QC steps
   - Address quality issues before proceeding

3. **Reproducibility**
   - Document all parameters
   - Save workflow XML from Galaxy
   - Version control input files

4. **Multiple Comparisons Correction**
   - Use adjusted p-values (Benjamini-Hochberg FDR)
   - Never use raw p-values

5. **Validation**
   - Run independent validation (qRT-PCR) on top hits
   - Verify enrichment results match biological hypothesis

---

## 📊 Expected Results

### Sample Output Tables

**Differential Expression Results (Top 10 Up-regulated):**

| Gene ID | Gene Name | Basemean | Log2FC | p-value | padj | Regulation |
|---------|-----------|----------|--------|---------|------|-----------|
| FBgn0263574 | *nsr* | 2341 | +3.45 | 2.3e-8 | 0.001 | ↑↑ UP |
| FBgn0028502 | *CG7897* | 1234 | +3.12 | 5.6e-7 | 0.003 | ↑↑ UP |
| FBgn0003012 | *hsr* | 3456 | +2.89 | 1.2e-6 | 0.005 | ↑↑ UP |

**Expected DE Gene Statistics:**

```
Dataset: Pasilla gene depletion vs. Wild-type
Total genes analyzed: 13,602
Significantly DE genes (padj < 0.05): ~2,200 (16%)
  ├─ Up-regulated: ~1,234
  └─ Down-regulated: ~987
```

### Quality Metrics Summary

**Pre-processing Results:**

```
Input reads (total):        94,567,234
After Cutadapt trimming:    92,456,123 (97.8% retained)
After mapping to genome:    85,234,567 (87% mapped, 75.7M uniquely)
```

**Post-Mapping Metrics:**

```
Duplicate reads:            22,340,567 (25.9% of mapped) ✓
PCR duplicates flagged:     5,234,123 (Picard output)
Exonic reads:               73,456,234 (86.2%) ✓
Intronic reads:             7,234,512 (8.5%) ✓
Intergenic reads:           4,543,821 (5.3%) ✓
```

---

## 🔧 Troubleshooting

### Common Issues and Solutions

**Issue: Very Low Mapping Rate (<70%)**

| Possible Cause | Solution |
|---|---|
| Wrong reference genome | Verify organism and genome version match |
| Contamination | Check FastQC for adapter/rRNA contamination |
| Bad sequencing quality | Examine per-base quality scores; retrim if needed |
| Library prep problem | Review sample preparation protocol; may need new library |

**Issue: High Duplicate Rate (>60%)**

| Possible Cause | Solution |
|---|---|
| PCR over-amplification | Check library prep kit; lower PCR cycles |
| Very low input RNA | Higher duplicates expected with low complexity |
| Highly expressed genes | Normal for small genome studies; not a problem |
| Actual duplicates | Use MarkDuplicates, but don't remove for DE analysis |

**Issue: Few/No DE Genes (<1%)**

| Possible Cause | Solution |
|---|---|
| No biological effect | The treatment truly had minimal impact |
| Insufficient replicates | Add more biological replicates for power |
| Small fold changes | Loosen log2FC threshold (e.g., >0.5 instead of >1) |
| Batch effects | Check for sequencing run or sample prep batch effects |
| Wrong comparison | Verify sample condition assignments are correct |

**Issue: Too Many DE Genes (>50%)**

| Possible Cause | Solution |
|---|---|
| Batch effects | Look for non-biological patterns in PCA plot |
| Incorrect condition assignment | Verify sample metadata is correct |
| Technical artifacts | Check post-mapping QC (strand detection, coverage) |
| Very large effect | Possible, but verify with GO enrichment |

**Issue: Unequal Library Sizes**

| Symptom | Solution |
|---|---|
| Some samples 10x more/fewer reads | Normal variation; DESeq2 handles normalization |
| Systematic size difference by condition | Possible sequencing depth batch; DESeq2 normalizes |

**Issue: Failed Mapping Step**

| Error Message | Cause | Solution |
|---|---|---|
| "STAR index error" | Reference genome not built | Let Galaxy auto-build; may take hours |
| "Paired-end mismatch" | R1 and R2 have different read counts | Use Cutadapt output directly; trimming mismatches reads |
| "GTF format error" | Annotation file corrupted or wrong format | Verify GTF is gzipped; check format with preview |

---

## 📚 References

### Original Publication & Datasets

Brooks, D. G., et al. 2011. RNA-Seq Data Analysis in Galaxy. Dataset GSM461176-GSM461182 (Pasilla gene depletion in Drosophila).

### Tool Publications

- **STAR:** Ultrafast universal RNA-seq aligner. Dobin, A., et al. (2013). Bioinformatics, 29(1), 15-21.
- **featureCounts:** Subread, Bioconductor RNA-Seq analysis package.
- **DESeq2:** Statistical analysis of differential gene expression. Love, M. I., et al. (2014). Genome Biology, 15(12), 550.
- **RSeQC:** Quality assessment of RNA-Seq data. Wang, L., et al. (2012). Bioinformatics, 28(16), 2184-2185.
- **Cutadapt:** Adapter and quality-based trimming. Martin, M. (2011). EMBnet.journal, 17(1), 10-12.

### Galaxy Training Materials

- [Galaxy Training Network - Reference-based RNA-Seq](https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/ref-based/tutorial.html)
- [Galaxy for RNA-Seq Analysis - Methods in Molecular Biology (2021)](https://doi.org/10.1007/978-1-0716-1307-8_20)

### Recommended Reading

- Soneson, C., & Delorenzi, M. (2013). "A comparison of methods for differential expression analysis of RNA-seq data." PLoS ONE, 8(10), e78563.
- Dillies, M. A., et al. (2013). "A comprehensive evaluation of normalization methods for Illumina high-throughput RNA sequencing data analysis." Briefings in Bioinformatics, 14(6), 671-683.
- Auer, P. L., & Doerge, R. W. (2010). "Statistical design and analysis of RNA sequencing data." Genetics, 185(2), 405-416.

### Gene Ontology Resources

- Gene Ontology: https://geneontology.org/
- g:Profiler: https://biit.cs.ut.ee/gprofiler/gost

---

## 📋 Dataset Requirements & System Specs

### Input Requirements

- **Format:** FASTQ (gzip compressed)
- **Read Type:** Single-end or paired-end (36-150 bp typical)
- **Minimum Replicates:** 3 per condition (5+ recommended)
- **Library Type:** mRNA-Seq (typically Illumina TrueSeq or similar)

### Computational Resources

- **RAM:** 8+ GB (for Drosophila), 32+ GB (for human/mammalian)
- **Storage:** 1-5 GB per sample (FASTQ → BAM → counts)
- **Processing Time:** 2-8 hours per sample (depends on size and system)
- **No Installation Needed:** Galaxy provides all tools

---

## 🎓 Learning Outcomes

After completing this tutorial, you will be able to:

✓ Assess RNA-Seq data quality using Falco/FastQC and MultiQC  
✓ Trim adapters and low-quality bases appropriately  
✓ Map reads to a reference genome using STAR  
✓ Evaluate mapping quality using RSeQC and Picard tools  
✓ Quantify gene expression from BAM files using featureCounts  
✓ Perform statistical differential expression analysis with DESeq2  
✓ Interpret volcano plots, MA-plots, and heatmaps  
✓ Identify biologically relevant gene sets via GO enrichment  
✓ Troubleshoot common RNA-Seq analysis issues  
✓ Document reproducible workflows in Galaxy  

---

## 📄 Citation

If you use this tutorial or Galaxy for your research, please cite:

```
Batut, B., Hiltemann, S., Bagnacani, A., et al. (2018).
Community-driven data analysis training for biology.
Nature Biotechnology, 36(10), 943-945.

Galaxy Training Network. "Reference-based RNA-Seq data analysis."
Available at: https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/ref-based/tutorial.html
```

---

## 📞 Support & Community

- **Galaxy Help Forum:** https://help.galaxyproject.org/
- **Galaxy Training Chat:** https://gitter.im/Galaxy-Training-Network/
- **GTN GitHub Issues:** https://github.com/galaxyproject/training-material/issues

---

## 📜 License

This material is licensed under the **Creative Commons Attribution 4.0 International License** (CC BY 4.0).

You are free to:
- Share and adapt this material
- Use it for educational purposes
- Include it in your own projects

Provided you:
- Give appropriate credit
- Indicate if changes were made
- Include a link to the license

---

**Last Updated:** April 2026  
**Version:** 2.0  
**Status:** Complete & Tested ✓

---

## 🙏 Acknowledgments

This README is based on the Galaxy Training Network's "Reference-based RNA-Seq data analysis" tutorial, developed and maintained by the Galaxy community including Bérénice Batut, Marius van den Beek, Maria Doyle, and Nicola Soranzo.

Original dataset from Brooks et al. (2011) *Pasilla* gene depletion study in *Drosophila melanogaster*.

