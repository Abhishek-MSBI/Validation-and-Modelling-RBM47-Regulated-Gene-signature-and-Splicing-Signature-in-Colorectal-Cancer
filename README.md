# Validation-and-Modelling-RBM47-Regulated-Gene-signature-and-Splicing-Signature-in-Colorectal-Cancer
This repository contains code, data, and workflow documentation for the validation and modeling of RBM47-regulated gene and splicing signatures in colorectal cancer.

This analysis serves as the foundation for a project to develop and validate a novel prognostic biomarker for colorectal cancer.

## Project Goal

The project is structured in several phases:
1.  *Analyze* the transcriptomic impact of RBM47 overexpression in a CRC cell line.
2.  *Identify* a high-confidence "gene signature" (Differentially Expressed Genes) and "splicing signature" (Alternative Splicing events) regulated by RBM47.
3.  *Validate* these signatures against a large clinical cohort (TCGA-COAD) to assess their prognostic power.
4.  *Model* the validated signatures to build a minimal, predictive biomarker for patient survival.

This repository contains the full workflow and results for *Phases 1 and 2*.

## Analysis Workflow

1.  *RNA-Seq Data Processing:* Raw FASTQ files (GSE282110) were aligned to the reference genome (e.g., using HISAT2) to produce BAM alignment files.
2.  *Gene Quantification:* The alignment files were processed with featureCounts to generate the raw gene count matrix (RBM47_counts.txt).
3.  *DGE Analysis:* The raw count matrix was analyzed for differential expression (e.g., using DESeq2). The raw results are in RBM47_DGE_results_corrected.csv. The filter_dge.py script was used to filter these for significance (padj < 0.05, |log2FoldChange| > 1).
4.  *Alternative Splicing Analysis:* The same BAM alignment files were analyzed with rMATS to detect differential alternative splicing. The BASH script run_rmats.sh was used to automate this step.
5.  *AS Filtering:* The raw rMATS output was filtered for significant events (FDR < 0.05, |IncLevelDifference| > 0.1) using the filter_rmats.R script.

## Key Findings: The RBM47 Signature

This analysis successfully identified a clear, two-part signature regulated by RBM47:

*Gene Signature:* *[XX]* significant Differentially Expressed Genes (DEGs).
    (Run python filter_dge.py and put the number from the "Significant genes found" output here)
*Splicing Signature:* *2,699* significant Differential Alternative Splicing events.

The complete, filtered results for this signature are available in this repository.

## Repository Contents

*/rMATS_Significant_Results_Publication/*: Folder containing the 5 filtered CSV files of all significant splicing events.
*RBM47_significant_DGE_results.csv*: The filtered list of significant DEGs.
*RBM47_DGE_results_corrected.csv*: The raw, unfiltered DGE output.
*RBM47_counts.txt*: The raw gene count matrix.
*run_rmats.sh*: BASH script to run the rMATS analysis.
*filter_rmats.R*: R script to filter the rMATS results.
*filter_dge.py*: Python script to filter the DGE results.

## Next Steps

The gene and splicing signatures identified here will now be used for *Phase 3*: clinical validation and survival analysis using the TCGA-COAD patient cohort.
