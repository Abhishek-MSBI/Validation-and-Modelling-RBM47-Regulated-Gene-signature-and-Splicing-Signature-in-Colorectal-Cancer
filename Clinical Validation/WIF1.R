# --- Load the main library ---
library(TCGAbiolinks)

# --- 1. Re-create the query (this is fast) ---
print("--- Step 1: Re-creating query ---")
query_coad <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts" 
  # Note: If you used "HTSeq - Counts" before, change it here.
)

# --- 2. Re-prepare the se_coad object (This may take 5-10 minutes) ---
print("--- Step 2: Re-creating se_coad object... (This is not downloading) ---")
# This will read from your 'GDCdata' folder, NOT re-download.
se_coad <- GDCprepare(query = query_coad)

print("--- 'se_coad' is back in your environment! ---")

# --- 3. Rebuild all helper objects ---
print("--- Step 3: Rebuilding helper objects... ---")
library(SummarizedExperiment)
library(edgeR)

# Recreate 'clinical_df' (with survival data)
clinical_df <- as.data.frame(colData(se_coad))
clinical_df$survival_time <- ifelse(
  !is.na(clinical_df$days_to_death),
  clinical_df$days_to_death,
  clinical_df$days_to_last_follow_up
)
clinical_df$vital_status_numeric <- ifelse(
  clinical_df$vital_status == "Dead",
  1,
  0
)

# Recreate 'gene_info'
gene_info <- as.data.frame(rowData(se_coad))

# Recreate 'log2_cpm_matrix'
# Check the assay name. "unstranded" is common for STAR - Counts
count_matrix <- assay(se_coad, "unstranded") 
cpm_matrix <- edgeR::cpm(count_matrix)
log2_cpm_matrix <- log2(cpm_matrix + 1)

print("--- Environment is READY. You can now run your function. ---")

library(survival)
library(survminer)

run_km_analysis("WIF1")

# --- Corrected Reusable KM Analysis Function ---

run_km_analysis <- function(gene_symbol) {
  
  print(paste("--- Starting Analysis for:", gene_symbol, "---"))
  
  # --- 1. Find the Gene ---
  ensembl_id <- gene_info$gene_id[gene_info$gene_name == gene_symbol]
  
  if (length(ensembl_id) == 0) {
    print(paste("ERROR: Gene symbol", gene_symbol, "not found. Skipping."))
    return(NULL)
  }
  
  # Get expression data (use [1] in case of multiple IDs)
  gene_expression <- log2_cpm_matrix[ensembl_id[1], ]
  
  # --- 2. Create Analysis DataFrame ---
  expression_df <- data.frame(
    full_barcode = names(gene_expression),
    expression = as.numeric(gene_expression)
  )
  
  survival_df <- data.frame(
    full_barcode = rownames(clinical_df),
    survival_time = clinical_df$survival_time,
    vital_status_numeric = clinical_df$vital_status_numeric
  )
  
  analysis_df <- merge(survival_df, expression_df, by = "full_barcode")
  
  # --- 3. Clean Data ---
  analysis_df <- na.omit(analysis_df)
  analysis_df <- analysis_df[analysis_df$survival_time > 0, ]
  
  if(nrow(analysis_df) == 0 || length(unique(analysis_df$vital_status_numeric)) < 2) {
    print(paste("ERROR: Not enough data or only one event type for", gene_symbol, ". Skipping."))
    return(NULL)
  }
  
  # --- 4. Find Cutpoint ---
  cutpoint_object <- surv_cutpoint(
    data = analysis_df,
    time = "survival_time",
    event = "vital_status_numeric",
    variables = "expression"
  )
  
  analysis_df <- surv_categorize(cutpoint_object)
  
  # --- 5. (CORRECTED) Fit & Plot Model ---
  # This is the standard, robust way to call survfit.
  # We build the Surv() object INSIDE the formula.
  fit_model <- survfit(
    Surv(survival_time, vital_status_numeric) ~ expression, 
    data = analysis_df
  )
  
  # --- 6. Create the Plot ---
  km_plot <- ggsurvplot(
    fit_model,
    data = analysis_df,
    pval = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    legend.title = paste(gene_symbol, "Expression"),
    legend.labs = c("High", "Low"),
    xlab = "Time (days)",
    ylab = "Overall Survival Probability",
    title = paste("Prognostic Value of", gene_symbol, "in TCGA-COAD"),
    ggtheme = theme_classic()
  )
  
  # Print the plot to your RStudio "Plots" pane
  print(km_plot)
  
  # --- 7. Save Files ---
  # Save the plot
  ggsave(paste0(gene_symbol, "_KM_Plot.png"), width = 8, height = 8, dpi = 300)
  
  # Save the data
  write.csv(analysis_df, paste0(gene_symbol, "_analysis_data.csv"), row.names = FALSE)
  
  print(paste("--- Analysis for", gene_symbol, "Complete. Plot and data saved. ---"))
  
  # Return the plot object
  return(km_plot)
}
run_km_analysis("WIF1")