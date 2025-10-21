# --- 1.1: Install BiocManager (if you don't have it) ---
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# --- 1.2: Install Bioconductor and CRAN packages ---
# This might take a few minutes
BiocManager::install("TCGAbiolinks")
install.packages(c("survminer", "survival"))

# --- 1.3: Load the libraries for this session ---
library(TCGAbiolinks)
library(SummarizedExperiment) # Will be loaded by TCGAbiolinks, but good to be explicit
library(survival)
library(survminer)

# --- 2.1: Query for the gene expression data ---
# We'll get Transcriptome Profiling (RNA-Seq)
# GDCquery will find all files that match these criteria
query_coad <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# --- 2.2: Download the data ---
# This is the BIG download. It will download all the raw files.
# It may take a very long time (30+ minutes) depending on your internet.
# It will create a "GDCdata" folder in your working directory.
tryCatch({
  GDCdownload(
    query = query_coad,
    method = "api",
    files.per.chunk = 100 # Helps with large downloads
  )
}, error = function(e) {
  message("Download failed, trying with 'client' method...")
  GDCdownload(
    query = query_coad,
    method = "client" # Fallback method
  )
})

# --- 2.3: Prepare the data into an R object ---
# This is a crucial step.
# TCGAbiolinks reads all those downloaded files and combines them
# into a single, powerful object: a `SummarizedExperiment`.
# This object contains:
#   1. The count matrix (assay)
#   2. The patient clinical data (colData)
#   3. The gene information (rowData)
# This step also takes a few minutes.
se_coad <- GDCprepare(query = query_coad)

# You can inspect the object:
print(se_coad)
# It will show you the dimensions (genes x patients)

# --- 3.1: Get the clinical data ---
# The clinical data is stored in the 'colData'.
# We can convert it to a standard data.frame.
clinical_df <- as.data.frame(colData(se_coad))

# --- 3.2: Create the survival variables (time and status) ---
# We need two columns:
# 1. 'time': Number of days to death or last follow-up.
# 2. 'status': The event (1 = Dead, 0 = Alive).

# Use 'days_to_death' if patient is dead, otherwise use 'days_to_last_follow_up'
clinical_df$survival_time <- ifelse(
  !is.na(clinical_df$days_to_death),
  clinical_df$days_to_death,
  clinical_df$days_to_last_follow_up
)

# Create the event status
clinical_df$vital_status_numeric <- ifelse(
  clinical_df$vital_status == "Dead",
  1,
  0
)

# --- 3.3: Get the gene expression data ---
# The raw counts are in the 'assay' part of the object.
# Rownames are Ensembl IDs (e.g., ENSG00000...).
count_matrix <- assay(se_coad, "unstranded") # or "tpm_unstranded" if available

if (!require("BiocManager", quietly=TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")


# --- 3.4: Normalize and filter for your gene (RBM47) ---
# We can't use raw counts. We'll do a simple log2-transform
# of "counts-per-million" (CPM) for visualization.
# For a more robust analysis, use vst() or voom() as planned for your model.
cpm_matrix <- edgeR::cpm(count_matrix)
log2_cpm_matrix <- log2(cpm_matrix + 1)

# The gene names (symbols) are in 'rowData'.
gene_info <- as.data.frame(rowData(se_coad))

# Find the Ensembl ID for RBM47
# (Replace 'RBM47' with any gene from your signature)
gene_of_interest <- "RBM47"
ensembl_id <- gene_info$gene_id[gene_info$gene_name == gene_of_interest]

# Get the log2-cpm expression values for just that one gene
# We use [1] in case the name maps to multiple Ensembl IDs
gene_expression <- log2_cpm_matrix[ensembl_id[1], ]

# --- 3.5: (Corrected) Create the final analysis data.frame ---

# --- Part A: Get Expression Data ---
# `gene_expression` is a named vector. The `names()` are the full barcodes.
expression_df <- data.frame(
  full_barcode = names(gene_expression),
  expression = as.numeric(gene_expression)
)

# --- Part B: Get Survival Data ---
# `clinical_df` is the data.frame from colData(se_coad).
# The `rownames(clinical_df)` are the full barcodes that match the expression data.
survival_df <- data.frame(
  full_barcode = rownames(clinical_df),
  survival_time = clinical_df$survival_time,
  vital_status_numeric = clinical_df$vital_status_numeric
)

# --- Part C: Merge ---
# NOW, we merge using 'full_barcode', which exists in both data frames.
analysis_df <- merge(survival_df, expression_df, by = "full_barcode")

# --- Part D: Clean up ---
# 1. Remove any patient with NA in time, status, or expression
analysis_df <- na.omit(analysis_df)

# 2. Remove any patient with a survival time of 0 or less (not valid)
analysis_df <- analysis_df[analysis_df$survival_time > 0, ]

# --- Part E: Check the data ---
# Let's run the checks again. This time you should see observations.
print("--- Checking Cleaned Data ---")
str(analysis_df)
summary(analysis_df)
print("--- Event Table (0=Alive, 1=Dead) ---")
table(analysis_df$vital_status_numeric, useNA = "ifany")
print("--- End of Check ---")

# --- 4.1: Find the optimal cutpoint ---
# This function tests all possible cutoffs and finds the one
# that gives the most significant (lowest p-value) separation.
cutpoint_object <- surv_cutpoint(
  data = analysis_df,
  time = "survival_time",
  event = "vital_status_numeric",
  variables = "expression" # The variable we're testing
)

# You can see the best cutoff value
print(cutpoint_object)

# --- 4.2: Categorize patients based on the cutpoint ---
# This adds a new column: 'High' or 'Low'
analysis_df <- surv_categorize(cutpoint_object)

# --- 5.1: Create the survival object ---
# This tells R which column is time and which is the event
surv_object <- Surv(
  time = analysis_df$survival_time,
  event = analysis_df$vital_status_numeric
)

# --- 5.2: Fit the survival model ---
# This fits the Kaplan-Meier model, comparing the groups
# (which we named 'expression' in step 4.2)
fit_model <- survfit(surv_object ~ expression, data = analysis_df)

# --- 5.3: Create the plot! ---
# ggsurvplot is a powerful function to make a publication-ready plot
ggsurvplot(
  fit_model,
  data = analysis_df,
  pval = TRUE,                 # Add the log-rank p-value
  conf.int = TRUE,             # Add confidence intervals
  risk.table = TRUE,           # Add the "at-risk" table below the plot
  legend.title = paste(gene_of_interest, "Expression"),
  legend.labs = c("High", "Low"),
  xlab = "Time (days)",
  ylab = "Overall Survival Probability",
  title = paste("Prognostic Value of", gene_of_interest, "in TCGA-COAD"),
  ggtheme = theme_classic()    # Use a clean theme
)
ggsave("RBM47_KM_Survival_Plot.png", width = 8, height = 6, dpi = 300)

# --- 1. Publication-Ready Plot (KM + All Tables) ---

# We will store the plot in a variable 'g' to customize it
g <- ggsurvplot(
  fit_model,
  data = analysis_df,
  
  # --- Main Plot Aesthetics ---
  title = paste("Prognostic Value of", gene_of_interest, "in TCGA-COAD"),
  xlab = "Time (days)",
  ylab = "Overall Survival Probability",
  
  # --- P-Value & Confidence Interval ---
  pval = TRUE,
  conf.int = TRUE,
  
  # --- Legend ---
  legend.title = paste(gene_of_interest, "Expression"),
  legend.labs = c("High", "Low"),
  
  # --- Risk Table ---
  risk.table = TRUE,
  tables.height = 0.2, # Adjust height of all tables
  tables.theme = theme_cleantable(), # A clean theme for the tables
  
  # --- NEW: Cumulative Events Table ---
  cumevents = TRUE,
  
  # --- NEW: Cumulative Censored Table ---
  cumcensor = TRUE,
  
  # --- Other Visual Tweaks ---
  ggtheme = theme_bw() # A clean publication-ready theme
)

# Print the plot
print


# --- 2. Cumulative Hazard Plot ---

ggsurvplot(
  fit_model,
  data = analysis_df,
  title = paste("Cumulative Hazard by", gene_of_interest, "Expression"),
  xlab = "Time (days)",
  ylab = "Cumulative Hazard",
  
  # --- The key change is this line: ---
  fun = "cumhaz",
  
  pval = TRUE,
  conf.int = TRUE,
  legend.title = paste(gene_of_interest, "Expression"),
  legend.labs = c("High", "Low"),
  ggtheme = theme_bw()
)

# --- 3. Cumulative Events Plot ---

ggsurvplot(
  fit_model,
  data = analysis_df,
  title = paste("Cumulative Events (Deaths) by", gene_of_interest, "Expression"),
  xlab = "Time (days)",
  ylab = "Cumulative Number of Events",
  
  # --- The key change is this line: ---
  fun = "event",
  
  legend.title = paste(gene_of_interest, "Expression"),
  legend.labs = c("High", "Low"),
  ggtheme = theme_bw()
)

# Make sure your 'gene_of_interest' variable is still set
# (If not, set it: gene_of_interest <- "RBM47")

# Define a filename based on your gene
csv_filename <- paste0(gene_of_interest, "_analysis_dataframe.csv")

# Save the data.frame as a CSV file
write.csv(
  analysis_df,
  file = csv_filename,
  row.names = FALSE  # We don't need the R-specific row numbers
)

print(paste("Successfully saved data to:", csv_filename))

# Define a filename
summary_filename <- paste0(gene_of_interest, "_statistical_summary.txt")

# Use sink() to divert all console output to a file
sink(summary_filename)

# --- Print all the important information ---
cat("====================================================\n")
cat("          STATISTICAL SUMMARY FOR:", gene_of_interest, "\n")
cat("====================================================\n\n")

cat("--- 1. Optimal Cutpoint Object ---\n")
print(cutpoint_object)
cat("\n\n")

cat("--- 2. Survival Fit Model (Log-Rank Test) ---\n")
# This print() is what shows the p-value
print(fit_model)
cat("\n\n")

cat("--- 3. Detailed Survival Model Summary (Event Tables) ---\n")
# This summary() shows the detailed tables (n.risk, n.events, etc.)
summary(fit_model)
cat("\n\n")

cat("--- 4. Final Cleaned Data Head ---\n")
print(head(analysis_df))

# --- Stop diverting the output ---
sink()

print(paste("Successfully saved summary to:", summary_filename))

# --- 3.1: Save the survfit model object ---
model_filename <- paste0(gene_of_interest, "_survfit_model.Rds")
saveRDS(fit_model, file = model_filename)

# --- 3.2: Save the cutpoint object ---
cutpoint_filename <- paste0(gene_of_interest, "_cutpoint_object.Rds")
saveRDS(cutpoint_object, file = cutpoint_filename)

# --- 3.3: Save the final ggsurvplot list-object (optional, but good) ---
# (Assuming your plot is stored in a variable 'g' from the last step)
if (exists("g")) {
  plot_obj_filename <- paste0(gene_of_interest, "_ggsurvplot_object.Rds")
  saveRDS(g, file = plot_obj_filename)
}

print("Successfully saved all .Rds model files.")

# --- 1. Define your signature ---
# We'll use the genes we've confirmed and a few more from your rMATS file.
# (You can add more gene symbols to this list!)
gene_signature_list <- c("RBM47", "WIF1", "FANCI", "FAAP20", "HLA-F-AS1")

# --- 2. Get the expression data for ALL these genes ---
# (We assume 'gene_info' and 'log2_cpm_matrix' are still loaded)

# Find the Ensembl IDs for all our genes
ensembl_ids <- gene_info$gene_id[gene_info$gene_name %in% gene_signature_list]

# Filter the big expression matrix to *only* our signature genes
signature_expression_matrix <- log2_cpm_matrix[ensembl_ids, ]

# --- 3. Clean up the matrix ---
# We need gene symbols as rownames, not Ensembl IDs.
# Create a small mapping data.frame
gene_map <- data.frame(
  id = gene_info$gene_id,
  symbol = gene_info$gene_name
)
# Keep only the genes in our matrix
gene_map <- gene_map[gene_map$id %in% rownames(signature_expression_matrix), ]
# Remove any duplicate gene symbols, keeping the first
gene_map <- gene_map[!duplicated(gene_map$symbol), ]

# Match the matrix rows to the map and assign new names
signature_expression_matrix <- signature_expression_matrix[gene_map$id, ]
rownames(signature_expression_matrix) <- gene_map$symbol

# --- 4. Transpose the matrix ---
# For machine learning, we need (patients x genes)
# (Currently, it is genes x patients)
ml_expression_df <- as.data.frame(t(signature_expression_matrix))

# Add a 'full_barcode' column for merging
ml_expression_df$full_barcode <- rownames(ml_expression_df)

# --- 5. Get the Survival Data ---
# (This object 'clinical_df' should still be in your R environment)
survival_df <- data.frame(
  full_barcode = rownames(clinical_df),
  survival_time = clinical_df$survival_time,
  vital_status_numeric = clinical_df$vital_status_numeric
)

# --- 6. Merge ML Expression Data + Survival Data ---
final_ml_df <- merge(survival_df, ml_expression_df, by = "full_barcode")

# --- 7. Final Clean ---
final_ml_df <- na.omit(final_ml_df)
final_ml_df <- final_ml_df[final_ml_df$survival_time > 0, ]

# --- 8. Inspect your final ML data! ---
print("--- Final Machine Learning Dataframe is Ready ---")
# This will show 'survival_time', 'vital_status_numeric', 'RBM47', 'WIF1', etc.
print(head(final_ml_df))





# --- START: FINAL CONSOLIDATED SAVING SCRIPT ---
# --- Add this entire block to the very end of your script ---

cat("\n\n--- STARTING FINAL SAVE PROCESS ---\n")

# --- 1. Save the FINAL Machine Learning Data (for Python) ---
# This is the most important part for our next step (LASSO-Cox).
# We will split your 'final_ml_df' into the two files we discussed.

if (exists("final_ml_df")) {
  # File 1: The Survival Data (y)
  survival_data_y <- final_ml_df[, c("full_barcode", "survival_time", "vital_status_numeric")]
  write.csv(survival_data_y, 
            file = "tcga_survival_data_y.csv", 
            row.names = FALSE)
  
  # File 2: The Feature Matrix (X)
  # This selects 'full_barcode' and all genes in your signature list
  feature_columns <- c("full_barcode", gene_signature_list) 
  feature_matrix_X <- final_ml_df[, feature_columns]
  write.csv(feature_matrix_X, 
            file = "tcga_feature_matrix_X.csv", 
            row.names = FALSE)
  
  cat("SUCCESS: Saved 'tcga_survival_data_y.csv' and 'tcga_feature_matrix_X.csv'.\n")
} else {
  cat("ERROR: 'final_ml_df' was not found. Could not save ML data files.\n")
}


# --- 2. Save the Core "se_coad" Object (The Time-Saver) ---
# This file is BIG, but it contains all the downloaded/prepared data.
# Next time, you can skip the GDCquery/GDCdownload/GDCprepare steps
# and just run: se_coad <- readRDS("se_coad_prepared.Rds")
if (exists("se_coad")) {
  saveRDS(se_coad, file = "se_coad_prepared.Rds")
  cat("SUCCESS: Saved main 'se_coad_prepared.Rds' object.\n")
} else {
  cat("WARNING: 'se_coad' object not found. Skipping save.\n")
}


# --- 3. Save the Missing Plots ---
# Your script generated these plots but didn't save them to files.
if (exists("fit_model") & exists("analysis_df") & exists("gene_of_interest")) {
  
  # Save the Cumulative Hazard plot
  plot_cumhaz <- ggsurvplot(
    fit_model, data = analysis_df, fun = "cumhaz", pval = TRUE,
    title = paste("Cumulative Hazard by", gene_of_interest, "Expression"),
    xlab = "Time (days)", ylab = "Cumulative Hazard",
    legend.title = paste(gene_of_interest, "Expression"),
    legend.labs = c("High", "Low"), ggtheme = theme_bw()
  )
  ggsave(paste0(gene_of_interest, "_Cumulative_Hazard_Plot.png"), 
         plot = print(plot_cumhaz), width = 8, height = 6, dpi = 300)
  
  # Save the Cumulative Events plot
  plot_event <- ggsurvplot(
    fit_model, data = analysis_df, fun = "event",
    title = paste("Cumulative Events by", gene_of_interest, "Expression"),
    xlab = "Time (days)", ylab = "Cumulative Number of Events",
    legend.title = paste(gene_of_interest, "Expression"),
    legend.labs = c("High", "Low"), ggtheme = theme_bw()
  )
  ggsave(paste0(gene_of_interest, "_Cumulative_Events_Plot.png"), 
         plot = print(plot_event), width = 8, height = 6, dpi = 300)
  
  cat("SUCCESS: Saved 'Cumulative_Hazard_Plot.png' and 'Cumulative_Events_Plot.png'.\n")
} else {
  cat("WARNING: 'fit_model' or 'analysis_df' not found. Skipping plot saves.\n")
}

cat("--- FINAL SAVE PROCESS COMPLETE ---\n")

# --- END: FINAL CONSOLIDATED SAVING SCRIPT ---