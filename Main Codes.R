####  ---------------------------------------------------------------------------
####   The code block below performs Cox proportional hazards 
####   regression on proteins for early and late atrial fibrillation (AF) cohorts 
####   using Bonferroni correction for multiple testing. It calculates hazard ratios (HR), 
####   confidence intervals, and p-values for each protein in both cohorts, adjusting for 
####   clinical factors. The significant proteins (Bonferroni-adjusted) are saved for 
####   further Reactome, GO, and KEGG enrichment analysis. 
####  ---------------------------------------------------------------------------

library(survival)
library(dplyr)

# Load data for early and late cohorts
Early_cohort <- readRDS("~/Early-onset AF/Data/Early_cohort.RDS")
Late_cohort <- readRDS("~/Early-onset AF/Data/Late_cohort.RDS")

# Define clinical factors
clinical_factors <- c("Age_at_recruitment", "Gender_Male1_Female0", 
                      "Ethnic_white", "Education_level_numeric", "Socioeconomic_deprivation", "Current_smoking",
                      "Alcohol_status", "Severe_pollution", "physical_inactivity", 
                      "No_coffee_intake", "Loneliness", "Broad_depression", "hypertension", 
                      "Low_LDL_C", "Low_Triglycerides", "Diabetes", "Overweight", "Elevated_CRP", 
                      "Clinical_comorbidities")

# Load filtered protein data, excluding the first column which is 'Participant_ID'
filtered_protein_data <- readRDS("~/Early-onset AF/Data/filtered_protein_data.RDS")
Proteins <- colnames(filtered_protein_data)[-1]

# Create an empty data frame to store results
results <- data.frame(Protein = character(),
                      HR = numeric(),
                      CI_lower = numeric(),
                      CI_upper = numeric(),
                      P_value = numeric(),
                      Bonferroni_significance = logical(),
                      stringsAsFactors = FALSE)

# Define total number of tests for Bonferroni correction
total_tests <- length(Proteins)
alpha_bonferroni <- 0.05 / total_tests

# Loop over each protein and perform Cox regression for Early cohort
for (protein in Proteins) {
  
  # Create formula for the Cox regression
  formula <- as.formula(paste("Surv(Follow_up_duration_years, AF) ~", 
                              paste(c(clinical_factors, protein), collapse = " + ")))
  
  # Fit the Cox model
  cox_model <- coxph(formula, data = Early_cohort)
  
  # Extract the summary of the Cox model
  cox_summary <- summary(cox_model)
  
  # Extract HR, confidence interval, and p-value for the protein
  hr <- cox_summary$coefficients[protein, "exp(coef)"]
  ci_lower <- cox_summary$conf.int[protein, "lower .95"]
  ci_upper <- cox_summary$conf.int[protein, "upper .95"]
  p_value <- cox_summary$coefficients[protein, "Pr(>|z|)"]
  
  # Apply Bonferroni correction
  significant <- p_value < alpha_bonferroni
  
  # Append results to the results data frame
  results <- rbind(results, data.frame(Protein = protein,
                                       HR = hr,
                                       CI_lower = ci_lower,
                                       CI_upper = ci_upper,
                                       P_value = p_value,
                                       Bonferroni_significance = significant))
}

# Save results for Early cohort
saveRDS(results, file = "~/Early-onset AF/Data/Cox_Protein_Early.RDS")
write.csv(results, "~/Early-onset AF/Data/Cox_Protein_Early.csv", row.names = FALSE)

# Load and check the saved Early cohort results
Cox_Protein_Early <- readRDS("~/Early-onset AF/Data/Cox_Protein_Early.RDS")
table(Cox_Protein_Early$Bonferroni_significance)

#-------------------------------  Late cohort  ----------------------------------------

# Recreate an empty data frame to store results for Late cohort
results <- data.frame(Protein = character(),
                      HR = numeric(),
                      CI_lower = numeric(),
                      CI_upper = numeric(),
                      P_value = numeric(),
                      Bonferroni_significance = logical(),
                      stringsAsFactors = FALSE)

# Loop over each protein and perform Cox regression for Late cohort
for (protein in Proteins) {
  
  # Create formula for the Cox regression
  formula <- as.formula(paste("Surv(Follow_up_duration_years, AF) ~", 
                              paste(c(clinical_factors, protein), collapse = " + ")))
  
  # Fit the Cox model
  cox_model <- coxph(formula, data = Late_cohort)
  
  # Extract the summary of the Cox model
  cox_summary <- summary(cox_model)
  
  # Extract HR, confidence interval, and p-value for the protein
  hr <- cox_summary$coefficients[protein, "exp(coef)"]
  ci_lower <- cox_summary$conf.int[protein, "lower .95"]
  ci_upper <- cox_summary$conf.int[protein, "upper .95"]
  p_value <- cox_summary$coefficients[protein, "Pr(>|z|)"]
  
  # Apply Bonferroni correction
  significant <- p_value < alpha_bonferroni
  
  # Append results to the results data frame
  results <- rbind(results, data.frame(Protein = protein,
                                       HR = hr,
                                       CI_lower = ci_lower,
                                       CI_upper = ci_upper,
                                       P_value = p_value,
                                       Bonferroni_significance = significant))
}

# Save results for Late cohort
saveRDS(results, file = "~/Early-onset AF/Data/Cox_Protein_Late.RDS")
write.csv(results, "~/Early-onset AF/Data/Cox_Protein_Late.csv", row.names = FALSE)

# Load and check the saved Late cohort results
Cox_Protein_Late <- readRDS("~/Early-onset AF/Data/Cox_Protein_Late.RDS")
table(Cox_Protein_Late$Bonferroni_significance)

####  ---------------------------------------------------------------------------
####   The code block below performs Reactome, GO, and KEGG enrichment analysis on 
####   Bonferroni-significant proteins from Cox regression for early and late 
####   atrial fibrillation cohorts. It maps protein symbols to Entrez gene IDs, 
####   conducts the enrichment analysis, and saves the results. 
####   A custom plotting function then visualizes the top enriched pathways, 
####   presenting results from Reactome, GO (BP, MF, CC), and KEGG 
####   with a sophisticated color scheme for clarity and contrast.
####  ---------------------------------------------------------------------------

# Load necessary libraries
library(dplyr)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)  # Human genome annotation package
library(ReactomePA)

# Load Cox regression results for early and late cohorts
Cox_Protein_Early <- readRDS("~/Early-onset AF/Data/Cox_Protein_Early.RDS")
Cox_Protein_Late <- readRDS("~/Early-onset AF/Data/Cox_Protein_Late.RDS")

# 1. Filter Bonferroni_significant proteins in Early cohort
significant_early <- Cox_Protein_Early %>%
  filter(Bonferroni_significance == TRUE) %>%  # Keep only significant proteins
  dplyr::select(Protein)  # Select the Protein column

# 2. Filter Bonferroni_significant proteins in Late cohort
significant_late <- Cox_Protein_Late %>%
  filter(Bonferroni_significance == TRUE) %>%  # Keep only significant proteins
  dplyr::select(Protein)  # Select the Protein column

# Map protein symbols to Entrez IDs for enrichment analysis
gene_entrez <- mapIds(org.Hs.eg.db, keys = significant_early$Protein, 
                      column = "ENTREZID", keytype = "SYMBOL", 
                      multiVals = "first")
gene_entrez_early <- data.frame(SYMBOL = names(gene_entrez), ENTREZID = gene_entrez, stringsAsFactors = FALSE)

gene_entrez <- mapIds(org.Hs.eg.db, keys = significant_late$Protein, 
                      column = "ENTREZID", keytype = "SYMBOL", 
                      multiVals = "first")
gene_entrez_late <- data.frame(SYMBOL = names(gene_entrez), ENTREZID = gene_entrez, stringsAsFactors = FALSE)

# Load merged protein data and map Entrez IDs for all proteins
merged_protein_data <- readRDS("~/房颤细胞层面的危险因素/Data/merged_protein_data.RDS")
Proteins <- colnames(merged_protein_data)[-1]
length(Proteins)

gene_entrez_universe <- mapIds(org.Hs.eg.db, keys = Proteins, 
                               column = "ENTREZID", keytype = "SYMBOL", 
                               multiVals = "first")
gene_entrez_universe <- data.frame(SYMBOL = names(gene_entrez_universe), ENTREZID = gene_entrez_universe, stringsAsFactors = FALSE)

# --------------------------- Reactome Analysis -----------------------------------------
Reactome_early <- enrichPathway(gene = gene_entrez_early$ENTREZID, pvalueCutoff = 0.05, readable = TRUE,
                                universe = gene_entrez_universe$ENTREZID, pAdjustMethod = "fdr")
saveRDS(Reactome_early, file = "~/Early-onset AF/Data/Reactome_early.RDS")

Reactome_late <- enrichPathway(gene = gene_entrez_late$ENTREZID, pvalueCutoff = 0.05, readable = TRUE,
                               universe = gene_entrez_universe$ENTREZID, pAdjustMethod = "fdr")
saveRDS(Reactome_late, file = "~/Early-onset AF/Data/Reactome_late.RDS")

# --------------------------- GO Analysis -----------------------------------------
ego_early <- enrichGO(gene = gene_entrez_early$ENTREZID,
                      universe = gene_entrez_universe$ENTREZID,
                      OrgDb = org.Hs.eg.db,
                      ont = "ALL",
                      pAdjustMethod = "fdr",
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.05,
                      readable = TRUE)
saveRDS(ego_early, file = "~/Early-onset AF/Data/ego_early.RDS")

ego_late <- enrichGO(gene = gene_entrez_late$ENTREZID,
                     universe = gene_entrez_universe$ENTREZID,
                     OrgDb = org.Hs.eg.db,
                     ont = "ALL",
                     pAdjustMethod = "fdr",
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.05,
                     readable = TRUE)
saveRDS(ego_late, file = "~/Early-onset AF/Data/ego_late.RDS")

#-------------------------- KEGG Analysis ---------------------------------------------
kk_early <- enrichKEGG(gene = gene_entrez_early$ENTREZID,
                       organism = 'hsa',
                       universe = gene_entrez_universe$ENTREZID,
                       pAdjustMethod = "fdr",
                       pvalueCutoff = 0.05)
saveRDS(kk_early, file = "~/Early-onset AF/Data/kk_early.RDS")

kk_late <- enrichKEGG(gene = gene_entrez_late$ENTREZID,
                      organism = 'hsa',
                      universe = gene_entrez_universe$ENTREZID,
                      pAdjustMethod = "fdr",
                      pvalueCutoff = 0.05)
saveRDS(kk_late, file = "~/Early-onset AF/Data/kk_late.RDS")

#----------------------------- Plotting ---------------------------------
# Load saved data
Reactome_early <- readRDS("~/Early-onset AF/Data/Reactome_early.RDS")
ego_early <- readRDS("~/Early-onset AF/Data/ego_early.RDS")
kk_early <- readRDS("~/Early-onset AF/Data/kk_early.RDS")

# Function to create enrichment plots
create_enrichment_plot <- function(Reactome_early, ego_early, kk_early, top_n = 6) {
  reactome_df <- as.data.frame(Reactome_early@result)[1:top_n, ]
  go_df <- as.data.frame(ego_early@result)
  kegg_df <- as.data.frame(kk_early@result)[1:top_n, ]
  
  go_mf_df <- go_df %>% filter(ONTOLOGY == "MF") %>% arrange(p.adjust) %>% head(top_n)  
  go_bp_df <- go_df %>% filter(ONTOLOGY == "BP") %>% arrange(p.adjust) %>% head(top_n)
  go_cc_df <- go_df %>% filter(ONTOLOGY == "CC") %>% arrange(p.adjust) %>% head(top_n)
  
  reactome_df$source <- "Reactome"
  go_bp_df$source <- "GO BP"
  go_mf_df$source <- "GO MF"
  go_cc_df$source <- "GO CC"
  kegg_df$source <- "KEGG"
  
  all_data <- bind_rows(reactome_df, go_bp_df, go_mf_df, go_cc_df, kegg_df) %>%
    mutate(logP = -log10(p.adjust))  
  
  all_data <- all_data %>% arrange(source, p.adjust)
  
  plot <- ggplot(all_data, aes(x = factor(Description, levels = Description), y = logP, fill = source)) +
    geom_bar(stat = "identity", width = 0.7) +  
    theme_minimal() +
    labs(title = "Enrichment Analysis",
         x = "Pathway",
         y = "-log10(FDR-adjusted p-value)") +
    scale_fill_manual(values = c("GO BP" = "#4c566a", 
                                 "GO MF" = "#88c0d0", 
                                 "GO CC" = "#a3be8c", 
                                 "KEGG" = "#d08770",   
                                 "Reactome" = "#b48ead")) +  
    geom_text(aes(label = Count), hjust = 0.5, vjust = -1, color = "black", family = "Helvetica") +  
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  
    theme(text = element_text(family = "Helvetica", face = "bold", color = "black"),  
          axis.text.x = element_text(angle = 45, hjust = 1, family = "Helvetica", face = "bold", color = "black"),  
          axis.text.y = element_text(size = 9, face = "bold", family = "Helvetica", color = "black"),  
          panel.grid.minor = element_blank(),
          panel.spacing = unit(0.5, "lines"),
          legend.position = "top") +
    ylim(0, max(all_data$logP) * 1.2)  
  
  return(plot)
}

create_enrichment_plot(Reactome_early, ego_early, kk_early)

####  ---------------------------------------------------------------------------
####   The code block below performs machine learning model training, tuning, 
####   and evaluation for atrial fibrillation (AF) prediction. It uses LightGBM 
####   as the classifier, with significant proteins identified via Cox regression 
####   as features. The dataset is split into an 80% training-validation set and 
####   a 20% independent external validation set. The model is tuned using cross-validation 
####   and random search to find optimal hyperparameters, and its performance is evaluated 
####   on the external test set. AUC and feature importance are also computed for model 
####   interpretation.
####  ---------------------------------------------------------------------------
# 1. Load necessary packages
library(mlr3)
library(mlr3tuning)
library(mlr3learners)
# remotes::install_github("mlr-org/mlr3extralearners@*release")
library(mlr3extralearners)
library(data.table)
library(lightgbm)
library(pROC)
library(dplyr)
library(survival)
library(paradox)

# 2. Load data
new_AF_data <- readRDS("~/Early-onset AF/Data/new_AF_data.RDS")
Early_cohort <- readRDS("~/Early-onset AF/Data/Early_cohort.RDS")
Cox_Protein_Early <- readRDS("~/Early-onset AF/Data/Cox_Protein_Early.RDS")

# Select significant proteins
significant_early <- Cox_Protein_Early %>%
  filter(Bonferroni_significance == TRUE) %>%
  dplyr::select(Protein)

# Match date data
matching_data <- new_AF_data[new_AF_data$`Participant ID` %in% Early_cohort$Participant_ID, ]
matching_dates <- matching_data[, c("Participant ID", "Date of consenting to join UK Biobank")]
colnames(matching_dates) <- c("Participant_ID", "Date_consenting_to_join")
Early_cohort <- merge(Early_cohort, matching_dates, by.x = "Participant_ID", by.y = "Participant_ID", all.x = TRUE)

# Ensure Date_consenting_to_join is in date format
Early_cohort$Date_consenting_to_join <- as.Date(Early_cohort$Date_consenting_to_join)
# Sort by date
Early_cohort <- Early_cohort[order(Early_cohort$Date_consenting_to_join), ]

# Calculate the cutoff point for splitting: first 80% as training+validation, last 20% as independent external validation
cutoff_index <- floor(0.8 * nrow(Early_cohort))

# Split the dataset
train_val_data <- Early_cohort[1:cutoff_index, ]  # First 80% as training+validation
external_test_data <- Early_cohort[(cutoff_index + 1):nrow(Early_cohort), ]  # Last 20% as external validation

# 3. Extract features and target variable
predictors <- significant_early$Protein  # Extract protein features

# Training and external test set
X_train <- train_val_data[, predictors]
y_train <- train_val_data$AF

X_external_test <- external_test_data[, predictors]
y_external_test <- external_test_data$AF

# Ensure AF column is factor type
train_val_data$AF <- as.factor(train_val_data$AF)
external_test_data$AF <- as.factor(external_test_data$AF)

# 4. Create tasks
task_train <- TaskClassif$new("AF_train", backend = train_val_data[,c(predictors,"AF")], target = "AF", positive = "1")
task_test <- TaskClassif$new("AF_external_test", backend = external_test_data[,c(predictors,"AF")], target = "AF", positive = "1")

# 5. Set up LightGBM learner and parameter space
lrn_lgb <- lrn("classif.lightgbm", predict_type = "prob")

param_set <- ps(
  num_leaves = p_int(lower = 15, upper = 63),
  learning_rate = p_dbl(lower = 0.01, upper = 0.1),
  feature_fraction = p_dbl(lower = 0.7, upper = 0.9),
  bagging_fraction = p_dbl(lower = 0.7, upper = 0.9)
)

resampling <- rsmp("cv", folds = 5)
terminator <- trm("evals", n_evals = 30)

# Use TuningInstanceBatchSingleCrit to create tuning instance
instance <- TuningInstanceBatchSingleCrit$new(
  task = task_train,
  learner = lrn_lgb,
  resampling = resampling,
  measure = msr("classif.auc"),
  search_space = param_set,
  terminator = terminator
)

# Tuning
tuner <- tnr("random_search")
tuner$optimize(instance)

# Check best parameters
best_params <- instance$result
print(best_params)

saveRDS(instance, file = "~/Early-onset AF/Data/instance.RDS")

# 9. Train model using best parameters
lrn_lgb$param_set$values <- instance$result_learner_param_vals
lrn_lgb$train(task_train)

saveRDS(lrn_lgb, file = "~/Early-onset AF/Data/lrn_lgb.RDS")

lrn_lgb <- readRDS("~/Early-onset AF/Data/lrn_lgb.RDS")

# 10. Predict on the external validation set and evaluate
prediction <- lrn_lgb$predict(task_test)

# 11. Calculate AUC
roc_obj_external <- roc(y_external_test, prediction$prob[, 2])  # Take the second column of probability values
auc_value_external <- auc(roc_obj_external)
print(paste("External Validation AUC: ", auc_value_external))

# 12. Calculate AUC 95% confidence interval
ci_auc <- ci.auc(roc_obj_external)
print(paste("AUC 95% CI: ", ci_auc[1], "-", ci_auc[3]))

# 13. Feature importance
importance <- lrn_lgb$importance()
print(importance)

####  ---------------------------------------------------------------------------
####   The code block below performs machine learning model training, tuning, 
####   and evaluation for atrial fibrillation (AF) prediction. It uses LightGBM 
####   as the classifier, with significant proteins identified via Cox regression 
####   as features. The dataset is split into an 80% training-validation set and 
####   a 20% independent external validation set. The model is tuned using cross-validation 
####   and random search to find optimal hyperparameters, and its performance is evaluated 
####   on the external test set. AUC and feature importance are also computed for model 
####   interpretation.
####  ---------------------------------------------------------------------------
# 1. Load necessary packages
library(mlr3)
library(mlr3tuning)
library(mlr3learners)
# remotes::install_github("mlr-org/mlr3extralearners@*release")
library(mlr3extralearners)
library(data.table)
library(lightgbm)
library(pROC)
library(dplyr)
library(survival)
library(paradox)

# 2. Load data
new_AF_data <- readRDS("~/Early-onset AF/Data/new_AF_data.RDS")
Early_cohort <- readRDS("~/Early-onset AF/Data/Early_cohort.RDS")
Cox_Protein_Early <- readRDS("~/Early-onset AF/Data/Cox_Protein_Early.RDS")

# Select significant proteins
significant_early <- Cox_Protein_Early %>%
  filter(Bonferroni_significance == TRUE) %>%
  dplyr::select(Protein)

# Match date data
matching_data <- new_AF_data[new_AF_data$`Participant ID` %in% Early_cohort$Participant_ID, ]
matching_dates <- matching_data[, c("Participant ID", "Date of consenting to join UK Biobank")]
colnames(matching_dates) <- c("Participant_ID", "Date_consenting_to_join")
Early_cohort <- merge(Early_cohort, matching_dates, by.x = "Participant_ID", by.y = "Participant_ID", all.x = TRUE)

# Ensure Date_consenting_to_join is in date format
Early_cohort$Date_consenting_to_join <- as.Date(Early_cohort$Date_consenting_to_join)
# Sort by date
Early_cohort <- Early_cohort[order(Early_cohort$Date_consenting_to_join), ]

# Calculate the cutoff point, with the first 80% as training+validation set, and the last 20% as independent external validation set
cutoff_index <- floor(0.8 * nrow(Early_cohort))

# Split the dataset
train_val_data <- Early_cohort[1:cutoff_index, ]  # First 80% as training+validation
external_test_data <- Early_cohort[(cutoff_index + 1):nrow(Early_cohort), ]  # Last 20% as external test set

# 3. Extract features and target variable
predictors <- significant_early$Protein  # Extract protein features

# Ensure AF column is a factor type
train_val_data$AF <- as.factor(train_val_data$AF)
external_test_data$AF <- as.factor(external_test_data$AF)

# Load pre-trained model
lrn_lgb <- readRDS("~/Early-onset AF/Data/lrn_lgb.RDS")

# Extract feature importance and convert to data frame, including protein names
importance <- lrn_lgb$importance()
protein_importance <- data.frame(
  Protein = names(importance),  # Extract protein names
  Importance = importance       # Extract importance values
)
# Sort by importance in descending order
protein_importance <- protein_importance[order(-protein_importance$Importance), ]

# Initialize
selected_proteins <- c()
results <- data.frame(
  Step = integer(),
  Protein = character(),
  AUC = numeric(),
  AUC_CI_lower = numeric(),
  AUC_CI_upper = numeric(),
  DeLong_p_value = numeric()
)
max_auc <- 0
last_roc_obj <- NULL

# Forward stepwise selection
for (protein in protein_importance$Protein) {  #i in seq(1, nrow(protein_importance), by = 3
  # Add current protein to the feature set
  current_proteins <- c(selected_proteins, protein)
  
  # Create tasks
  task_train <- TaskClassif$new("AF_train", backend = train_val_data[, c(current_proteins, "AF")], target = "AF", positive = "1")
  task_test <- TaskClassif$new("AF_external_test", backend = external_test_data[, c(current_proteins, "AF")], target = "AF", positive = "1")
  
  # Train the model
  lrn_lgb$train(task_train)
  
  # Predict on the external test set
  prediction <- lrn_lgb$predict(task_test)
  
  # Calculate AUC
  roc_obj <- roc(external_test_data$AF, prediction$prob[, 2])  # Use second column for probability values
  auc_value <- auc(roc_obj)
  ci_auc <- ci.auc(roc_obj)  # Calculate 95% confidence interval for AUC
  
  # Calculate DeLong p-value, if there's a previous ROC object
  p_value <- NA
  if (!is.null(last_roc_obj)) {
    p_value <- roc.test(last_roc_obj, roc_obj, method = "delong")$p.value
  }
  
  # Store the current results
  results <- rbind(results, data.frame(
    Step = nrow(results) + 1,
    Protein = protein,
    AUC = auc_value,
    AUC_CI_lower = ci_auc[1],
    AUC_CI_upper = ci_auc[3],
    DeLong_p_value = p_value
  ))
  
  # Update
  selected_proteins <- current_proteins
  max_auc <- auc_value
  last_roc_obj <- roc_obj
}

# Output all results
print(results)

# Save results
saveRDS(results, file = "~/Early-onset AF/Data/external_test_selection_results.RDS")


####  ---------------------------------------------------------------------------
####   The code block below compares the performance of multiple LightGBM models 
####   for atrial fibrillation prediction on both the training and external validation 
####   datasets. The models include protein-only, clinical factors, combined protein and 
####   clinical factors, and simplified clinical factors. It computes the AUC and confidence 
####   intervals for each model, compares them pairwise using the DeLong test, 
####   and applies Bonferroni correction to adjust for multiple comparisons. The results, 
####   including p-values and significance after Bonferroni correction, are saved for both 
####   training and validation datasets.
####  ---------------------------------------------------------------------------
library(pROC)
library(dplyr)

# Load data
new_AF_data <- readRDS("~/Early-onset AF/Data/new_AF_data.RDS")
Early_cohort <- readRDS("~/Early-onset AF/Data/Early_cohort.RDS")

# Match date data
matching_data <- new_AF_data[new_AF_data$`Participant ID` %in% Early_cohort$Participant_ID, ]
matching_dates <- matching_data[, c("Participant ID", "Date of consenting to join UK Biobank")]
colnames(matching_dates) <- c("Participant_ID", "Date_consenting_to_join")
Early_cohort <- merge(Early_cohort, matching_dates, by.x = "Participant_ID", by.y = "Participant_ID", all.x = TRUE)

# Ensure Date_consenting_to_join is in date format
Early_cohort$Date_consenting_to_join <- as.Date(Early_cohort$Date_consenting_to_join)
Early_cohort <- Early_cohort[order(Early_cohort$Date_consenting_to_join), ]

# Split the dataset
cutoff_index <- floor(0.8 * nrow(Early_cohort))
train_val_data <- Early_cohort[1:cutoff_index, ]
external_test_data <- Early_cohort[(cutoff_index + 1):nrow(Early_cohort), ]

# Load models
models <- list(
  model1 = readRDS("~/Early-onset AF/Data/模型比较/纯蛋白/lrn_lgb.RDS"),
  model2 = readRDS("~/Early-onset AF/Data/模型比较/临床因素/lrn_lgb.RDS"),
  model3 = readRDS("~/Early-onset AF/Data/模型比较/蛋白和临床因素/lrn_lgb.RDS"),
  model4 = readRDS("~/Early-onset AF/Data/模型比较/蛋白和简易临床因素/lrn_lgb.RDS")
)

# Ensure AF column is factor type
train_val_data$AF <- factor(train_val_data$AF, levels = c("0", "1"))
external_test_data$AF <- factor(external_test_data$AF, levels = c("0", "1"))

# Result data frame
results <- data.frame(Model = character(),
                      Dataset = character(),
                      AUC_modelA = numeric(),
                      CI_lower_modelA = numeric(),
                      CI_upper_modelA = numeric(),
                      AUC_modelB = numeric(),
                      CI_lower_modelB = numeric(),
                      CI_upper_modelB = numeric(),
                      P_value = numeric(),
                      Bonferroni = character(),
                      stringsAsFactors = FALSE)

# Calculate AUC for training and validation datasets
roc_list_train <- lapply(models, function(model) {
  pred_train <- model$predict_newdata(train_val_data)
  roc(train_val_data$AF, pred_train$prob[, 1])
})

roc_list_test <- lapply(models, function(model) {
  pred_test <- model$predict_newdata(external_test_data)
  roc(external_test_data$AF, pred_test$prob[, 1])
})

# Calculate AUC and confidence intervals
auc_train <- sapply(roc_list_train, auc)
ci_train <- sapply(roc_list_train, ci.auc)
auc_test <- sapply(roc_list_test, auc)
ci_test <- sapply(roc_list_test, ci.auc)

# Pairwise AUC comparison and record results
for (i in 1:(length(models) - 1)) {
  for (j in (i + 1):length(models)) {
    # Training set comparison
    p_value_train <- roc.test(roc_list_train[[i]], roc_list_train[[j]], method = "delong")$p.value
    results <- rbind(results, data.frame(Model = paste("Model", i, "vs", "Model", j),
                                         Dataset = "Training",
                                         AUC_modelA = auc_train[i],
                                         CI_lower_modelA = ci_train[1, i],
                                         CI_upper_modelA = ci_train[3, i],
                                         AUC_modelB = auc_train[j],
                                         CI_lower_modelB = ci_train[1, j],
                                         CI_upper_modelB = ci_train[3, j],
                                         P_value = p_value_train,
                                         Bonferroni = NA,
                                         stringsAsFactors = FALSE))
    
    # Validation set comparison
    p_value_test <- roc.test(roc_list_test[[i]], roc_list_test[[j]], method = "delong")$p.value
    results <- rbind(results, data.frame(Model = paste("Model", i, "vs", "Model", j),
                                         Dataset = "Validation",
                                         AUC_modelA = auc_test[i],
                                         CI_lower_modelA = ci_test[1, i],
                                         CI_upper_modelA = ci_test[3, i],
                                         AUC_modelB = auc_test[j],
                                         CI_lower_modelB = ci_test[1, j],
                                         CI_upper_modelB = ci_test[3, j],
                                         P_value = p_value_test,
                                         Bonferroni = NA,
                                         stringsAsFactors = FALSE))
  }
}

# Bonferroni correction
for (i in 1:nrow(results)) {
  if (!is.na(results$P_value[i])) {
    if (results$Dataset[i] == "Training") {
      bonferroni_criteria <- results$P_value[i] < 0.05 / choose(length(models), 2)
    } else {  # Validation
      bonferroni_criteria <- results$P_value[i] < 0.05 / choose(length(models), 2)
    }
    results$Bonferroni[i] <- ifelse(bonferroni_criteria, "True", "False")
  }
}

# Print results
print(results)
saveRDS(results, file = "~/Early-onset AF/Data/多模型比较.RDS")

# Create data frames for training and validation results
train_results <- data.frame(
  Model = names(models),
  Dataset = "Training",
  AUC = auc_train,
  CI_lower = ci_train[1, ],
  CI_upper = ci_train[3, ],
  stringsAsFactors = FALSE
)

saveRDS(train_results, file = "~/Early-onset AF/Data/多模型比较训练集结果.RDS")

validation_results <- data.frame(
  Model = names(models),
  Dataset = "Validation",
  AUC = auc_test,
  CI_lower = ci_test[1, ],
  CI_upper = ci_test[3, ],
  stringsAsFactors = FALSE
)

saveRDS(validation_results, file = "~/Early-onset AF/Data/多模型比较验证集结果.RDS")

# Print results
print(train_results)
print(validation_results)

####  ---------------------------------------------------------------------------
####   The code block below performs a comprehensive single-cell RNA sequencing 
####   (scRNA-seq) analysis on human atrial samples, including quality control (QC), 
####   filtering of cells, removal of doublets, and integration of multiple datasets 
####   from different patients. It also includes visualizing cell types, identifying 
####   key genes from machine learning models, and performing differential gene expression 
####   analysis between atrial fibrillation (AF) and control samples.
####  ---------------------------------------------------------------------------
library(Seurat)
library(DoubletFinder)
library(ggplot2)
library(reshape2)

# Load and preprocess multiple datasets (AF1 to AF7, C1 to C5)
setwd("~/AF_atlas/raw_data/human_atrial_data")
counts <- Read10X(data.dir = "AF1")
Object <- CreateSeuratObject(counts = counts, project = "AF1", min.features = 200)
Object <- RenameCells(Object, add.cell.id = "AF1")
Object[["percent.mt"]] <- PercentageFeatureSet(Object, pattern = "^MT-")
VlnPlot(Object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

Object_QC <- subset(Object, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & 
                      nCount_RNA > 10^(mean(log10(Object$nCount_RNA)) - 2*sd(log10(Object$nCount_RNA))) &
                      nCount_RNA < 10^(mean(log10(Object$nCount_RNA)) + 2*sd(log10(Object$nCount_RNA))) & 
                      percent.mt < 25)
Object_QC[["condition"]] <- "AF1"
AF1 <- Object_QC

# Repeat similar steps for AF2 to AF7 and C1 to C5
# For brevity, details for these datasets are summarized in this format

# Combine all datasets
AF_combine <- merge(AF1, c(AF2, AF3, AF4, AF5, AF6, AF7, C1, C2, C3, C4, C5))
Idents(AF_combine) <- AF_combine$condition

# Define QC function for removing mitochondrial (MT) and ERCC genes and filtering doublets
single_cell_qc <- function(patient_id, combined_data) {
  patient_data <- subset(combined_data, idents = patient_id)
  mt_ercc_genes <- grep("^(MT-|ERCC)", rownames(patient_data), value = TRUE)
  patient_data_filtered <- patient_data[!rownames(patient_data) %in% mt_ercc_genes, ]
  patient_data_SCT <- SCTransform(patient_data_filtered)
  patient_data_SCT <- RunPCA(patient_data_SCT)
  patient_data_SCT <- RunUMAP(patient_data_SCT, dims = 1:10)
  
  # Remove doublets using DoubletFinder
  sweep.res.list <- paramSweep_v3(patient_data_SCT, PCs = 1:10, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  nExp_poi <- round(0.075 * nrow(patient_data_SCT@meta.data))
  
  patient_data_filtered <- doubletFinder_v3(patient_data_SCT, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, sct = TRUE)
  df_classification_colname <- grep("^DF.classifications", colnames(patient_data_filtered@meta.data), value = TRUE)
  patient_data_filtered <- patient_data_filtered[, patient_data_filtered@meta.data[, df_classification_colname] == "Singlet"]
  
  return(patient_data_filtered)
}

# Apply QC function to all datasets
AF_QC <- NULL
for (patient_id in unique(Idents(AF_combine))) {
  filtered_data <- single_cell_qc(patient_id, AF_combine)
  AF_QC <- if (is.null(AF_QC)) filtered_data else merge(AF_QC, filtered_data)
}

# Save the filtered dataset
saveRDS(AF_QC, file = "~/AF_atlas/raw_data/human_atrial_data/AF_QC.RDS")

# Split dataset by condition and merge patients with low cell count
Object.list <- SplitObject(AF_QC, split.by = "condition")
Object.list <- lapply(Object.list, function(patient_data) {
  if (ncol(patient_data) < 500) merge_small_patients(patient_data, Object.list)
  return(patient_data)
})

# Perform SCTransform, integration, and dimensionality reduction
for (patient_id in names(Object.list)) {
  Object.list[[patient_id]] <- SCTransform(Object.list[[patient_id]], verbose = TRUE)
}

features <- SelectIntegrationFeatures(object.list = Object.list, nfeatures = 3000)
Object.anchors <- FindIntegrationAnchors(object.list = Object.list, normalization.method = 'SCT', anchor.features = features)
Object.combined.sct <- IntegrateData(anchorset = Object.anchors, normalization.method = 'SCT')

# Save integrated data
saveRDS(Object.combined.sct, file = "~/AF_atlas/Data/single_cell_analysis/Object.combined.sct.RDS")

# Run PCA, UMAP, clustering, and marker analysis
Object.combined.sct <- RunPCA(Object.combined.sct)
Object.combined.sct <- RunUMAP(Object.combined.sct, dims = 1:20)
Object.combined.sct <- FindNeighbors(Object.combined.sct, dims = 1:20)
Object.combined.sct <- FindClusters(Object.combined.sct, resolution = 0.5)

# Rename clusters to cell types
Object.combined.sct <- RenameIdents(Object.combined.sct,
                                    "11" = "EC", "3" = "DC", "14" = "SMC", 
                                    "8" = "FB", "2" = "T_cell", "9" = "MP", 
                                    "7" = "Neutrophil")

# Save final data
saveRDS(Object.combined.sct, file = "~/Early-onset AF/Data/AF_renamed.RDS")

# Extract top proteins from a LightGBM model and visualize
lrn_lgb <- readRDS("~/Early-onset AF/Data/model_comparison/protein_only/lrn_lgb.RDS")
importance_data <- lrn_lgb$importance()
top_9_proteins <- names(sort(importance_data, decreasing = TRUE))[1:9]

# Dot plot and violin plot for top proteins across conditions
DotPlot(Object.combined.sct, features = top_9_proteins)
VlnPlot(Object.combined.sct, features = top_9_proteins, stack = TRUE, split.by = "condition", assay = "SCT", flip = TRUE) +
  scale_fill_manual(values = c("AF" = "#B22222", "Control" = "#4682B4"))

# Differential gene expression in FB and EC cells
FB_cells <- subset(Object.combined.sct, idents = "FB")
FB_diff_genes <- FindMarkers(FB_cells, ident.1 = "AF", ident.2 = "Control", group.by = "condition", assay = "SCT")
print(FB_diff_genes["IGFBP7", ])

EC_cells <- subset(Object.combined.sct, idents = "EC")
EC_diff_genes <- FindMarkers(EC_cells, ident.1 = "AF", ident.2 = "Control", group.by = "condition", assay = "SCT")
print(EC_diff_genes["IGFBP7", ])





















