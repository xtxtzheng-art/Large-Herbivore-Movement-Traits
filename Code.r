################################################################################
# Project: Broad-scale patterns of large herbivore species richness are linked to movement-related traits
# Description: Multi-step workflow including VIF selection, Spatial COM-Poisson 
#              modeling (spaMM), Pseudo-R2 calculation, and Visualization.
# Author: Xueting Zheng/Nanjing University
# Year: 2026
################################################################################

# Load Required Libraries
library(car)
library(MASS)
library(spaMM)
library(spdep)
library(dplyr)
library(parallel)
library(ggplot2)
library(reshape2)

# ==============================================================================
# 0. Global Configurations
# ==============================================================================
# Change this path to your local repository root
project_root <- getwd() 
data_path    <- file.path(project_root, "Processed_Data", "VIF_Input")
output_root  <- file.path(project_root, "Results")

# Model Configuration Mapping
# A=local, D=context, C=macro_fixed, E=move_trait, F_var=bodymass, MS=macro_scale
model_configs <- list(
  "lcmxfMS" = c("A", "D", "C", "E", "F_var", "MS"),
  "lcmfMS"  = c("A", "D", "C", "F_var", "MS"),
  "lmfMS"   = c("A", "C", "F_var", "MS"),
  "lmxfMS"  = c("A", "C", "E", "F_var", "MS")
)

# ==============================================================================
# STEP 1: VIF (Variance Inflation Factor) Selection
# Purpose: Remove multicollinearity while protecting key environmental variables.
# ==============================================================================

PROTECTED_VARS <- c("NPP", "PSN", "TSN")

# Helper: Calculate VIF for a dataframe
calculate_vif <- function(df) {
  df_scaled <- as.data.frame(scale(df))
  model <- lm(SR ~ ., data = df_scaled)
  vif_values <- vif(model)
  return(data.frame(col = names(vif_values), vif = as.numeric(vif_values)))
}

# Helper: Identify protected variables using pattern matching
find_protected_vars <- function(variable_names) {
  found <- c()
  for (var in variable_names) {
    var_upper <- toupper(var)
    for (p in PROTECTED_VARS) {
      pattern <- paste0("^", p, "$|^", p, "[._]|[._]", p, "$|_", p, "_|", p, "\\.")
      if (grepl(pattern, var_upper)) {
        found <- c(found, var)
        break
      }
    }
  }
  return(unique(found))
}

# Run Step 1 Processing
flist <- list.files(data_path, pattern = ".csv", full.names = TRUE)
lapply(flist, function(name) {
  df <- read.csv(name)
  id_col <- df$Id
  df <- df[, !names(df) %in% c("Id")]
  
  protected <- find_protected_vars(colnames(df))
  
  # Recursive VIF removal for non-protected variables
  df_others <- df[, !colnames(df) %in% protected, drop = FALSE]
  df_prot   <- df[, colnames(df) %in% protected, drop = FALSE]
  
  if (ncol(df_others) > 0) {
    repeat {
      vifs <- calculate_vif(cbind(df_prot, df_others))
      non_prot_vifs <- vifs[!vifs$col %in% protected, ]
      
      if (nrow(non_prot_vifs) == 0 || max(non_prot_vifs$vif) <= 3) break
      
      var_to_remove <- non_prot_vifs$col[which.max(non_prot_vifs$vif)]
      df_others <- df_others[, colnames(df_others) != var_to_remove, drop = FALSE]
    }
  }
  
  df_final <- cbind(Id = id_col, df_prot, df_others)
  write.csv(df_final, file.path(data_path, paste0("vif_cleaned_", basename(name))), row.names = FALSE)
})

# ==============================================================================
# STEP 2: Spatial COM-Poisson Model Fitting & Pseudo-R2
# Purpose: Estimate coefficients using spatial correlation structures.
# ==============================================================================

# Internal function for parallel execution
fit_single_rep <- function(rep_idx, data_source, coords_source, id_source, formula_obj, sample_size, dmax_val) {
  library(spaMM)
  library(spdep)
  set.seed(1000 + rep_idx)
  
  indices <- sample(seq_len(nrow(data_source)), sample_size)
  s_data <- data_source[indices, ]
  s_coords <- coords_source[indices, ]
  
  # Create Spatial Weight Matrix
  nb <- dnearneigh(s_coords, 0, dmax_val)
  W <- nb2mat(nb, style = "B", zero.policy = TRUE)
  rownames(W) <- colnames(W) <- id_source[indices]
  
  tryCatch({
    # Fit Spatial COM-Poisson
    fit <- corrHLfit(formula = formula_obj, data = s_data, family = COMPoisson(),
                     adjMatrix = W, method = "ML", control.HLfit = list(LevenbergM = TRUE))
    
    # Calculate Pseudo-R2 (1 - Deviance/Null Deviance)
    nu_val <- environment(fit$family$linkfun)$nu
    null_glm <- glm(formula_obj, data = s_data, family = COMPoisson(nu = nu_val))
    pseudo_r2 <- 1 - (null_glm$deviance / null_glm$null.deviance)
    
    # Export Coefficients and Stats
    res_row <- as.data.frame(t(c(fixef(fit), PseudoR2 = pseudo_r2, LogLik = as.numeric(logLik(fit)))))
    return(list(combined_results = res_row, rep_idx = rep_idx, success = TRUE))
  }, error = function(e) list(success = FALSE, error = e$message))
}

# --- Parallel Execution Configuration ---
n_reps <- 100          # Number of iterations (bootstrap/sampling)
# Automatically detect CPU cores and reserve 2 for system stability
n_cores <- parallel::detectCores() - 2  
sample_size <- 2772    # Sample size for each iteration
dmax_val <- 200000     # Max distance for spatial weight matrix (adjust units as needed)

# Helper function to dynamically generate formulas based on model configurations
# Maps symbols (A, D, C, E) to actual variable groups in your dataset
get_formula <- function(config) {
  # Base variables (Local, Macro, BodyMass)
  # Ensure these variable names match your column headers in Processed_Data
  base_vars <- c("Local_Var", "Macro_Var", "BodyMass") 
  
  if ("D" %in% config) base_vars <- c(base_vars, "Context_Var")
  if ("E" %in% config) base_vars <- c(base_vars, "Trait_Var")
  
  # Construct spaMM formula with spatial term (adjMatrix)
  formula_str <- paste("SR ~", paste(base_vars, collapse = " + "), "+ (1|adjMatrix)")
  return(as.formula(formula_str))
}

# --- Main Modeling Loop ---
# Iterates through all model configurations (lcmxfMS, etc.) and spatial scales
for (model_name in names(model_configs)) {
  config <- model_configs[[model_name]]
  formula_obj <- get_formula(config)
  
  for (scale in scales) {
    message(sprintf("Processing: Model [%s] at Scale [%s km]...", model_name, scale))
    
    # 1. Load data for the specific scale
    # Expected path: project_root/Processed_Data/scale_data_10km.csv
    data_file <- file.path(project_root, "Processed_Data", paste0("scale_data_", scale, "km.csv"))
    
    if (!file.exists(data_file)) {
      warning(paste("Data file not found, skipping:", data_file))
      next
    }
    
    current_data <- read.csv(data_file)
    current_coords <- current_data[, c("X", "Y")] # Ensure X and Y coordinates exist
    
    # 2. Setup Parallel Cluster
    cl <- makeCluster(n_cores)
    
    # Export necessary objects and libraries to each worker node
    clusterExport(cl, c("fit_single_rep", "current_data", "current_coords", 
                        "formula_obj", "sample_size", "dmax_val"))
    clusterEvalQ(cl, {
      library(spaMM)
      library(spdep)
    })
    
    # 3. Execute Parallel Computing via parLapply
    results <- parLapply(cl, 1:n_reps, function(i) {
      fit_single_rep(i, current_data, current_coords, current_data$Id, 
                     formula_obj, sample_size, dmax_val)
    })
    
    stopCluster(cl) # Free up system resources
    
    # 4. Extract and Save Successful Results
    success_results <- do.call(rbind, lapply(results, function(x) if(x$success) x$combined_results))
    
    if (!is.null(success_results)) {
      success_results$Rep <- 1:nrow(success_results)
      
      # Create output directory: output_root/model_scale/
      out_path <- file.path(output_root, paste0(model_name, "_", scale, "km"))
      if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)
      
      write.csv(success_results, file.path(out_path, "R2_results.csv"), row.names = FALSE)
      message(paste("Successfully saved modeling results to:", out_path))
    }
  }
}

# ==============================================================================
# STEP 3: Independent Explanatory Power (IEP) Calculation
# Purpose: Calculate Delta Pseudo-R2 between full and reduced models.
# ==============================================================================

# This section iterates through the output folders of Step 2,
# aligns "Rep" IDs, and subtracts R2 values (e.g., Full - Reduced).
# Example: Calculating IEP of "Landscape Context" 
# by comparing lcmxfMS (Full) and lmxfMS (Reduced without context)
comparison_pairs <- list(
  "Landscape_Context_IEP" = c(full = "lcmxfMS", reduced = "lmxfMS"),
  "Movement_Traits_IEP"   = c(full = "lcmxfMS", reduced = "lcmfMS")
)

scales <- c("10", "20", "40", "100", "200", "400") # Define scales analyzed

for (pair_name in names(comparison_pairs)) {
  pair <- comparison_pairs[[pair_name]]
  results_list <- list()
  
  for (scale in scales) {
    # Construct folder paths based on Step 2 output structure
    full_path <- file.path(output_root, paste0(pair["full"], "_", scale, "km"), "R2_results.csv")
    red_path  <- file.path(output_root, paste0(pair["reduced"], "_", scale, "km"), "R2_results.csv")
    
    if (file.exists(full_path) && file.exists(red_path)) {
      full_df <- read.csv(full_path)
      red_df  <- read.csv(red_path)
      
      # Ensure reps are aligned by ID (assuming 'Rep' column exists)
      merged <- merge(full_df[, c("Rep", "PseudoR2")], 
                      red_df[, c("Rep", "PseudoR2")], 
                      by = "Rep", suffixes = c("_full", "_red"))
      
      # Calculate IEP: (Full R2 - Reduced R2) * 100 for percentage
      merged[[paste0("Scale_", scale)]] <- (merged$PseudoR2_full - merged$PseudoR2_red) * 100
      
      results_list[[scale]] <- merged[, c("Rep", paste0("Scale_", scale))]
    } else {
      warning(paste("Missing data for", pair_name, "at scale", scale))
    }
  }
  
  # Combine all scales into one table
  if (length(results_list) > 0) {
    final_iep_df <- Reduce(function(x, y) merge(x, y, by = "Rep", all = TRUE), results_list)
    
    # Save the IEP results for Visualization (Step 4)
    iep_output_path <- file.path(output_root, paste0(pair_name, "_Summary.csv"))
    write.csv(final_iep_df, iep_output_path, row.names = FALSE)
    message(paste("IEP calculation completed for:", pair_name))
    
    # Optional: Automatically trigger Step 4 plot
    # plot_iep_distribution(iep_output_path, output_root)
  }
}

# ==============================================================================
# STEP 4: Data Visualization (Boxplots)
# Purpose: Visualizing the distribution of IEP across different scales.
# ==============================================================================

# Theme and Color Mapping
color_map <- c(
  "all" = "#595959", "dist_homerange" = "#8B82A2", "homerange" = "#762A83", 
  "maxspeed" = "#E8B32C", "dist" = "#4AADAC"
)

plot_iep_distribution <- function(input_csv, out_dir) {
  data <- read.csv(input_csv)
  
  # Reshape for ggplot
  df_long <- melt(data, id.vars = "Rep", variable.name = "Scale", value.name = "IEP")
  df_long$Scale <- factor(df_long$Scale, labels = c("10", "20", "40", "100", "200", "400"))
  
  p <- ggplot(df_long, aes(x = Scale, y = IEP)) +
    stat_boxplot(geom ='errorbar', width = 0.25) +
    geom_boxplot(fill = NA, notch = TRUE, outlier.shape = NA, linewidth = 1) +
    theme_classic() +
    labs(title = basename(input_csv), x = "Scale (km)", y = "Independent Explanatory Power (%)") +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
    coord_cartesian(ylim = c(0, 25))
  
  ggsave(file.path(out_dir, paste0(basename(input_csv), "_plot.png")), p, width = 8, height = 8)
}

# Final message

message("Script setup complete. Ensure paths are correctly configured before running full loops.")

