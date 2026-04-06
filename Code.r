################################################################################
# Project: Broad-scale patterns of large herbivore species richness are linked to movement-related traits
# Description: Integrated workflow: VIF selection -> Spatial COM-Poisson (spaMM) 
#              -> IEP Calculation -> Visualization.
# Author: Xueting Zheng / Nanjing University
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
# Set working directory to project root (assumes .Rproj location)
project_root <- getwd() 
data_path    <- file.path(project_root, "Processed_Data", "VIF_Input")
output_root  <- file.path(project_root, "Results")

# Create output directory if it doesn't exist
if (!dir.exists(output_root)) dir.create(output_root, recursive = TRUE)

# Model Configuration Mapping
# Symbols: A=local, D=context, C=macro_fixed, E=move_trait, F_var=bodymass, MS=macro_scale
model_configs <- list(
  "lcmxfMS" = c("A", "D", "C", "E", "F_var", "MS"),
  "lcmfMS"  = c("A", "D", "C", "F_var", "MS"),
  "lmfMS"   = c("A", "C", "F_var", "MS"),
  "lmxfMS"  = c("A", "C", "E", "F_var", "MS")
)

# Define spatial scales analyzed in the study
scales <- c("10", "20", "40", "100", "200", "400") 

# ==============================================================================
# STEP 1: VIF (Variance Inflation Factor) Selection
# ==============================================================================
PROTECTED_VARS <- c("NPP", "PSN", "TSN") # Do not remove these key predictors

calculate_vif <- function(df) {
  df_scaled <- as.data.frame(scale(df))
  model <- lm(SR ~ ., data = df_scaled)
  vif_values <- vif(model)
  return(data.frame(col = names(vif_values), vif = as.numeric(vif_values)))
}

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

# Run Step 1: Automated VIF cleaning
message("Starting Step 1: VIF Selection...")
flist <- list.files(data_path, pattern = ".csv", full.names = TRUE)
lapply(flist, function(name) {
  df <- read.csv(name)
  id_col <- df$Id
  df_vars <- df[, !names(df) %in% c("Id")]
  protected <- find_protected_vars(colnames(df_vars))
  df_others <- df_vars[, !colnames(df_vars) %in% protected, drop = FALSE]
  df_prot   <- df_vars[, colnames(df_vars) %in% protected, drop = FALSE]
  
  if (ncol(df_others) > 0) {
    repeat {
      vifs <- calculate_vif(cbind(df_prot, df_others))
      non_prot_vifs <- vifs[!vifs$col %in% protected, ]
      if (nrow(non_prot_vifs) == 0 || max(non_prot_vifs$vif) <= 3) break
      var_to_remove <- non_prot_vifs$col[which.max(non_prot_vifs$vif)]
      df_others <- df_others[, colnames(df_others) != var_to_remove, drop = FALSE]
    }
  }
  write.csv(cbind(Id = id_col, df_prot, df_others), 
            file.path(data_path, paste0("vif_cleaned_", basename(name))), row.names = FALSE)
})

# ==============================================================================
# STEP 2: Spatial COM-Poisson Modeling 
# ==============================================================================

# Core function for single iteration
fit_single_rep <- function(rep_idx, data_source, coords_source, id_source, formula_obj, sample_size, dmax_val) {
  library(spaMM)
  library(spdep)
  set.seed(1000 + rep_idx)
  indices <- sample(seq_len(nrow(data_source)), sample_size)
  s_data <- data_source[indices, ]
  s_coords <- coords_source[indices, ]
  
  nb <- dnearneigh(s_coords, 0, dmax_val)
  W <- nb2mat(nb, style = "B", zero.policy = TRUE)
  rownames(W) <- colnames(W) <- id_source[indices]
  
  tryCatch({
    fit <- corrHLfit(formula = formula_obj, data = s_data, family = COMPoisson(),
                     adjMatrix = W, method = "ML", control.HLfit = list(LevenbergM = TRUE))
    nu_val <- environment(fit$family$linkfun)$nu
    null_glm <- glm(formula_obj, data = s_data, family = COMPoisson(nu = nu_val))
    pseudo_r2 <- 1 - (null_glm$deviance / null_glm$null.deviance)
    res_row <- as.data.frame(t(c(fixef(fit), PseudoR2 = pseudo_r2)))
    return(list(combined_results = res_row, success = TRUE))
  }, error = function(e) list(success = FALSE, error = e$message))
}

n_reps <- 99           # Set to 2 or 3 for a quick test run before full analysis 
n_cores <- parallel::detectCores() - 2  
sample_size <- 2772    # Sample size for each iteration
dmax_val <- 200000     # Max distance for spatial weight matrix (adjust units as needed)           
           

# Formula builder
get_formula <- function(config) {
  # !!! IMPORTANT: Replace these strings with your actual column names from CSV !!!
  vars <- c("Local_Env", "Macro_Climate", "Body_Mass") 
  if ("D" %in% config) vars <- c(vars, "Landscape_Context")
  if ("E" %in% config) vars <- c(vars, "Movement_Trait")
  
  as.formula(paste("SR ~", paste(vars, collapse = " + "), "+ (1|adjMatrix)"))
}

for (model_name in names(model_configs)) {
  formula_obj <- get_formula(model_configs[[model_name]])
  for (scale in scales) {
    message(sprintf("Fitting: %s at %s km...", model_name, scale))
    data_file <- file.path(project_root, "Processed_Data", paste0("vif_cleaned_data_", scale, "km.csv"))
    if (!file.exists(data_file)) next
    
    current_data <- read.csv(data_file)
    cl <- makeCluster(n_cores)
    clusterExport(cl, c("fit_single_rep", "current_data", "formula_obj", "sample_size", "dmax_val"))
    clusterEvalQ(cl, { library(spaMM); library(spdep) })
    
    results <- parLapply(cl, 1:n_reps, function(i) {
      fit_single_rep(i, current_data, current_data[,c("X","Y")], current_data$Id, formula_obj, sample_size, dmax_val)
    })
    stopCluster(cl)
    
    success_results <- do.call(rbind, lapply(results, function(x) if(x$success) x$combined_results))
    if (!is.null(success_results)) {
      success_results$Rep <- 1:nrow(success_results)
      out_path <- file.path(output_root, paste0(model_name, "_", scale, "km"))
      if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)
      write.csv(success_results, file.path(out_path, "R2_results.csv"), row.names = FALSE)
    }
  }
}

# ==============================================================================
# STEP 3: IEP Calculation 
# ==============================================================================
comparison_pairs <- list(
  "Landscape_Context_IEP" = c(full = "lcmxfMS", reduced = "lmxfMS"),
  "Movement_Traits_IEP"   = c(full = "lcmxfMS", reduced = "lcmfMS")
)

for (pair_name in names(comparison_pairs)) {
  pair <- comparison_pairs[[pair_name]]
  results_list <- list()
  for (scale in scales) {
    full_p <- file.path(output_root, paste0(pair["full"], "_", scale, "km"), "R2_results.csv")
    red_p  <- file.path(output_root, paste0(pair["reduced"], "_", scale, "km"), "R2_results.csv")
    if (file.exists(full_p) && file.exists(red_p)) {
      merged <- merge(read.csv(full_p)[,c("Rep","PseudoR2")], read.csv(red_p)[,c("Rep","PseudoR2")], by="Rep", suffixes=c("_f","_r"))
      merged[[paste0("Scale_", scale)]] <- (merged$PseudoR2_f - merged$PseudoR2_r) * 100
      results_list[[scale]] <- merged[, c("Rep", paste0("Scale_", scale))]
    }
  }
  if (length(results_list) > 0) {
    final_iep <- Reduce(function(x, y) merge(x, y, by = "Rep", all = TRUE), results_list)
    write.csv(final_iep, file.path(output_root, paste0(pair_name, "_Summary.csv")), row.names = FALSE)
  }
}

# ==============================================================================
# STEP 4: Visualization with Robust Data Cleaning (MAD)
# ==============================================================================

# Robust cleaning function using MAD
clean_iep_data <- function(vec, mad_threshold = 3) {
  # 1. Remove NA and basic illogical values (IEP shouldn't be negative or physically impossible)
  vec <- vec[!is.na(vec) & vec >= 0 & vec <= 100]
  
  if (length(vec) < 3) return(numeric(0))
  
  # 2. MAD Calculation
  median_val <- median(vec)
  mad_val <- mad(vec, constant = 1.4826) # standard constant for normal distribution consistency
  
  # 3. Define bounds
  lower_bound <- median_val - mad_threshold * mad_val
  upper_bound <- median_val + mad_threshold * mad_val
  
  # 4. Filter
  cleaned_vec <- vec[vec >= lower_bound & vec <= upper_bound]
  return(cleaned_vec)
}

plot_iep <- function(input_csv, out_dir) {
  data <- read.csv(input_csv)
  
  # Prepare a list to collect cleaned data per scale
  cleaned_list <- list()
  
  # Identify scale columns (assuming names like "Scale_10", "Scale_20"...)
  scale_cols <- grep("Scale_", colnames(data), value = TRUE)
  
  for (col in scale_cols) {
    raw_values <- data[[col]]
    cleaned_values <- clean_iep_data(raw_values)
    
    # 5. Quality Control: Only keep scales with at least 3 valid samples
    if (length(cleaned_values) >= 3) {
      cleaned_list[[col]] <- data.frame(
        Rep = 1:length(cleaned_values),
        Scale = col,
        IEP = cleaned_values
      )
    } else {
      message(sprintf("Warning: Scale [%s] in [%s] skipped due to insufficient valid samples.", 
                      col, basename(input_csv)))
    }
  }
  
  if (length(cleaned_list) == 0) {
    message(sprintf("Skipping plot for [%s]: No valid data remaining after cleaning.", basename(input_csv)))
    return(NULL)
  }
  
  # Combine cleaned data back to long format
  df_long <- do.call(rbind, cleaned_list)
  df_long$Scale <- factor(gsub("Scale_", "", df_long$Scale), levels = scales)
  
  # 6. Plotting
  p <- ggplot(df_long, aes(x = Scale, y = IEP)) +
    stat_boxplot(geom = 'errorbar', width = 0.25) +
    geom_boxplot(fill = "white", 
                 notch = TRUE, 
                 outlier.shape = 16, # Show outliers as dots after MAD filtering
                 outlier.size = 1,
                 linewidth = 0.8) +
    theme_classic() +
    labs(title = gsub("_Summary.csv", "", basename(input_csv)), 
         x = "Spatial Scale (km)", 
         y = "Independent Explanatory Power (%)") +
    coord_cartesian(ylim = c(0, 25)) # Focus on realistic range
  
  ggsave(file.path(out_dir, paste0(basename(input_csv), "_cleaned.png")), p, width = 7, height = 6)
}

# Execute Step 4 automatically
iep_summaries <- list.files(output_root, pattern = "_Summary.csv", full.names = TRUE)
lapply(iep_summaries, function(f) plot_iep(f, output_root))

# Execute Step 4 automatically
iep_summaries <- list.files(output_root, pattern = "_Summary.csv", full.names = TRUE)
lapply(iep_summaries, function(f) plot_iep(f, output_root))

message("Full Workflow Complete. Results and Plots are in: ", output_root)
