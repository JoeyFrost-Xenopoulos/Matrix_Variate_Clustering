# test_extensive_k_selection.R
# Comprehensive testing of automatic k selection across many scenarios

source("../R/Matrix_Init.R")
source("../R/Matrix.R")
source("../R/Matrix_KS_Score.R")
source("../R/Matrix_Noise.R")
source("../R/Matrix_Noise_BR.R")
source("../R/Matrix_Utils.R")

# Enhanced data generator with more control
generate_complex_data <- function(
    n_per_group, rows, cols, 
    group_means, row_covs, col_covs,
    noise_prop = 0,
    noise_type = "uniform",
    overlap = 0,  # Cluster overlap (0-1)
    outliers = 0  # Proportion of extreme outliers
) {
  g <- length(group_means)
  total_clean <- sum(n_per_group)
  n_noise <- floor(total_clean * noise_prop)
  n_outliers <- floor(total_clean * outliers)
  total_n <- total_clean + n_noise + n_outliers
  
  x_list <- vector("list", total_n)
  true_labels <- integer(total_n)
  
  # Generate clean clusters
  idx <- 1
  for (group in 1:g) {
    n_group <- n_per_group[group]
    for (i in 1:n_group) {
      # Add overlap by mixing means
      if (overlap > 0 && group < g) {
        mix_prob <- runif(1)
        if (mix_prob < overlap) {
          # Blend with next group
          actual_mean <- (group_means[[group]] + group_means[[group + 1]]) / 2
        } else {
          actual_mean <- group_means[[group]]
        }
      } else {
        actual_mean <- group_means[[group]]
      }
      
      Z <- matrix(rnorm(rows * cols), rows, cols)
      eigen_row <- eigen(row_covs[[group]])
      eigen_col <- eigen(col_covs[[group]])
      
      A <- eigen_row$vectors %*% diag(sqrt(eigen_row$values)) %*% t(eigen_row$vectors)
      B <- eigen_col$vectors %*% diag(sqrt(eigen_col$values)) %*% t(eigen_col$vectors)
      
      x_list[[idx]] <- actual_mean + A %*% Z %*% B
      true_labels[idx] <- group
      idx <- idx + 1
    }
  }
  
  # Add standard noise
  if (n_noise > 0) {
    all_data <- do.call(rbind, lapply(x_list[1:(total_clean - n_outliers)], as.vector))
    data_range <- range(all_data, na.rm = TRUE)
    data_sd <- sd(all_data, na.rm = TRUE)
    
    for (i in 1:n_noise) {
      if (noise_type == "uniform") {
        noise_mat <- matrix(runif(rows * cols, data_range[1], data_range[2]), rows, cols)
      } else if (noise_type == "gaussian") {
        noise_mat <- matrix(rnorm(rows * cols, mean(all_data), data_sd), rows, cols)
      } else if (noise_type == "cauchy") {
        noise_mat <- matrix(rcauchy(rows * cols, 0, data_sd), rows, cols)
      }
      x_list[[idx]] <- noise_mat
      true_labels[idx] <- 0
      idx <- idx + 1
    }
  }
  
  # Add extreme outliers
  if (n_outliers > 0) {
    all_data <- do.call(rbind, lapply(x_list[1:idx-1], as.vector))
    data_range <- range(all_data, na.rm = TRUE)
    extreme_range <- data_range + c(-10, 10) * sd(all_data)
    
    for (i in 1:n_outliers) {
      x_list[[idx]] <- matrix(runif(rows * cols, extreme_range[1], extreme_range[2]), rows, cols)
      true_labels[idx] <- -1  # Special label for outliers
      idx <- idx + 1
    }
  }
  
  list(x_list = x_list, true_labels = true_labels)
}

# Comprehensive evaluation function
evaluate_comprehensive <- function(result, true_labels, x_list, test_name) {
  # Basic metrics
  n_points <- length(true_labels)
  noise_true <- true_labels == 0
  noise_pred <- result$cluster == 0
  
  # Noise detection metrics
  noise_accuracy <- mean(noise_true == noise_pred)
  noise_precision <- if(sum(noise_pred) > 0) sum(noise_true & noise_pred) / sum(noise_pred) else NA
  noise_recall <- if(sum(noise_true) > 0) sum(noise_true & noise_pred) / sum(noise_true) else NA
  noise_f1 <- if(!is.na(noise_precision) && !is.na(noise_recall) && (noise_precision + noise_recall) > 0) {
    2 * noise_precision * noise_recall / (noise_precision + noise_recall)
  } else { NA }
  
  # Cluster purity for non-noise points
  non_noise_idx <- which(true_labels > 0 & result$cluster > 0)
  if (length(non_noise_idx) > 0) {
    true_non_noise <- true_labels[non_noise_idx]
    pred_non_noise <- result$cluster[non_noise_idx]
    
    # Adjusted Rand Index approximation
    n <- length(non_noise_idx)
    agreements <- 0
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        same_true <- (true_non_noise[i] == true_non_noise[j])
        same_pred <- (pred_non_noise[i] == pred_non_noise[j])
        if (same_true == same_pred) agreements <- agreements + 1
      }
    }
    total_pairs <- n * (n - 1) / 2
    ari <- if(total_pairs > 0) agreements / total_pairs else 1
    
    # Cluster completeness
    unique_true <- unique(true_non_noise)
    unique_pred <- unique(pred_non_noise)
    cluster_mapping <- table(true_non_noise, pred_non_noise)
    cluster_purity <- sum(apply(cluster_mapping, 1, max)) / n
  } else {
    ari <- NA
    cluster_purity <- NA
  }
  
  # Timing and convergence
  iterations <- result$iterations
  converged <- result$converged
  final_loglik <- tail(result$logLik, 1)
  
  # K selection info
  selected_k <- if(!is.null(result$k_selection)) result$k_selection$selected_k else NA
  ks_stat <- if(!is.null(result$k_selection)) min(result$k_selection$ks_scores, na.rm = TRUE) else NA
  ks_pval <- if(!is.null(result$k_selection)) max(result$k_selection$ks_pvalues, na.rm = TRUE) else NA
  
  return(data.frame(
    test_name = test_name,
    selected_k = selected_k,
    noise_proportion_estimated = result$noise$pi,
    noise_proportion_true = mean(noise_true),
    noise_accuracy = noise_accuracy,
    noise_precision = noise_precision,
    noise_recall = noise_recall,
    noise_f1 = noise_f1,
    ari = ari,
    cluster_purity = cluster_purity,
    iterations = iterations,
    converged = converged,
    final_loglik = final_loglik,
    ks_statistic = ks_stat,
    ks_pvalue = ks_pval,
    stringsAsFactors = FALSE
  ))
}

# Run extensive tests
cat("\n" + paste(rep("=", 80), collapse = "") + "\n")
cat("EXTENSIVE TEST SUITE 1: Comprehensive K Selection\n")
cat(paste(rep("=", 80), collapse = "") + "\n\n")

results_all <- list()

# Test 1: Varying noise proportions
cat("Test 1: Varying Noise Proportions (0% to 30%)\n")
cat(paste(rep("-", 60), collapse = "") + "\n")

noise_levels <- c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30)
rows <- 4
cols <- 3
g <- 2

group_means <- list(matrix(0, rows, cols), matrix(3, rows, cols))
row_covs <- list(diag(rows), diag(rows))
col_covs <- list(diag(cols), diag(cols))

for (noise_level in noise_levels) {
  cat(sprintf("  Noise proportion: %.0f%%", noise_level * 100))
  
  test_data <- generate_complex_data(
    n_per_group = c(50, 50),
    rows = rows, cols = cols,
    group_means = group_means,
    row_covs = row_covs,
    col_covs = col_covs,
    noise_prop = noise_level,
    noise_type = "uniform"
  )
  
  result <- tryCatch({
    matrix_variate_noise_fit(
      x_list = test_data$x_list,
      g = g,
      noise_type = "hc",
      estimate_k = TRUE,
      verbose = FALSE,
      nstart = 10,
      max_iter = 100
    )
  }, error = function(e) {
    cat(" ✗ Error:", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(result)) {
    eval <- evaluate_comprehensive(result, test_data$true_labels, test_data$x_list, 
                                   paste0("noise_", noise_level*100, "pct"))
    results_all[[length(results_all)+1]] <- eval
    cat(sprintf(" ✓ Noise detected: %.1f%% (true: %.0f%%), ARI: %.3f\n", 
                eval$noise_proportion_estimated * 100,
                noise_level * 100,
                eval$ari))
  }
}

# Test 2: Different noise distributions
cat("\nTest 2: Different Noise Distributions\n")
cat(paste(rep("-", 60), collapse = "") + "\n")

noise_distributions <- c("uniform", "gaussian", "cauchy")
for (dist in noise_distributions) {
  cat(sprintf("  Noise distribution: %s", dist))
  
  test_data <- generate_complex_data(
    n_per_group = c(50, 50),
    rows = rows, cols = cols,
    group_means = group_means,
    row_covs = row_covs,
    col_covs = col_covs,
    noise_prop = 0.15,
    noise_type = dist
  )
  
  result <- tryCatch({
    matrix_variate_noise_fit(
      x_list = test_data$x_list,
      g = g,
      noise_type = "hc",
      estimate_k = TRUE,
      verbose = FALSE,
      nstart = 10,
      max_iter = 100
    )
  }, error = function(e) {
    cat(" ✗ Error:", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(result)) {
    eval <- evaluate_comprehensive(result, test_data$true_labels, test_data$x_list, 
                                   paste0("dist_", dist))
    results_all[[length(results_all)+1]] <- eval
    cat(sprintf(" ✓ Noise F1: %.3f, ARI: %.3f\n", eval$noise_f1, eval$ari))
  }
}

# Test 3: Varying cluster separation
cat("\nTest 3: Varying Cluster Separation\n")
cat(paste(rep("-", 60), collapse = "") + "\n")

separations <- c(0.5, 1, 2, 3, 4, 5)
for (sep in separations) {
  cat(sprintf("  Separation between clusters: %.1f", sep))
  
  group_means_sep <- list(matrix(0, rows, cols), matrix(sep, rows, cols))
  
  test_data <- generate_complex_data(
    n_per_group = c(50, 50),
    rows = rows, cols = cols,
    group_means = group_means_sep,
    row_covs = row_covs,
    col_covs = col_covs,
    noise_prop = 0.10
  )
  
  result <- tryCatch({
    matrix_variate_noise_fit(
      x_list = test_data$x_list,
      g = g,
      noise_type = "hc",
      estimate_k = TRUE,
      verbose = FALSE,
      nstart = 10,
      max_iter = 100
    )
  }, error = function(e) {
    cat(" ✗ Error:", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(result)) {
    eval <- evaluate_comprehensive(result, test_data$true_labels, test_data$x_list, 
                                   paste0("sep_", sep))
    results_all[[length(results_all)+1]] <- eval
    cat(sprintf(" ✓ ARI: %.3f, Cluster purity: %.3f\n", eval$ari, eval$cluster_purity))
  }
}

# Save comprehensive results
if (length(results_all) > 0) {
  final_df <- do.call(rbind, results_all)
  write.csv(final_df, "extensive_k_selection_results.csv", row.names = FALSE)
  cat("\n✓ Results saved to: extensive_k_selection_results.csv\n")
  
  # Print summary statistics
  cat("\n" + paste(rep("=", 80), collapse = "") + "\n")
  cat("SUMMARY STATISTICS\n")
  cat(paste(rep("=", 80), collapse = "") + "\n")
  cat(sprintf("\nNoise Detection Accuracy: Mean = %.3f, SD = %.3f\n", 
              mean(final_df$noise_accuracy, na.rm = TRUE),
              sd(final_df$noise_accuracy, na.rm = TRUE)))
  cat(sprintf("ARI (Cluster Recovery): Mean = %.3f, SD = %.3f\n", 
              mean(final_df$ari, na.rm = TRUE),
              sd(final_df$ari, na.rm = TRUE)))
  cat(sprintf("Average Iterations: %.1f\n", mean(final_df$iterations)))
  cat(sprintf("Convergence Rate: %.1f%%\n", mean(final_df$converged) * 100))
}

cat("\n✓ Extensive K Selection Tests Complete\n")