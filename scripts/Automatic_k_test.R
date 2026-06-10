# test_automatic_k_selection.R
# Tests the automatic noise parameter (k) selection scheme for HC noise type

# Load the algorithm files
source("../R/Matrix_Init.R")
source("../R/Matrix.R")
source("../R/Matrix_KS_Score.R")
source("../R/Matrix_Noise.R")
source("../R/Matrix_Noise_BR.R")
source("../R/Matrix_Utils.R")

# Helper function to generate synthetic matrix-variate data
generate_test_data <- function(n_per_group, rows, cols, 
                               group_means, row_covs, col_covs,
                               noise_prop = 0) {
  total_n <- sum(n_per_group)
  g <- length(group_means)
  
  x_list <- vector("list", total_n)
  true_labels <- integer(total_n)
  noise_idx <- sample(total_n, floor(total_n * noise_prop))
  
  idx <- 1
  for (group in 1:g) {
    n_group <- n_per_group[group]
    for (i in 1:n_group) {
      # Generate matrix from matrix-variate normal
      Z <- matrix(rnorm(rows * cols), rows, cols)
      eigen_row <- eigen(row_covs[[group]])
      eigen_col <- eigen(col_covs[[group]])
      
      A <- eigen_row$vectors %*% diag(sqrt(eigen_row$values)) %*% t(eigen_row$vectors)
      B <- eigen_col$vectors %*% diag(sqrt(eigen_col$values)) %*% t(eigen_col$vectors)
      
      x_list[[idx]] <- group_means[[group]] + A %*% Z %*% B
      true_labels[idx] <- group
      idx <- idx + 1
    }
  }
  
  # Add noise observations
  if (noise_prop > 0) {
    for (i in noise_idx) {
      # Generate uniform noise within range of data
      all_data <- do.call(rbind, lapply(x_list, as.vector))
      data_range <- range(all_data, na.rm = TRUE)
      x_list[[i]] <- matrix(runif(rows * cols, data_range[1], data_range[2]), rows, cols)
      true_labels[i] <- 0  # 0 indicates noise
    }
  }
  
  list(x_list = x_list, true_labels = true_labels)
}

# Function to save results to CSV
save_results_to_csv <- function(results_list, filename) {
  # Convert results to data frame
  df <- data.frame(
    test_name = character(),
    selected_k = numeric(),
    noise_proportion = numeric(),
    iterations = integer(),
    converged = logical(),
    ks_statistic = numeric(),
    ks_pvalue = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (name in names(results_list)) {
    result <- results_list[[name]]
    if (!is.null(result$k_selection)) {
      df <- rbind(df, data.frame(
        test_name = name,
        selected_k = result$k_selection$selected_k,
        noise_proportion = result$noise$pi,
        iterations = result$iterations,
        converged = result$converged,
        ks_statistic = min(result$k_selection$ks_scores, na.rm = TRUE),
        ks_pvalue = max(result$k_selection$ks_pvalues, na.rm = TRUE),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  write.csv(df, filename, row.names = FALSE)
  cat("Results saved to:", filename, "\n")
  return(df)
}

# Save detailed grid search results
save_grid_results <- function(result, filename) {
  if (!is.null(result$k_selection)) {
    grid_df <- data.frame(
      k = result$k_selection$k_grid,
      ks_score = result$k_selection$ks_scores,
      ks_pvalue = result$k_selection$ks_pvalues,
      n_used = result$k_selection$n_used
    )
    write.csv(grid_df, filename, row.names = FALSE)
    cat("Grid search results saved to:", filename, "\n")
    return(grid_df)
  }
  return(NULL)
}

cat("\n========================================\n")
cat("TEST 1: Automatic K Selection\n")
cat("========================================\n")

# Store results
all_test_results <- list()
grid_results <- list()

# Test 1.1: Basic functionality
cat("\n--- Test 1.1: Basic functionality with well-separated clusters ---\n")
set.seed(123)

rows <- 4
cols <- 3
g <- 2

group_means <- list(
  matrix(0, rows, cols),
  matrix(2, rows, cols)
)

row_covs <- list(diag(rows), diag(rows))
col_covs <- list(diag(cols), diag(cols))

test_data <- generate_test_data(
  n_per_group = c(30, 30),
  rows = rows, cols = cols,
  group_means = group_means,
  row_covs = row_covs,
  col_covs = col_covs,
  noise_prop = 0
)

cat("Running automatic k selection...\n")
result1 <- matrix_variate_noise_fit(
  x_list = test_data$x_list,
  g = g,
  noise_type = "hc",
  estimate_k = TRUE,
  adaptive_grid = TRUE,
  verbose = TRUE,
  nstart = 10,
  max_iter = 100
)

all_test_results[["basic_clusters"]] <- result1
grid_results[["basic_clusters"]] <- save_grid_results(result1, "grid_search_basic.csv")

cat("\n--- Results for Test 1.1 ---\n")
if (!is.null(result1$k_selection)) {
  cat("✓ Selected k:", result1$k_selection$selected_k, "\n")
  cat("✓ KS scores range:", range(result1$k_selection$ks_scores), "\n")
  cat("✓ Grid size:", length(result1$k_selection$k_grid), "\n")
  cat("✓ n_used (non-list check):", class(result1$k_selection$n_used), "\n")
  # Check if n_used is a list and convert if needed
  if (is.list(result1$k_selection$n_used)) {
    cat("⚠ n_used is a list, converting to numeric...\n")
    result1$k_selection$n_used <- unlist(result1$k_selection$n_used)
  }
  cat("✓ n_used values:", paste(result1$k_selection$n_used[1:min(5, length(result1$k_selection$n_used))], collapse=", "), "\n")
}
cat("✓ Noise proportion:", result1$noise$pi, "\n")
cat("✓ Iterations:", result1$iterations, "\n")
cat("✓ Converged:", result1$converged, "\n")
cat("✓ Cluster assignments range:", range(result1$cluster), "\n")

# Test 1.2: With noise
cat("\n--- Test 1.2: Data with 10% noise ---\n")
set.seed(456)

test_data_noisy <- generate_test_data(
  n_per_group = c(50, 50),
  rows = rows, cols = cols,
  group_means = group_means,
  row_covs = row_covs,
  col_covs = col_covs,
  noise_prop = 0.1
)

cat("Running automatic k selection on noisy data...\n")
result2 <- matrix_variate_noise_fit(
  x_list = test_data_noisy$x_list,
  g = g,
  noise_type = "hc",
  estimate_k = TRUE,
  adaptive_grid = TRUE,
  verbose = TRUE,
  nstart = 10,
  max_iter = 100
)

all_test_results[["noisy_data"]] <- result2
grid_results[["noisy_data"]] <- save_grid_results(result2, "grid_search_noisy.csv")

cat("\n--- Results for Test 1.2 ---\n")
cat("✓ Noise proportion (true=0.1):", result2$noise$pi, "\n")
cat("✓ Selected k:", result2$k_selection$selected_k, "\n")
cat("✓ Converged:", result2$converged, "\n")

# Test 1.3: Different grid types
cat("\n--- Test 1.3: Adaptive vs Custom grid ---\n")
set.seed(789)

rows <- 3
cols <- 3

test_data_small <- generate_test_data(
  n_per_group = c(40, 40),
  rows = rows, cols = cols,
  group_means = list(matrix(0, rows, cols), matrix(3, rows, cols)),
  row_covs = list(diag(rows), diag(rows)),
  col_covs = list(diag(cols), diag(cols)),
  noise_prop = 0.05
)

cat("Running with adaptive grid...\n")
result_adaptive <- matrix_variate_noise_fit(
  x_list = test_data_small$x_list,
  g = 2,
  noise_type = "hc",
  estimate_k = TRUE,
  adaptive_grid = TRUE,
  verbose = FALSE,
  nstart = 5,
  max_iter = 50
)

custom_grid <- 10^seq(-12, -2, length.out = 20)
cat("Running with custom grid...\n")
result_custom <- matrix_variate_noise_fit(
  x_list = test_data_small$x_list,
  g = 2,
  noise_type = "hc",
  estimate_k = TRUE,
  k_grid = custom_grid,
  adaptive_grid = FALSE,
  verbose = FALSE,
  nstart = 5,
  max_iter = 50
)

all_test_results[["adaptive_grid"]] <- result_adaptive
all_test_results[["custom_grid"]] <- result_custom

cat("\nAdaptive grid selected k:", result_adaptive$k_selection$selected_k, "\n")
cat("Custom grid selected k:", result_custom$k_selection$selected_k, "\n")

# Save all results to CSV
cat("\n--- Saving Results ---\n")
results_df <- save_results_to_csv(all_test_results, "automatic_k_selection_results.csv")

# Create summary table
cat("\n--- Summary Table ---\n")
print(results_df)

cat("\n=== Automatic K Selection Tests Complete ===\n")