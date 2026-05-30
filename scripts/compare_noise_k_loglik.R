#!/usr/bin/env Rscript

project_root <- normalizePath(getwd(), winslash = "/")
if (!file.exists(file.path(project_root, "R", "Matrix.R"))) {
	parent_root <- normalizePath(file.path(project_root, ".."), winslash = "/")
	if (file.exists(file.path(parent_root, "R", "Matrix.R"))) {
		project_root <- parent_root
	} else {
		stop("Could not locate the project root. Run this script from the project root or the scripts directory.")
	}
}

source(file.path(project_root, "R", "Matrix.R"))
source(file.path(project_root, "R", "Matrix_Utils.R"))
source(file.path(project_root, "R", "Matrix_Init.R"))
source(file.path(project_root, "R", "HC_Noise_Search_Utils.R"))
source(file.path(project_root, "R", "HC_Noise_Search.R"))
source(file.path(project_root, "R", "Matrix_Noise.R"))

set.seed(123)

simulate_noisy_matrix <- function(mean_matrix, noise_sd = 0.25) {
	noise <- matrix(rnorm(length(mean_matrix), mean = 0, sd = noise_sd), nrow = nrow(mean_matrix))
	mean_matrix + noise
}

simulate_x_list <- function(n = 30, r = 4, p = 5, noise_sd = 0.25) {
	base_mean <- matrix(seq_len(r * p), nrow = r, ncol = p)
	base_mean <- base_mean / max(base_mean)
	lapply(seq_len(n), function(i) {
		shift <- (i - 1) / max(1, n - 1)
		mean_matrix <- base_mean + shift * 0.35
		simulate_noisy_matrix(mean_matrix, noise_sd = noise_sd)
	})
}

x_list <- simulate_x_list(n = 30, r = 4, p = 5, noise_sd = 0.3)

auto_grid <- matrix_noise_hc_search_grid(
	noise_k_grid = 10^seq(-12, -2, length.out = 25),
	x_list = x_list
)

candidate_grid <- sort(unique(auto_grid))

results <- data.frame(
	noise_k = candidate_grid,
	log10_noise_k = log10(candidate_grid),
	logLik = NA_real_,
	converged = NA,
	stringsAsFactors = FALSE
)

for (i in seq_along(candidate_grid)) {
	fit <- matrix_variate_noise_fit(
		x_list = x_list,
		g = 2,
		noise_type = "hc",
		noise_k = candidate_grid[i],
		select_noise_k = FALSE,
		max_iter = 100,
		tol = 1e-06,
		nstart = 5,
		verbose = FALSE
	)
	results$logLik[i] <- tail(fit$logLik, 1)
	results$converged[i] <- isTRUE(fit$converged)
}

best_idx <- which.max(results$logLik)
best_row <- results[best_idx, , drop = FALSE]

print(results)
cat("\nBest candidate:\n")
print(best_row)

plot(
	results$log10_noise_k,
	results$logLik,
	type = "b",
	pch = 19,
	xlab = expression(log[10](noise[k])),
	ylab = "Observed-data log-likelihood",
	main = "HC noise-k scan: log-likelihood vs log10(noise_k)"
)
abline(v = best_row$log10_noise_k, col = "red", lty = 2)
points(best_row$log10_noise_k, best_row$logLik, col = "red", pch = 19)
