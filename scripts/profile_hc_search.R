# Profiling helper for HC noise search (serial vs parallel)
set.seed(123)

x_list <- replicate(30, matrix(rnorm(20), 4, 5), simplify = FALSE)

cat("Serial run timing:\n")
print(system.time({
  fit_serial <- matrix_variate_noise_fit(x_list, g = 2, select_noise_k = TRUE, use_parallel = FALSE, verbose = FALSE)
}))

cat("Parallel run timing (2 cores):\n")
print(system.time({
  fit_parallel <- matrix_variate_noise_fit(x_list, g = 2, select_noise_k = TRUE, use_parallel = TRUE, n_cores = 2, verbose = FALSE)
}))

cat("Selected k (serial):", fit_serial$noise$search$selected_k, "\n")
if (!is.null(fit_parallel)) cat("Selected k (parallel):", fit_parallel$noise$search$selected_k, "\n")
