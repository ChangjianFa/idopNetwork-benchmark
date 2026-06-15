setwd("d:/data/idopnetwork")
source("simulation/config.R")
source("simulation/01_simulate_lv.R")
source("simulation/03_run_idopnetwork.R")
library(idopNetwork)

B_true <- make_ground_truth(N_SPECIES, SPARSITY_FRAC, EDGE_STRENGTH, MASTER_SEED)
set.seed(MASTER_SEED)
r_vec <- runif(N_SPECIES, 0.1, 0.5)
names(r_vec) <- rownames(B_true)
abund <- simulate_lv(B_true, r_vec, N_POSITIONS, NOISE_SD, seed = MASTER_SEED + 100)

cat("Running run_idopnetwork (maxit=200)...\n")
cfg_test <- list(IDOP_N_INTERP=30L, IDOP_THREADS=1L, IDOP_MAXIT=200L)
res <- run_idopnetwork(abund, cfg_test)
cat("Status:", res$status, "\n")
cat("Runtime:", round(res$runtime/60, 1), "min\n")
if (!is.null(res$score_mat)) {
  cat("score_mat dim:", dim(res$score_mat), "\n")
  cat("Non-zero edges:", sum(res$score_mat > 0), "\n")
}
