# ============================================================
# 参数敏感性扫描：噪声方差 × 样本量
# ============================================================
# 扫描线1：NOISE_SD in [0.01, 0.03, 0.05, 0.07, 0.10]，固定 N_POSITIONS=30
# 扫描线2：N_POSITIONS in [15, 20, 30, 50, 80]，固定 NOISE_SD=0.05
# 每个参数值运行 20 次 Monte Carlo 重复
# ============================================================

source("simulation/config.R")
source("simulation/01_simulate_lv.R")
source("simulation/03_run_idopnetwork.R")
source("simulation/04_run_wgcna.R")
source("simulation/05_run_genie3.R")
source("simulation/06_evaluate.R")

if (!requireNamespace("idopNetwork", quietly = TRUE))
  devtools::install_local(PKG_DIR, quiet = TRUE)
library(idopNetwork)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

N_SWEEP_REPS <- 20L
SWEEP_SEED   <- 123L

# 固定真实网络（与 run_all.R 相同）
set.seed(MASTER_SEED)
B_true <- make_ground_truth(N_SPECIES, SPARSITY_FRAC, EDGE_STRENGTH, MASTER_SEED)
r_vec  <- runif(N_SPECIES, 0.1, 0.5)

cfg_base <- list(
  IDOP_N_INTERP    = IDOP_N_INTERP,
  IDOP_MAXIT       = IDOP_MAXIT,
  WGCNA_SOFT_POWER = WGCNA_SOFT_POWER,
  GENIE3_NTREES    = GENIE3_NTREES
)

# 运行单次模拟并评估
run_one <- function(noise_sd, n_pos, rep_id, seed_val) {
  abund <- tryCatch(
    simulate_lv(B_true, r_vec, n_pos, noise_sd, seed = seed_val),
    error = function(e) NULL
  )
  if (is.null(abund)) {
    return(data.frame(
      method = c("idopNetwork", "WGCNA", "GENIE3"),
      AUROC = NA, AUPRC = NA, TPR = NA, FPR = NA, MCC = NA,
      dir_acc = NA, sign_acc = NA, runtime = NA, status = "simulate_lv failed",
      noise_sd = noise_sd, n_pos = n_pos, rep_id = rep_id,
      stringsAsFactors = FALSE
    ))
  }
  rownames(abund) <- paste0("sp", seq_len(N_SPECIES))

  r_idop  <- run_idopnetwork(abund, cfg_base)
  r_wgcna <- run_wgcna(abund, cfg_base)
  r_genie <- run_genie3(abund, cfg_base)

  df <- rbind(
    evaluate_method(r_idop,  B_true),
    evaluate_method(r_wgcna, B_true),
    evaluate_method(r_genie, B_true)
  )
  df$noise_sd <- noise_sd
  df$n_pos    <- n_pos
  df$rep_id   <- rep_id
  df
}

# ---- 扫描线1：噪声方差（固定 N_POSITIONS=30）----
noise_vals <- c(0.01, 0.03, 0.05, 0.07, 0.10)
noise_rows <- vector("list", length(noise_vals) * N_SWEEP_REPS)
idx <- 0L

message("\n========== 噪声方差扫描（N_POSITIONS=30）==========")
for (nv in noise_vals) {
  for (r in seq_len(N_SWEEP_REPS)) {
    idx <- idx + 1L
    message(sprintf("[噪声扫描] noise=%.2f  rep=%d/%d  (总进度 %d/%d)",
                    nv, r, N_SWEEP_REPS, idx, length(noise_vals) * N_SWEEP_REPS))
    noise_rows[[idx]] <- run_one(nv, 30L, r, SWEEP_SEED + idx)
  }
  # 每个噪声值完成后保存一次检查点
  saveRDS(noise_rows[seq_len(idx)],
          file.path(OUT_DIR, "sweep_noise_checkpoint.rds"))
}

# ---- 扫描线2：样本量（固定 NOISE_SD=0.05）----
npos_vals <- c(15L, 20L, 30L, 50L, 80L)
npos_rows <- vector("list", length(npos_vals) * N_SWEEP_REPS)
idx2 <- 0L

message("\n========== 样本量扫描（NOISE_SD=0.05）==========")
for (np in npos_vals) {
  for (r in seq_len(N_SWEEP_REPS)) {
    idx2 <- idx2 + 1L
    message(sprintf("[样本量扫描] n_pos=%d  rep=%d/%d  (总进度 %d/%d)",
                    np, r, N_SWEEP_REPS, idx2, length(npos_vals) * N_SWEEP_REPS))
    npos_rows[[idx2]] <- run_one(0.05, np, r, SWEEP_SEED + 10000L + idx2)
  }
  saveRDS(npos_rows[seq_len(idx2)],
          file.path(OUT_DIR, "sweep_npos_checkpoint.rds"))
}

# ---- 合并并保存 ----
noise_df        <- do.call(rbind, noise_rows)
noise_df$sweep  <- "noise"
npos_df         <- do.call(rbind, npos_rows)
npos_df$sweep   <- "npos"

sweep_results <- rbind(noise_df, npos_df)
saveRDS(sweep_results, file.path(OUT_DIR, "sweep_results.rds"))
write.csv(sweep_results, file.path(OUT_DIR, "sweep_results.csv"), row.names = FALSE)
message("\n扫描完成！结果保存至 sweep_results.csv")
message(sprintf("共 %d 行记录", nrow(sweep_results)))
