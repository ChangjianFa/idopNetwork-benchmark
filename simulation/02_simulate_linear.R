# ============================================================
# 场景B：线性 VAR(1) 模型模拟（对 WGCNA/GENIE3 更公平的场景）
# ============================================================
# X(t) = B_scaled * X(t-1) + epsilon
# B_true[i, j] 含义：物种 j 调控物种 i（j → i）
# ============================================================

# 模拟单次 VAR(1) 轨迹，返回 n_species × n_pos 矩阵
simulate_var1 <- function(B_true, n_positions, noise_sd, seed) {
  set.seed(seed)
  n <- nrow(B_true)

  # 稳定化：缩放 B 使谱半径 < 0.9，防止发散
  sr <- max(abs(eigen(B_true, only.values = TRUE)$values))
  B_scaled <- if (sr >= 1) B_true * 0.9 / sr else B_true

  X <- matrix(0.0, n, n_positions)
  X[, 1] <- abs(rnorm(n, 5, 1))

  for (t in 2:n_positions) {
    X[, t] <- pmax(B_scaled %*% X[, t - 1] + rnorm(n, 0, noise_sd), 0)
  }

  rownames(X) <- rownames(B_true)
  colnames(X) <- paste0("t", seq_len(n_positions))
  return(X)
}

# 批量生成所有重复并保存（复用 LV 的 B_true 保证可比性）
generate_all_var1 <- function(cfg, B_true) {
  replicates <- lapply(seq_len(cfg$N_REPLICATES), function(rep_id) {
    abund <- simulate_var1(B_true, cfg$N_POSITIONS, cfg$NOISE_SD,
                           seed = cfg$MASTER_SEED + rep_id * 200)
    list(rep_id = rep_id, abund = abund, scenario = "VAR1")
  })

  saveRDS(list(B_true = B_true, replicates = replicates),
          file.path(cfg$OUT_DIR, "sim_var1.rds"))
  message(sprintf("场景B：已生成 %d 次重复，保存至 sim_var1.rds", cfg$N_REPLICATES))
  invisible(list(B_true = B_true, replicates = replicates))
}
