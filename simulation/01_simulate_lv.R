# ============================================================
# 场景A：Lotka-Volterra ODE 模拟数据生成
# ============================================================
# B_true[i, j] 含义：物种 j 调控物种 i（j → i）
# 对角线：负自调控（保证系统稳定）
# ============================================================

library(deSolve)

# 生成稀疏有向真实网络交互矩阵
make_ground_truth <- function(n, sparsity, strength_range, seed) {
  set.seed(seed)
  B <- matrix(0, n, n)
  n_edges <- round(n * (n - 1) * sparsity)
  off_diag_idx <- which(row(B) != col(B))
  selected <- sample(off_diag_idx, n_edges)
  B[selected] <- runif(n_edges, strength_range[1], strength_range[2])
  diag(B) <- runif(n, -0.5, -0.1)  # 负自调控确保局部稳定性
  rownames(B) <- colnames(B) <- paste0("sp", seq_len(n))
  return(B)
}

# Lotka-Volterra ODE 右端项
lv_ode <- function(t, state, parms) {
  x <- pmax(state, 0)  # 非负约束
  dx <- parms$r * x + as.vector(parms$B %*% x) * x
  list(dx)
}

# 模拟单次 LV 轨迹，返回 n_species × n_pos 矩阵
simulate_lv <- function(B, r_vec, n_pos, noise_sd, seed) {
  set.seed(seed)
  n <- nrow(B)
  x0 <- runif(n, 0.5, 2.0)
  names(x0) <- rownames(B)

  # 在时间轴 [0.1, n_pos] 上取 n_pos 个等距点
  times <- seq(0.1, n_pos, length.out = n_pos)
  out <- tryCatch(
    deSolve::ode(y = x0, times = times, func = lv_ode,
                 parms = list(B = B, r = r_vec), method = "lsoda"),
    error = function(e) NULL
  )
  if (is.null(out)) return(NULL)

  abund <- t(out[, -1])   # 行=物种，列=时间点
  abund[abund < 0] <- 0

  # 添加 log-normal 乘性噪声（模拟测序计数误差）
  noise_factor <- exp(matrix(rnorm(length(abund), 0, noise_sd), nrow(abund)))
  abund_noisy <- abund * noise_factor
  abund_noisy[abund_noisy < 0] <- 0

  rownames(abund_noisy) <- rownames(B)
  colnames(abund_noisy) <- paste0("t", seq_len(n_pos))
  return(abund_noisy)
}

# 批量生成所有重复的模拟数据并保存
generate_all_lv <- function(cfg) {
  B_true <- make_ground_truth(cfg$N_SPECIES, cfg$SPARSITY_FRAC,
                               cfg$EDGE_STRENGTH, cfg$MASTER_SEED)
  set.seed(cfg$MASTER_SEED)
  r_vec <- runif(cfg$N_SPECIES, 0.1, 0.5)  # 全部正增长率，确保 colSums 沿梯度递增
  names(r_vec) <- rownames(B_true)

  replicates <- lapply(seq_len(cfg$N_REPLICATES), function(rep_id) {
    abund <- simulate_lv(B_true, r_vec, cfg$N_POSITIONS, cfg$NOISE_SD,
                         seed = cfg$MASTER_SEED + rep_id * 100)
    list(rep_id = rep_id, abund = abund, scenario = "LV")
  })

  saveRDS(list(B_true = B_true, r_vec = r_vec, replicates = replicates),
          file.path(cfg$OUT_DIR, "sim_lv.rds"))
  message(sprintf("场景A：已生成 %d 次重复，保存至 sim_lv.rds", cfg$N_REPLICATES))
  invisible(list(B_true = B_true, r_vec = r_vec, replicates = replicates))
}
