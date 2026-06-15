# 调试脚本：逐步测试 idopNetwork 管道，定位 $ 操作符错误
setwd("d:/data/idopnetwork")
source("simulation/config.R")
source("simulation/01_simulate_lv.R")
library(idopNetwork)

cat("=== 生成测试数据 ===\n")
B_true <- make_ground_truth(N_SPECIES, SPARSITY_FRAC, EDGE_STRENGTH, MASTER_SEED)
set.seed(MASTER_SEED)
r_vec <- runif(N_SPECIES, 0.1, 0.5)
names(r_vec) <- rownames(B_true)
abund <- simulate_lv(B_true, r_vec, N_POSITIONS, NOISE_SD, seed = MASTER_SEED + 100)
cat("abund 维度:", dim(abund), "\n")
cat("colSums 范围:", range(colSums(abund)), "\n")
cat("有零值的物种数:", sum(rowSums(abund == 0) > 0), "\n")

cat("\n=== Step 1: data.frame 转换 ===\n")
df <- data.frame(name = rownames(abund), abund, check.names = FALSE)
cat("df 维度:", dim(df), "列名头部:", head(colnames(df), 4), "\n")

cat("\n=== Step 2: data_cleaning ===\n")
x_thresh <- round(ncol(abund) * 0.5)
cat("过滤阈值 x =", x_thresh, "\n")
cleaned <- tryCatch(
  data_cleaning(df, x = x_thresh),
  error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(cleaned)) cat("cleaned 维度:", dim(cleaned), "rownames:", head(rownames(cleaned)), "\n")

cat("\n=== Step 3: power_equation_fit (单线程) ===\n")
pfit <- tryCatch(
  power_equation_fit(cleaned, n = IDOP_N_INTERP, trans = log10, thread = 1),
  error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(pfit)) {
  cat("pfit 元素:", names(pfit), "\n")
  cat("power_fit 维度:", dim(pfit$power_fit), "\n")
  cat("power_par 维度:", dim(pfit$power_par), "\n")
  cat("Time 长度:", length(pfit$Time), "\n")
}

cat("\n=== Step 4: get_interaction (物种1) ===\n")
rel1 <- tryCatch(
  get_interaction(pfit$original_data, 1),
  error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(rel1)) {
  cat("ind.name:", rel1$ind.name, "\n")
  cat("dep.name:", paste(rel1$dep.name, collapse=","), "\n")
  cat("coefficient:", paste(round(rel1$coefficient, 3), collapse=","), "\n")
}

cat("\n=== Step 5: qdODE_all (物种1, maxit=50) ===\n")
ode1 <- tryCatch(
  qdODE_all(result = pfit, relationship = list(rel1), i = 1, maxit = 50),
  error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(ode1) && !all(is.na(ode1))) {
  cat("ode1 元素:", names(ode1), "\n")
  cat("fit 维度:", dim(ode1$fit), "\n")
} else {
  cat("ode1 返回 NA 或 NULL\n")
}

cat("\n=== Step 6: network_conversion ===\n")
if (!is.null(ode1) && !all(is.na(ode1))) {
  nc <- tryCatch(
    network_conversion(ode1),
    error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(nc)) cat("edges:\n"); print(nc$edge)
}

cat("\n=== 调试完成 ===\n")
