# ============================================================
# 主控脚本：运行完整模拟实验
# ============================================================
# 建议先以 Pilot 模式验证全流程（约 30 分钟）：
#   source("simulation/run_all.R")  # 使用默认参数
# 再修改 config.R 为完整参数运行（约 8-20 小时）
# ============================================================

# ---- 0. 加载配置和依赖 ----
source("simulation/config.R")
source("simulation/01_simulate_lv.R")
source("simulation/03_run_idopnetwork.R")
source("simulation/04_run_wgcna.R")
source("simulation/05_run_genie3.R")
source("simulation/06_evaluate.R")

# 安装并加载 idopNetwork（从本地子文件夹）
if (!requireNamespace("idopNetwork", quietly = TRUE)) {
  devtools::install_local(PKG_DIR, quiet = TRUE)
}
library(idopNetwork)

# 确保输出目录存在
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cfg <- list(
  N_SPECIES        = N_SPECIES,
  N_POSITIONS      = N_POSITIONS,
  N_REPLICATES     = N_REPLICATES,
  NOISE_SD         = NOISE_SD,
  SPARSITY_FRAC    = SPARSITY_FRAC,
  EDGE_STRENGTH    = EDGE_STRENGTH,
  IDOP_THREADS     = IDOP_THREADS,
  IDOP_MAXIT       = IDOP_MAXIT,
  IDOP_N_INTERP    = IDOP_N_INTERP,
  WGCNA_SOFT_POWER = WGCNA_SOFT_POWER,
  GENIE3_NTREES    = GENIE3_NTREES,
  MASTER_SEED      = MASTER_SEED,
  OUT_DIR          = OUT_DIR,
  FIG_DIR          = FIG_DIR
)

# ---- 1. 生成模拟数据 ----
message("\n========== 生成模拟数据 ==========")

sim_lv_data <- generate_all_lv(cfg)
B_true <- sim_lv_data$B_true

message(sprintf(
  "真实网络：%d 物种，%d 条真实边（稀疏度 %.1f%%）",
  nrow(B_true),
  sum(B_true[row(B_true) != col(B_true)] != 0),
  mean(B_true[row(B_true) != col(B_true)] != 0) * 100
))

# ---- 2. 主循环：逐重复运行三种方法 ----
checkpoint_file <- file.path(OUT_DIR, "checkpoint.rds")

# 支持断点续跑：若 checkpoint 存在则加载
if (file.exists(checkpoint_file)) {
  all_eval <- readRDS(checkpoint_file)
  start_rep <- length(all_eval) + 1L
  message(sprintf("从断点恢复：从第 %d 次重复继续", start_rep))
} else {
  all_eval <- list()
  start_rep <- 1L
}

for (rep_id in seq(start_rep, N_REPLICATES)) {
  message(sprintf("\n========== 重复 %d / %d ==========", rep_id, N_REPLICATES))

  abund <- sim_lv_data$replicates[[rep_id]]$abund

  if (is.null(abund)) {
    message("  [LV] 模拟失败，跳过")
    next
  }

  message("  [LV] 运行 idopNetwork ...")
  r_idop <- run_idopnetwork(abund, cfg)
  message(sprintf("    状态: %s  用时: %.1f 秒", r_idop$status,
                  ifelse(is.na(r_idop$runtime), 0, r_idop$runtime)))

  message("  [LV] 运行 WGCNA ...")
  r_wgcna <- run_wgcna(abund, cfg)
  message(sprintf("    状态: %s  用时: %.1f 秒", r_wgcna$status,
                  ifelse(is.na(r_wgcna$runtime), 0, r_wgcna$runtime)))

  message("  [LV] 运行 GENIE3 ...")
  r_genie <- run_genie3(abund, cfg)
  message(sprintf("    状态: %s  用时: %.1f 秒", r_genie$status,
                  ifelse(is.na(r_genie$runtime), 0, r_genie$runtime)))

  eval_df <- rbind(
    evaluate_method(r_idop,  B_true),
    evaluate_method(r_wgcna, B_true),
    evaluate_method(r_genie, B_true)
  )
  eval_df$rep_id   <- rep_id
  eval_df$scenario <- "LV"

  message(sprintf("  AUROC: idop=%.3f  WGCNA=%.3f  GENIE3=%.3f",
                  eval_df$AUROC[eval_df$method == "idopNetwork"],
                  eval_df$AUROC[eval_df$method == "WGCNA"],
                  eval_df$AUROC[eval_df$method == "GENIE3"]))

  all_eval[[rep_id]] <- eval_df
  saveRDS(all_eval, checkpoint_file)  # 每次重复后保存断点
}

# ---- 3. 汇总并保存最终结果 ----
message("\n========== 汇总结果 ==========")
results_df <- do.call(rbind, all_eval)
saveRDS(results_df, file.path(OUT_DIR, "final_results.rds"))
write.csv(results_df, file.path(OUT_DIR, "final_results.csv"), row.names = FALSE)
message("结果已保存：final_results.rds / final_results.csv")

# ---- 4. 生成图表 ----
message("\n========== 生成图表 ==========")
source("simulation/07_visualize.R")
make_all_figures(cfg)

message("\n========== 实验完成 ==========")
print(aggregate(cbind(AUROC, AUPRC, TPR, MCC, dir_acc) ~ method,
                data = results_df,
                FUN  = function(x) round(mean(x, na.rm = TRUE), 3)))
