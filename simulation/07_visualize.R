# ============================================================
# 结果可视化
# ============================================================

library(ggplot2)

# 加载并合并所有评估结果
load_results <- function(out_dir) {
  f <- file.path(out_dir, "final_results.rds")
  if (!file.exists(f)) stop("找不到 final_results.rds，请先运行 run_all.R")
  readRDS(f)
}

# 图1：AUROC 和 AUPRC 箱线图（按场景分面）
plot_auroc_auprc <- function(results_df, fig_dir) {
  # 随机基准线
  positive_rate <- 0.25 * 19 / 20  # ≈ 0.2375

  df_long <- reshape(results_df[, c("method", "scenario", "AUROC", "AUPRC")],
                     varying   = c("AUROC", "AUPRC"),
                     v.names   = "value",
                     timevar   = "metric",
                     times     = c("AUROC", "AUPRC"),
                     direction = "long")

  # 添加基准线数据
  baselines <- data.frame(
    metric   = c("AUROC", "AUPRC"),
    baseline = c(0.5, positive_rate)
  )

  p <- ggplot(df_long, aes(x = method, y = value, fill = method)) +
    geom_boxplot(outlier.shape = 21, alpha = 0.8, width = 0.6) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 1.2) +
    geom_hline(data = baselines, aes(yintercept = baseline),
               linetype = "dashed", color = "gray40", linewidth = 0.8) +
    facet_grid(metric ~ scenario, scales = "free_y") +
    scale_fill_manual(values = c("idopNetwork" = "#E74C3C",
                                 "WGCNA"       = "#3498DB",
                                 "GENIE3"      = "#2ECC71")) +
    labs(title = "网络推断性能比较（AUROC / AUPRC）",
         subtitle = "虚线 = 随机基准；场景A = LV模拟；场景B = 线性VAR(1)模拟",
         x = "方法", y = "指标值", fill = "方法") +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom",
          strip.background = element_rect(fill = "#F0F0F0"),
          plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave(file.path(fig_dir, "fig1_auroc_auprc.pdf"), p,
         width = 10, height = 7)
  ggsave(file.path(fig_dir, "fig1_auroc_auprc.png"), p,
         width = 10, height = 7, dpi = 300)
  message("已保存 fig1_auroc_auprc")
  invisible(p)
}

# 图2：运行时间小提琴图（log10 轴）
plot_runtime <- function(results_df, fig_dir) {
  df <- results_df[!is.na(results_df$runtime), ]

  p <- ggplot(df, aes(x = method, y = runtime, fill = method)) +
    geom_violin(alpha = 0.7, trim = FALSE) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9) +
    scale_y_log10(labels = function(x) paste0(round(x/60, 1), " min")) +
    facet_wrap(~scenario) +
    scale_fill_manual(values = c("idopNetwork" = "#E74C3C",
                                 "WGCNA"       = "#3498DB",
                                 "GENIE3"      = "#2ECC71")) +
    labs(title = "运行时间比较（log10 轴）",
         x = "方法", y = "运行时间（秒，log10）", fill = "方法") +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom",
          strip.background = element_rect(fill = "#F0F0F0"),
          plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave(file.path(fig_dir, "fig2_runtime.pdf"), p, width = 9, height = 5)
  ggsave(file.path(fig_dir, "fig2_runtime.png"), p, width = 9, height = 5, dpi = 300)
  message("已保存 fig2_runtime")
  invisible(p)
}

# 图3：方向准确率 + 符号准确率柱状图
plot_direction_sign <- function(results_df, fig_dir) {
  df_dir <- results_df[results_df$method %in% c("idopNetwork", "GENIE3") &
                       !is.na(results_df$dir_acc), ]
  df_sign <- results_df[results_df$method == "idopNetwork" &
                        !is.na(results_df$sign_acc), ]

  # 计算均值和标准误
  agg_dir <- aggregate(dir_acc ~ method + scenario, df_dir,
                       FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))
  agg_dir <- do.call(data.frame, agg_dir)
  colnames(agg_dir)[3:4] <- c("mean", "se")
  agg_dir$metric <- "方向准确率"

  if (nrow(df_sign) > 0) {
    agg_sign <- aggregate(sign_acc ~ method + scenario, df_sign,
                          FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))
    agg_sign <- do.call(data.frame, agg_sign)
    colnames(agg_sign)[3:4] <- c("mean", "se")
    agg_sign$metric <- "符号准确率"
    plot_df <- rbind(agg_dir, agg_sign)
  } else {
    plot_df <- agg_dir
  }

  p <- ggplot(plot_df, aes(x = method, y = mean, fill = method)) +
    geom_col(alpha = 0.85, width = 0.6) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                  width = 0.25, linewidth = 0.8) +
    geom_hline(yintercept = 0.5, linetype = "dashed",
               color = "gray40", linewidth = 0.8) +
    facet_grid(metric ~ scenario) +
    scale_fill_manual(values = c("idopNetwork" = "#E74C3C",
                                 "GENIE3"      = "#2ECC71")) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(title = "方向准确率与符号准确率（均值 ± 标准误）",
         subtitle = "虚线 = 0.5（随机猜测基准）",
         x = "方法", y = "准确率", fill = "方法") +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom",
          strip.background = element_rect(fill = "#F0F0F0"),
          plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave(file.path(fig_dir, "fig3_direction_sign.pdf"), p, width = 9, height = 6)
  ggsave(file.path(fig_dir, "fig3_direction_sign.png"), p, width = 9, height = 6, dpi = 300)
  message("已保存 fig3_direction_sign")
  invisible(p)
}

# 汇总统计表
summarize_results <- function(results_df) {
  metrics <- c("AUROC", "AUPRC", "TPR", "FPR", "MCC", "dir_acc", "sign_acc", "runtime")
  out <- lapply(metrics, function(m) {
    vals <- results_df[[m]]
    combos <- unique(results_df[, c("method", "scenario")])
    rows <- lapply(seq_len(nrow(combos)), function(i) {
      mth <- combos$method[i]; scn <- combos$scenario[i]
      v <- vals[results_df$method == mth & results_df$scenario == scn]
      data.frame(method   = mth,
                 scenario = scn,
                 metric   = m,
                 mean     = round(mean(v, na.rm = TRUE), 3),
                 sd       = round(sd(v,   na.rm = TRUE), 3),
                 n        = sum(!is.na(v)),
                 stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
  })
  do.call(rbind, out)
}

# 主入口
make_all_figures <- function(cfg) {
  results_df <- load_results(cfg$OUT_DIR)
  plot_auroc_auprc(results_df, cfg$FIG_DIR)
  plot_runtime(results_df, cfg$FIG_DIR)
  plot_direction_sign(results_df, cfg$FIG_DIR)

  summary_tbl <- summarize_results(results_df)
  write.csv(summary_tbl,
            file.path(cfg$OUT_DIR, "summary_table.csv"),
            row.names = FALSE)
  message("汇总统计表已保存为 summary_table.csv")
  print(summary_tbl)
  invisible(summary_tbl)
}
