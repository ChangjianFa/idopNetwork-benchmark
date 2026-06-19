# ============================================================
# 结果可视化（LV 场景）
# ============================================================

library(ggplot2)

load_results <- function(out_dir) {
  f <- file.path(out_dir, "final_results.rds")
  if (!file.exists(f)) stop("找不到 final_results.rds，请先运行 run_all.R")
  readRDS(f)
}

# 图1：AUROC / TPR / FPR / MCC 综合精度图（2×2，含 95% CI 误差棒）
plot_metrics_combined <- function(results_df, fig_dir) {
  df <- results_df[results_df$scenario == "LV",
                   c("method", "AUROC", "TPR", "FPR", "MCC")]

  df_long <- reshape(df,
                     varying   = c("AUROC", "TPR", "FPR", "MCC"),
                     v.names   = "value",
                     timevar   = "metric",
                     times     = c("AUROC", "TPR", "FPR", "MCC"),
                     direction = "long")
  df_long$metric <- factor(df_long$metric,
                           levels = c("AUROC", "TPR", "FPR", "MCC"))

  baselines <- data.frame(
    metric   = factor(c("AUROC", "TPR", "FPR", "MCC"),
                      levels = c("AUROC", "TPR", "FPR", "MCC")),
    baseline = c(0.5, 0.244, 0.244, 0)
  )

  ci_fun <- function(x) {
    x   <- x[!is.na(x)]
    n   <- length(x)
    m   <- mean(x)
    se  <- if (n > 1) sd(x) / sqrt(n) else 0
    data.frame(y = m, ymin = m - 1.96 * se, ymax = m + 1.96 * se)
  }

  p <- ggplot(df_long, aes(x = method, y = value, fill = method)) +
    geom_boxplot(outlier.shape = 21, alpha = 0.45, width = 0.5,
                 color = "gray30") +
    geom_jitter(width = 0.12, alpha = 0.3, size = 1.0, color = "gray50") +
    stat_summary(fun.data = ci_fun, geom = "errorbar",
                 width = 0.28, linewidth = 1.0, color = "black") +
    stat_summary(fun = mean, geom = "point",
                 shape = 23, size = 3, color = "black", fill = "white") +
    geom_hline(data = baselines, aes(yintercept = baseline),
               linetype = "dashed", color = "gray40", linewidth = 0.8,
               inherit.aes = FALSE) +
    facet_wrap(~metric, scales = "free_y", nrow = 2) +
    scale_fill_manual(values = c("idopNetwork" = "#E74C3C",
                                 "WGCNA"       = "#3498DB",
                                 "GENIE3"      = "#2ECC71")) +
    labs(title    = "LV 场景网络推断精度（30 次 Monte Carlo 重复）",
         subtitle = "箱线图 + 均值（菱形）+ 95% CI 误差棒；虚线 = 随机基准",
         x = "方法", y = "指标值", fill = "方法") +
    theme_bw(base_size = 13) +
    theme(legend.position    = "bottom",
          strip.background   = element_rect(fill = "#F0F0F0"),
          strip.text         = element_text(face = "bold"),
          plot.title         = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle      = element_text(hjust = 0.5, color = "gray40"))

  ggsave(file.path(fig_dir, "fig1_metrics_combined.pdf"), p, width = 10, height = 8)
  ggsave(file.path(fig_dir, "fig1_metrics_combined.png"), p, width = 10, height = 8,
         dpi = 300)
  message("已保存 fig1_metrics_combined")
  invisible(p)
}

# 图2：运行时间小提琴图（log10 轴，LV 场景）
plot_runtime <- function(results_df, fig_dir) {
  df <- results_df[results_df$scenario == "LV" & !is.na(results_df$runtime), ]

  p <- ggplot(df, aes(x = method, y = runtime, fill = method)) +
    geom_violin(alpha = 0.7, trim = FALSE) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9) +
    scale_y_log10(labels = function(x) {
      ifelse(x < 60, paste0(round(x, 0), " s"),
             paste0(round(x / 60, 1), " min"))
    }) +
    scale_fill_manual(values = c("idopNetwork" = "#E74C3C",
                                 "WGCNA"       = "#3498DB",
                                 "GENIE3"      = "#2ECC71")) +
    labs(title = "运行时间比较（LV 场景，log10 轴）",
         x = "方法", y = "运行时间（log10）", fill = "方法") +
    theme_bw(base_size = 13) +
    theme(legend.position  = "bottom",
          strip.background = element_rect(fill = "#F0F0F0"),
          plot.title       = element_text(hjust = 0.5, face = "bold"))

  ggsave(file.path(fig_dir, "fig2_runtime.pdf"), p, width = 7, height = 5)
  ggsave(file.path(fig_dir, "fig2_runtime.png"), p, width = 7, height = 5, dpi = 300)
  message("已保存 fig2_runtime")
  invisible(p)
}

# 图3：方向准确率 + 符号准确率柱状图（LV 场景）
plot_direction_sign <- function(results_df, fig_dir) {
  df_lv <- results_df[results_df$scenario == "LV", ]

  df_dir <- df_lv[df_lv$method %in% c("idopNetwork", "GENIE3") &
                  !is.na(df_lv$dir_acc), ]
  df_sign <- df_lv[df_lv$method == "idopNetwork" & !is.na(df_lv$sign_acc), ]

  agg_fun <- function(df, col, metric_label) {
    agg <- aggregate(df[[col]] ~ method, df,
                     FUN = function(x) c(mean = mean(x), se = sd(x) / sqrt(length(x))))
    agg <- do.call(data.frame, agg)
    colnames(agg)[2:3] <- c("mean", "se")
    agg$metric <- metric_label
    agg
  }

  plot_df <- agg_fun(df_dir, "dir_acc", "方向准确率")
  if (nrow(df_sign) > 0)
    plot_df <- rbind(plot_df, agg_fun(df_sign, "sign_acc", "符号准确率"))

  p <- ggplot(plot_df, aes(x = method, y = mean, fill = method)) +
    geom_col(alpha = 0.85, width = 0.6) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                  width = 0.25, linewidth = 0.8) +
    geom_hline(yintercept = 0.5, linetype = "dashed",
               color = "gray40", linewidth = 0.8) +
    facet_wrap(~metric, nrow = 1) +
    scale_fill_manual(values = c("idopNetwork" = "#E74C3C",
                                 "GENIE3"      = "#2ECC71")) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(title    = "方向准确率与符号准确率（LV 场景，均值 ± SE）",
         subtitle = "虚线 = 0.5（随机猜测基准）",
         x = "方法", y = "准确率", fill = "方法") +
    theme_bw(base_size = 13) +
    theme(legend.position  = "bottom",
          strip.background = element_rect(fill = "#F0F0F0"),
          strip.text       = element_text(face = "bold"),
          plot.title       = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle    = element_text(hjust = 0.5, color = "gray40"))

  ggsave(file.path(fig_dir, "fig3_direction_sign.pdf"), p, width = 8, height = 5)
  ggsave(file.path(fig_dir, "fig3_direction_sign.png"), p, width = 8, height = 5,
         dpi = 300)
  message("已保存 fig3_direction_sign")
  invisible(p)
}

# 汇总统计表
summarize_results <- function(results_df) {
  metrics <- c("AUROC", "AUPRC", "TPR", "FPR", "MCC", "dir_acc", "sign_acc", "runtime")
  df_lv   <- results_df[results_df$scenario == "LV", ]
  out <- lapply(metrics, function(m) {
    lapply(unique(df_lv$method), function(mth) {
      v <- df_lv[[m]][df_lv$method == mth]
      data.frame(method = mth, metric = m,
                 mean = round(mean(v, na.rm = TRUE), 3),
                 sd   = round(sd(v,   na.rm = TRUE), 3),
                 n    = sum(!is.na(v)),
                 stringsAsFactors = FALSE)
    })
  })
  do.call(rbind, unlist(out, recursive = FALSE))
}

# 主入口
make_all_figures <- function(cfg) {
  results_df <- load_results(cfg$OUT_DIR)
  plot_metrics_combined(results_df, cfg$FIG_DIR)
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
