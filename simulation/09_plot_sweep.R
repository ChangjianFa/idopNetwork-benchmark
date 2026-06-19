# ============================================================
# 参数敏感性分析可视化
# 2×2 布局：行=AUROC/MCC，列=噪声方差/样本量
# 均值折线 + 95% CI 带状区间
# ============================================================

library(ggplot2)

plot_sweep <- function(cfg) {
  rds_path <- file.path(cfg$OUT_DIR, "sweep_results.rds")
  if (!file.exists(rds_path))
    stop("找不到 sweep_results.rds，请先运行 run_sweep.R")

  df <- readRDS(rds_path)

  method_colors <- c("idopNetwork" = "#E74C3C",
                     "WGCNA"       = "#3498DB",
                     "GENIE3"      = "#2ECC71")

  # 95% CI 聚合函数
  agg_ci <- function(sub_df, metric, x_col) {
    sub_df <- sub_df[!is.na(sub_df[[metric]]), ]
    do.call(rbind, lapply(unique(sub_df$method), function(mth) {
      do.call(rbind, lapply(unique(sub_df[[x_col]]), function(xv) {
        v <- sub_df[[metric]][sub_df$method == mth & sub_df[[x_col]] == xv]
        n <- length(v)
        if (n < 2) return(NULL)
        m  <- mean(v)
        se <- sd(v) / sqrt(n)
        data.frame(method = mth, x = xv,
                   mean = m, lo = m - 1.96 * se, hi = m + 1.96 * se,
                   stringsAsFactors = FALSE)
      }))
    }))
  }

  make_panel <- function(agg, x_label, y_label, baseline, x_log = FALSE) {
    p <- ggplot(agg, aes(x = x, y = mean, color = method, fill = method)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
      geom_line(linewidth = 1.0) +
      geom_point(size = 2.5) +
      geom_hline(yintercept = baseline, linetype = "dashed",
                 color = "gray40", linewidth = 0.8) +
      scale_color_manual(values = method_colors) +
      scale_fill_manual(values  = method_colors) +
      labs(x = x_label, y = y_label, color = "方法", fill = "方法") +
      theme_bw(base_size = 12) +
      theme(legend.position  = "none",
            strip.background = element_rect(fill = "#F0F0F0"),
            panel.grid.minor = element_blank())
    if (x_log) p <- p + scale_x_log10()
    p
  }

  # 噪声扫描数据
  df_noise <- df[df$sweep == "noise", ]
  # 样本量扫描数据
  df_npos  <- df[df$sweep == "npos", ]

  # 聚合
  ag_noise_auroc <- agg_ci(df_noise, "AUROC", "noise_sd")
  ag_noise_mcc   <- agg_ci(df_noise, "MCC",   "noise_sd")
  ag_npos_auroc  <- agg_ci(df_npos,  "AUROC", "n_pos")
  ag_npos_mcc    <- agg_ci(df_npos,  "MCC",   "n_pos")

  # 四个面板
  p_noise_auroc <- make_panel(ag_noise_auroc,
                              "噪声标准差 (NOISE_SD)", "AUROC", 0.5)
  p_npos_auroc  <- make_panel(ag_npos_auroc,
                              "样本量 (N_POSITIONS)",  "AUROC", 0.5)
  p_noise_mcc   <- make_panel(ag_noise_mcc,
                              "噪声标准差 (NOISE_SD)", "MCC",   0.0)
  p_npos_mcc    <- make_panel(ag_npos_mcc,
                              "样本量 (N_POSITIONS)",  "MCC",   0.0)

  # 拼图（用 patchwork）
  if (!requireNamespace("patchwork", quietly = TRUE))
    install.packages("patchwork")
  library(patchwork)

  # 提取图例（从任意面板）
  p_legend <- ggplot(
    data.frame(method = names(method_colors),
               x = 1, y = 1),
    aes(x = x, y = y, color = method, fill = method)
  ) +
    geom_point(size = 3) +
    scale_color_manual(values = method_colors, name = "方法") +
    scale_fill_manual(values  = method_colors, name = "方法") +
    theme_void() +
    theme(legend.position = "bottom",
          legend.title    = element_text(face = "bold"))
  legend_grob <- cowplot::get_legend(p_legend +
                                       theme(legend.direction = "horizontal"))

  combined <- (p_noise_auroc | p_npos_auroc) /
              (p_noise_mcc   | p_npos_mcc) +
    plot_annotation(
      title    = "参数敏感性分析（LV 场景，idopNetwork vs WGCNA vs GENIE3）",
      subtitle = sprintf("均值 ± 95%% CI（每点 20 次 Monte Carlo 重复）；虚线 = 随机基准"),
      theme = theme(
        plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 11)
      )
    ) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  out_png <- file.path(cfg$FIG_DIR, "fig4_sensitivity.png")
  out_pdf <- file.path(cfg$FIG_DIR, "fig4_sensitivity.pdf")
  ggsave(out_png, combined, width = 12, height = 8, dpi = 300)
  ggsave(out_pdf, combined, width = 12, height = 8)
  message("已保存 fig4_sensitivity.png / .pdf")
  invisible(combined)
}
