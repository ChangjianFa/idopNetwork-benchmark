# ============================================================
# 可视化模拟生成的物种丰度曲线
# ============================================================

library(ggplot2)

plot_abundance_curves <- function(cfg) {
  lv_path   <- file.path(cfg$OUT_DIR, "sim_lv.rds")
  var1_path  <- file.path(cfg$OUT_DIR, "sim_var1.rds")

  if (!file.exists(lv_path) || !file.exists(var1_path)) {
    stop("找不到 sim_lv.rds 或 sim_var1.rds，请先运行 run_all.R")
  }

  lv   <- readRDS(lv_path)
  var1 <- readRDS(var1_path)

  # 将矩阵展开为 long format，取前 3 个重复
  make_long <- function(data, scenario_label, n_rep = 3) {
    reps <- seq_len(min(n_rep, length(data$replicates)))
    do.call(rbind, lapply(reps, function(r) {
      mat <- data$replicates[[r]]$abund
      sp  <- rownames(mat)
      n_pos <- ncol(mat)
      data.frame(
        position  = rep(seq_len(n_pos), each = nrow(mat)),
        abundance = as.vector(mat),
        species   = rep(sp, times = n_pos),
        replicate = paste0("rep", r),
        scenario  = scenario_label,
        stringsAsFactors = FALSE
      )
    }))
  }

  df_lv   <- make_long(lv,   "场景A: Lotka-Volterra ODE")
  df_var1 <- make_long(var1, "场景B: 线性 VAR(1)")
  df      <- rbind(df_lv, df_var1)

  # 过滤掉 <= 0 的值（log 轴需要正数）
  df <- df[df$abundance > 0, ]

  # ── 图1：LV 与 VAR1 并排，前 3 次重复叠加 ──────────────────────────
  p1 <- ggplot(df,
               aes(x     = position,
                   y     = abundance,
                   color = species,
                   group = interaction(species, replicate),
                   alpha = replicate)) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~scenario, scales = "free_y", ncol = 1) +
    scale_y_log10() +
    scale_alpha_manual(values = c(rep1 = 1.0, rep2 = 0.6, rep3 = 0.35),
                       guide  = "none") +
    scale_color_brewer(palette = "Set3") +
    labs(title    = "模拟物种丰度曲线（前 3 次重复，透明度递减）",
         subtitle = "实线=rep1，半透明=rep2/rep3；y 轴为 log10 尺度",
         x = "梯度位置",
         y = "丰度（log10）",
         color = "物种") +
    theme_bw(base_size = 13) +
    theme(legend.position    = "right",
          strip.background   = element_rect(fill = "#F0F0F0"),
          strip.text         = element_text(face = "bold"),
          plot.title         = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle      = element_text(hjust = 0.5, color = "gray40"))

  out1 <- file.path(cfg$FIG_DIR, "fig0_abundance_curves.png")
  ggsave(out1, p1, width = 10, height = 8, dpi = 300)
  message("已保存 ", out1)

  # ── 图2：仅 rep1，两场景并排，线性 y 轴（看绝对丰度范围）────────────
  df_rep1 <- df[df$replicate == "rep1", ]

  p2 <- ggplot(df_rep1,
               aes(x = position, y = abundance,
                   color = species, group = species)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.2, alpha = 0.7) +
    facet_wrap(~scenario, scales = "free_y", ncol = 2) +
    scale_color_brewer(palette = "Set3") +
    labs(title  = "第 1 次重复的物种丰度曲线（线性尺度）",
         x = "梯度位置",
         y = "丰度",
         color = "物种") +
    theme_bw(base_size = 13) +
    theme(legend.position  = "right",
          strip.background = element_rect(fill = "#F0F0F0"),
          strip.text       = element_text(face = "bold"),
          plot.title       = element_text(hjust = 0.5, face = "bold"))

  out2 <- file.path(cfg$FIG_DIR, "fig0b_abundance_rep1.png")
  ggsave(out2, p2, width = 12, height = 5, dpi = 300)
  message("已保存 ", out2)

  invisible(list(p_multi = p1, p_rep1 = p2))
}
