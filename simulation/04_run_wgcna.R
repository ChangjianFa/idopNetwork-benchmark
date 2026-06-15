# ============================================================
# WGCNA 封装函数
# ============================================================
# 输入：abund_matrix（行=物种，列=梯度位置）
# 输出：score_mat（TOM 相似性矩阵，对称，无向）
#
# 注：WGCNA 方向准确率固定约 0.5（属无向方法的结构特性）
# ============================================================

library(WGCNA)

run_wgcna <- function(abund_matrix, cfg) {
  tryCatch({
    t0 <- proc.time()

    # WGCNA 需要：行=样本（梯度位置），列=基因/物种
    datExpr <- t(abund_matrix)

    # 检查样本/基因质量
    gsg <- goodSamplesGenes(datExpr, verbose = 0)
    if (!gsg$allOK) {
      datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
    }

    if (nrow(datExpr) < 5 || ncol(datExpr) < 3) {
      stop("样本或物种数量不足，无法运行 WGCNA")
    }

    # 软阈值选择（自动或使用配置值）
    soft_power <- cfg$WGCNA_SOFT_POWER
    if (is.null(soft_power)) {
      powers <- c(1:10, seq(12, 20, by = 2))
      sft <- suppressWarnings(
        pickSoftThreshold(datExpr, powerVector = powers,
                          networkType = "signed", verbose = 0)
      )
      soft_power <- sft$powerEstimate
      if (is.na(soft_power) || is.null(soft_power)) soft_power <- 6L
    }

    # 计算邻接矩阵与拓扑重叠矩阵（TOM）
    adj <- adjacency(datExpr, power = soft_power, type = "signed")
    TOM <- TOMsimilarity(adj, TOMType = "signed", verbose = 0)
    rownames(TOM) <- colnames(TOM) <- colnames(datExpr)

    elapsed <- (proc.time() - t0)["elapsed"]
    list(method = "WGCNA", score_mat = TOM, sign_mat = NULL,
         runtime = as.numeric(elapsed), status = "ok",
         soft_power = soft_power)

  }, error = function(e) {
    list(method = "WGCNA", score_mat = NULL, sign_mat = NULL,
         runtime = NA_real_, status = paste("ERROR:", conditionMessage(e)))
  })
}
