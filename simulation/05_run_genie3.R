# ============================================================
# GENIE3 封装函数
# ============================================================
# 输入：abund_matrix（行=物种，列=梯度位置）
# 输出：score_mat[调控者, 目标]（有向，无符号）
#
# GENIE3 返回 weight_mat[target, regulator]（目标行，调控者列）
# 转置后得到 score_mat[调控者, 目标]，与 B_true 约定一致
# ============================================================

library(GENIE3)

run_genie3 <- function(abund_matrix, cfg) {
  tryCatch({
    t0 <- proc.time()

    # GENIE3 接受：行=基因/物种，列=样本（与 idopNetwork 相同方向）
    # 转为 matrix（GENIE3 要求数值矩阵）
    expr_mat <- as.matrix(abund_matrix)
    storage.mode(expr_mat) <- "double"

    if (nrow(expr_mat) < 3 || ncol(expr_mat) < 5) {
      stop("物种数或样本数不足，无法运行 GENIE3")
    }

    # 运行随机森林网络推断
    weight_mat <- GENIE3(expr_mat,
                         nTrees   = cfg$GENIE3_NTREES,
                         nCores   = cfg$IDOP_THREADS,
                         verbose  = FALSE)

    # weight_mat[target, regulator] → 转置为 score_mat[regulator, target]
    score_mat <- t(weight_mat)

    elapsed <- (proc.time() - t0)["elapsed"]
    list(method = "GENIE3", score_mat = score_mat, sign_mat = NULL,
         runtime = as.numeric(elapsed), status = "ok")

  }, error = function(e) {
    list(method = "GENIE3", score_mat = NULL, sign_mat = NULL,
         runtime = NA_real_, status = paste("ERROR:", conditionMessage(e)))
  })
}
