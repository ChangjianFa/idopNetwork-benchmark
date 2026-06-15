# ============================================================
# idopNetwork 封装函数（顺序执行版）
# ============================================================
# 不使用 qdODE_parallel()，改用顺序循环调用 qdODE_all()。
# 原因：qdODE_parallel() 内 parLapply 的闭包序列化失败，
#       导致 worker 中 result$ 访问出现 "$ operator is invalid
#       for atomic vectors" 错误。单线程调用完全正常。
# ============================================================

library(idopNetwork)

run_idopnetwork <- function(abund_matrix, cfg) {
  tryCatch({
    t0 <- proc.time()

    # Step1: data_cleaning 要求 data.frame，第一列为物种 ID
    df <- data.frame(name = rownames(abund_matrix),
                     abund_matrix,
                     check.names = FALSE)
    cleaned <- data_cleaning(df, x = round(ncol(abund_matrix) * 0.5))

    if (nrow(cleaned) < 3) stop("data_cleaning 后物种数不足（< 3）")

    # Step2: 幂律曲线拟合（用单线程，避免嵌套 cluster 冲突）
    pfit <- power_equation_fit(cleaned,
                               n      = cfg$IDOP_N_INTERP,
                               trans  = log10,
                               thread = 1L)

    if (nrow(pfit$original_data) < 3) stop("power_equation_fit 后物种数不足")

    # Step3: 变量选择（所有物种）
    sp_names <- rownames(pfit$original_data)
    n_sp     <- length(sp_names)

    relationship <- lapply(seq_len(n_sp), function(i) {
      tryCatch(
        get_interaction(pfit$original_data, i),
        error = function(e) list(ind.name = sp_names[i], dep.name = NA, coefficient = 0)
      )
    })
    names(relationship) <- sp_names

    # Step4: 顺序求解每个物种的 ODE（替代 qdODE_parallel）
    ode_results <- lapply(seq_len(n_sp), function(i) {
      tryCatch(
        qdODE_all(result       = pfit,
                  relationship = relationship,
                  i            = i,
                  maxit        = cfg$IDOP_MAXIT),
        error = function(e) NA
      )
    })
    names(ode_results) <- sp_names

    # Step5: 将 edge 列表转换为 N×N 得分矩阵
    score_mat <- matrix(0.0, n_sp, n_sp, dimnames = list(sp_names, sp_names))
    sign_mat  <- matrix(0.0, n_sp, n_sp, dimnames = list(sp_names, sp_names))

    for (i in seq_len(n_sp)) {
      x <- ode_results[[i]]
      if (is.null(x) || (length(x) == 1 && is.na(x))) next
      if (is.null(x$fit)) next

      nc <- tryCatch(network_conversion(x), error = function(e) NULL)
      if (is.null(nc)) next

      # network_conversion 在物种只有1个调控者时返回 t(data.frame(...))
      # 即矩阵而非 data.frame，$访问会失败，故强制转为 data.frame
      edges <- as.data.frame(nc$edge, stringsAsFactors = FALSE)
      edges$Effect <- as.numeric(edges$Effect)
      for (k in seq_len(nrow(edges))) {
        fr <- edges$From[k]
        to <- edges$To[k]
        if (!is.na(fr) && !is.na(to) && fr %in% sp_names && to %in% sp_names) {
          score_mat[fr, to] <- abs(edges$Effect[k])
          sign_mat[fr, to]  <- sign(edges$Effect[k])
        }
      }
    }

    elapsed <- as.numeric((proc.time() - t0)["elapsed"])
    list(method    = "idopNetwork",
         score_mat = score_mat,
         sign_mat  = sign_mat,
         runtime   = elapsed,
         status    = "ok")

  }, error = function(e) {
    list(method    = "idopNetwork",
         score_mat = NULL,
         sign_mat  = NULL,
         runtime   = NA_real_,
         status    = paste("ERROR:", conditionMessage(e)))
  })
}
