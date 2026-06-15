# ============================================================
# 性能评估函数
# ============================================================
# B_true[i, j] 含义：j 调控 i（j → i）
# score_mat[from=调控者, to=目标] 对应 B_true[to, from]
#
# 评估指标：
#   - AUROC：边存在检测（所有方法）
#   - AUPRC：边存在检测，适合不平衡数据（阳性率约25%）
#   - dir_acc：有向方法的方向准确率（idopNetwork、GENIE3）
#   - sign_acc：边符号准确率（仅 idopNetwork）
# ============================================================

library(pROC)
library(PRROC)

# 为所有有序物种对生成真实标签
make_pair_labels <- function(B_true) {
  sp <- rownames(B_true)
  pairs <- expand.grid(from = sp, to = sp, stringsAsFactors = FALSE)
  pairs <- pairs[pairs$from != pairs$to, ]

  # B_true[to, from] != 0 → 存在边 from → to
  pairs$true_edge <- as.integer(
    mapply(function(fr, to) B_true[to, fr] != 0, pairs$from, pairs$to)
  )
  pairs$true_sign <- mapply(function(fr, to) sign(B_true[to, fr]),
                            pairs$from, pairs$to)
  return(pairs)
}

# 对单个方法的输出进行评估
evaluate_method <- function(result, B_true) {
  method <- result$method

  if (is.null(result$score_mat)) {
    return(data.frame(method = method,
                      AUROC    = NA_real_,
                      AUPRC    = NA_real_,
                      dir_acc  = NA_real_,
                      sign_acc = NA_real_,
                      runtime  = result$runtime,
                      status   = result$status,
                      stringsAsFactors = FALSE))
  }

  sp_true <- rownames(B_true)
  sm      <- result$score_mat
  sp_pred <- intersect(sp_true, rownames(sm))

  if (length(sp_pred) < 3) {
    return(data.frame(method = method,
                      AUROC    = NA_real_, AUPRC = NA_real_,
                      dir_acc  = NA_real_, sign_acc = NA_real_,
                      runtime  = result$runtime,
                      status   = "物种重叠不足",
                      stringsAsFactors = FALSE))
  }

  pairs <- make_pair_labels(B_true[sp_pred, sp_pred, drop = FALSE])
  pairs$score <- mapply(function(fr, to) sm[fr, to], pairs$from, pairs$to)

  # 移除得分为 NA 的行
  pairs <- pairs[!is.na(pairs$score), ]

  # AUROC
  roc_obj <- tryCatch(
    pROC::roc(pairs$true_edge, pairs$score, quiet = TRUE, direction = "<"),
    error = function(e) NULL
  )
  auroc <- if (!is.null(roc_obj)) as.numeric(pROC::auc(roc_obj)) else NA_real_

  # AUPRC（正类 = true_edge == 1）
  pos_scores <- pairs$score[pairs$true_edge == 1]
  neg_scores <- pairs$score[pairs$true_edge == 0]
  auprc <- tryCatch({
    pr_obj <- PRROC::pr.curve(scores.class0 = pos_scores,
                               scores.class1 = neg_scores)
    pr_obj$auc.integral
  }, error = function(e) NA_real_)

  # 方向准确率（仅有向方法：score[from,to] > score[to,from] 为预测方向正确）
  dir_acc <- NA_real_
  if (method %in% c("idopNetwork", "GENIE3")) {
    true_pos <- pairs[pairs$true_edge == 1, ]
    true_pos$score_rev <- mapply(function(fr, to) sm[to, fr],
                                 true_pos$from, true_pos$to)
    valid <- !is.na(true_pos$score) & !is.na(true_pos$score_rev)
    if (sum(valid) > 0) {
      dir_acc <- mean(true_pos$score[valid] > true_pos$score_rev[valid])
    }
  }

  # 符号准确率（仅 idopNetwork）
  sign_acc <- NA_real_
  if (method == "idopNetwork" && !is.null(result$sign_mat)) {
    sign_sm  <- result$sign_mat
    true_pos <- pairs[pairs$true_edge == 1, ]
    pred_sign <- mapply(function(fr, to) {
      if (fr %in% rownames(sign_sm) && to %in% colnames(sign_sm))
        sign_sm[fr, to] else NA_real_
    }, true_pos$from, true_pos$to)
    valid <- !is.na(pred_sign) & pred_sign != 0
    if (sum(valid) > 0) {
      sign_acc <- mean(pred_sign[valid] == true_pos$true_sign[valid])
    }
  }

  data.frame(method   = method,
             AUROC    = auroc,
             AUPRC    = auprc,
             dir_acc  = dir_acc,
             sign_acc = sign_acc,
             runtime  = result$runtime,
             status   = result$status,
             stringsAsFactors = FALSE)
}
