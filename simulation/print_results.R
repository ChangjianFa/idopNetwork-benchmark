setwd("d:/data/idopnetwork")
source("simulation/config.R")
source("simulation/07_visualize.R")

cfg <- list(OUT_DIR = "d:/data/idopnetwork/simulation/results",
            FIG_DIR = "d:/data/idopnetwork/simulation/figures")

results_df <- load_results(cfg$OUT_DIR)

cat("\n========== 原始结果 ==========\n")
print(results_df[, c("method","scenario","rep_id","AUROC","AUPRC","dir_acc","sign_acc","runtime")])

cat("\n========== 汇总统计 ==========\n")
tbl <- summarize_results(results_df)
key_metrics <- tbl[tbl$metric %in% c("AUROC","AUPRC","dir_acc","sign_acc"), ]
print(key_metrics)

write.csv(tbl, file.path(cfg$OUT_DIR, "summary_table.csv"), row.names = FALSE)
cat("\n汇总表已保存\n")
