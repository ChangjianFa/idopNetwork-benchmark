# ============================================================
# 全局参数配置
# ============================================================

N_SPECIES      <- 10L       # 物种数（行）
N_POSITIONS    <- 30L       # 梯度位置数（列）
N_REPLICATES   <- 30L       # Monte Carlo 重复次数
NOISE_SD       <- 0.10      # log 尺度加性噪声标准差
SPARSITY_FRAC  <- 0.25      # 真实网络稀疏度（非对角非零比例）
EDGE_STRENGTH  <- c(-0.5, 0.5)  # 边权重均匀采样范围

IDOP_THREADS   <- 4L
IDOP_MAXIT     <- 1000L
IDOP_N_INTERP  <- 30L

WGCNA_SOFT_POWER <- NULL    # NULL = 自动选择
GENIE3_NTREES    <- 1000L

MASTER_SEED    <- 42L

# 路径（相对于项目根目录 d:/data/idopnetwork/）
PKG_DIR  <- "d:/data/idopnetwork/idopnetwork"
OUT_DIR  <- "d:/data/idopnetwork/simulation/results"
FIG_DIR  <- "d:/data/idopnetwork/simulation/figures"
