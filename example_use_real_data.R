############################################################
# 実データを使った解析の実行例
############################################################

# データローダーを読み込み
source("data_loader.R")

############################################################
# ステップ1: サンプルデータを生成（初回のみ）
############################################################

# 3つの形式でサンプルデータを生成
create_sample_data(output_dir = "sample_data", format = "wide")
create_sample_data(output_dir = "sample_data", format = "separate")

cat("\n=== サンプルデータが生成されました ===\n\n")

############################################################
# ステップ2: データを読み込み
############################################################

# 方法1: ワイド形式（因子と時系列が1ファイル）
cat("【方法1】ワイド形式のデータ読み込み\n")
data <- load_data_wide(
  file_path = "sample_data/data_wide.csv",
  factor_cols = c("F1", "F2", "F3"),
  time_prefix = "t_"
)

factors <- data$factors
Y_matrix <- data$Y_matrix
time_grid <- data$time_grid

cat("\n因子データの先頭:\n")
print(head(factors))

cat("\n時系列データのサイズ:\n")
cat(sprintf("  %d 条件 x %d 時間点\n", nrow(Y_matrix), ncol(Y_matrix)))

# 方法2: 分離形式（因子と時系列が別ファイル）
cat("\n\n【方法2】分離形式のデータ読み込み\n")
data2 <- load_data_separate(
  factor_file = "sample_data/factors.csv",
  timeseries_file = "sample_data/timeseries.csv",
  time_prefix = "t_"
)

cat("\n読み込み完了！\n")

############################################################
# ステップ3: MFPCA.Rの解析コードを実行
############################################################

cat("\n=== FPCA解析を実行 ===\n\n")

# 以下、MFPCA.Rのコードをそのまま実行可能
library(fda)
library(ggplot2)
library(dplyr)
library(tidyr)

set.seed(123)

# データ確認
N_cond <- nrow(Y_matrix)
n_time <- ncol(Y_matrix)

cat(sprintf("データサイズ: %d条件 x %d時間点\n", N_cond, n_time))

############################################################
# B-spline基底とスムージング
############################################################

nbasis <- 25
norder <- 4
basis <- create.bspline.basis(rangeval = range(time_grid), nbasis = nbasis, norder = norder)

cat("\nスムージング実行中（GCV最適化）...\n")

Y_t <- t(Y_matrix)
lambda_candidates <- 10^seq(-6, 2, by = 0.5)

gcv_values <- numeric(length(lambda_candidates))
for (i in seq_along(lambda_candidates)) {
  fdPar_temp <- fdPar(basis, Lfdobj = 2, lambda = lambda_candidates[i])
  smooth_temp <- smooth.basis(time_grid, Y_t, fdPar_temp)
  gcv_values[i] <- mean(smooth_temp$gcv)
}

optimal_lambda <- lambda_candidates[which.min(gcv_values)]
cat(sprintf("最適λ = %.2e\n", optimal_lambda))

fdPar_obj <- fdPar(basis, Lfdobj = 2, lambda = optimal_lambda)
smooth_result <- smooth.basis(time_grid, Y_t, fdPar_obj)
func_data <- smooth_result$fd

############################################################
# FPCA実行
############################################################

M <- 4  # 主成分数

cat("\nFPCA実行中...\n")
fpca_res <- pca.fd(func_data, nharm = M, centerfns = TRUE)

eigenvalues <- fpca_res$values
variance_explained <- eigenvalues / sum(eigenvalues)
cumulative_variance <- cumsum(variance_explained)

cat("\n分散説明率:\n")
for (i in 1:M) {
  cat(sprintf("  PC%d: %.2f%% (累積: %.2f%%)\n",
              i, variance_explained[i] * 100, cumulative_variance[i] * 100))
}

scores <- fpca_res$scores
colnames(scores) <- paste0("PC", 1:M)

############################################################
# 重回帰分析
############################################################

data_reg <- cbind(factors, scores)

# 因子列を取得（condを除く）
factor_cols <- setdiff(colnames(factors), "cond")

cat("\n重回帰分析実行中...\n")

# 標準化
data_scaled <- data_reg
for (col in c(factor_cols, paste0("PC", 1:M))) {
  data_scaled[[col]] <- scale(data_scaled[[col]])
}

# 各主成分に対する回帰
regression_results <- list()
for (i in 1:M) {
  pc_name <- paste0("PC", i)
  formula_str <- paste(pc_name, "~", paste(factor_cols, collapse = " + "))
  model <- lm(as.formula(formula_str), data = data_scaled)
  regression_results[[pc_name]] <- model
}

# 係数を抽出
coef_matrix <- matrix(NA, nrow = length(factor_cols), ncol = M)
rownames(coef_matrix) <- factor_cols
colnames(coef_matrix) <- paste0("PC", 1:M)

for (i in 1:M) {
  pc_name <- paste0("PC", i)
  coefs <- coef(regression_results[[pc_name]])
  coef_matrix[, i] <- coefs[factor_cols]
}

cat("\n標準化回帰係数:\n")
print(round(coef_matrix, 3))

############################################################
# 可視化
############################################################

cat("\nヒートマップ生成中...\n")

coef_df <- as.data.frame(coef_matrix) %>%
  mutate(Factor = rownames(coef_matrix)) %>%
  pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "Beta")

p_heatmap <- ggplot(coef_df, aes(x = PC, y = Factor, fill = Beta)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f", Beta)), size = 4, fontface = "bold") +
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white",
    midpoint = 0, name = "標準化係数"
  ) +
  labs(
    title = "因子の寄与度ヒートマップ（実データ解析）",
    x = "主成分",
    y = "因子"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p_heatmap)

cat("\n=== 解析完了！ ===\n")
cat("\n次のステップ:\n")
cat("  1. sample_data/data_wide.csv を自分のデータに置き換える\n")
cat("  2. factor_cols を実際の因子名に変更する\n")
cat("  3. このスクリプトを再実行する\n")
cat("\n詳細は README_data_format.md を参照してください\n")
