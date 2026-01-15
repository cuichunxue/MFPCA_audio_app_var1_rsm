############################################################
# 簡単な確認テスト：因子調整の動作確認
############################################################

# 簡易版（依存パッケージなし）
set.seed(123)

# パラメータ
n_time <- 100
time_grid <- seq(0, 1, length.out = n_time)

# 波形生成の仕組み
base_wave <- sin(2 * pi * time_grid) + 0.3 * sin(6 * pi * time_grid)
shape_F1 <- sin(4 * pi * time_grid)
shape_F2 <- (time_grid - 0.5)
shape_F3 <- (1 - (time_grid - 0.5)^2)

beta_F1 <- 0.5
beta_F2 <- 0.4
beta_F3 <- 0.3

generate_wave <- function(F1, F2, F3) {
  base_wave + F1 * beta_F1 * shape_F1 +
    F2 * beta_F2 * shape_F2 +
    F3 * beta_F3 * shape_F3
}

# テストケース
cat("==============================================\n")
cat("因子調整テスト\n")
cat("==============================================\n\n")

# 1. 目標波形を設定
target_F1 <- 0.7
target_F2 <- -0.5
target_F3 <- 0.3
target_wave <- generate_wave(target_F1, target_F2, target_F3)

cat("【目標】\n")
cat(sprintf("  因子: F1=%.2f, F2=%.2f, F3=%.2f\n", target_F1, target_F2, target_F3))
cat(sprintf("  波形の平均値: %.4f\n", mean(target_wave)))
cat(sprintf("  波形の標準偏差: %.4f\n\n", sd(target_wave)))

# 2. 初期状態（異なる因子）
initial_F1 <- -0.2
initial_F2 <- 0.6
initial_F3 <- -0.4
initial_wave <- generate_wave(initial_F1, initial_F2, initial_F3)

initial_rmse <- sqrt(mean((initial_wave - target_wave)^2))

cat("【初期状態】\n")
cat(sprintf("  因子: F1=%.2f, F2=%.2f, F3=%.2f\n", initial_F1, initial_F2, initial_F3))
cat(sprintf("  RMSE（目標との差）: %.6f\n\n", initial_rmse))

# 3. 最適化
cat("【最適化実行】\n")

objective <- function(x) {
  wave <- generate_wave(x[1], x[2], x[3])
  mean((wave - target_wave)^2)
}

result <- optim(
  par = c(initial_F1, initial_F2, initial_F3),
  fn = objective,
  method = "L-BFGS-B",
  lower = c(-1, -1, -1),
  upper = c(1, 1, 1)
)

opt_wave <- generate_wave(result$par[1], result$par[2], result$par[3])
opt_rmse <- sqrt(mean((opt_wave - target_wave)^2))

cat(sprintf("  最適化因子: F1=%.2f, F2=%.2f, F3=%.2f\n",
            result$par[1], result$par[2], result$par[3]))
cat(sprintf("  RMSE: %.6f\n\n", opt_rmse))

# 4. 結果の評価
cat("【結果】\n")
cat("  因子の推定精度:\n")
cat(sprintf("    F1: 目標=%.2f, 推定=%.2f, 誤差=%.3f (%.1f%%)\n",
            target_F1, result$par[1],
            abs(target_F1 - result$par[1]),
            abs(target_F1 - result$par[1]) / abs(target_F1) * 100))
cat(sprintf("    F2: 目標=%.2f, 推定=%.2f, 誤差=%.3f (%.1f%%)\n",
            target_F2, result$par[2],
            abs(target_F2 - result$par[2]),
            abs(target_F2 - result$par[2]) / abs(target_F2) * 100))
cat(sprintf("    F3: 目標=%.2f, 推定=%.2f, 誤差=%.3f (%.1f%%)\n",
            target_F3, result$par[3],
            abs(target_F3 - result$par[3]),
            abs(target_F3 - result$par[3]) / abs(target_F3) * 100))

cat(sprintf("\n  RMSE改善: %.6f → %.6f (%.1f%% 改善)\n",
            initial_rmse, opt_rmse,
            (initial_rmse - opt_rmse) / initial_rmse * 100))

cat("\n==============================================\n")
cat("✓ 因子調整により目標波形に合わせることができました！\n")
cat("==============================================\n")

# 数値データのサンプルを表示
cat("\n【波形データのサンプル（最初の10点）】\n")
comparison <- data.frame(
  時間 = time_grid[1:10],
  目標 = round(target_wave[1:10], 4),
  初期 = round(initial_wave[1:10], 4),
  最適化 = round(opt_wave[1:10], 4)
)
print(comparison)
