############################################################
# 因子調整デモ：目標波形に合わせる実証
############################################################

library(fda)
library(ggplot2)
library(dplyr)

set.seed(123)

# データ生成パラメータ
N_cond <- 60
n_time <- 100
time_grid <- seq(0, 1, length.out = n_time)

# 因子データ生成
factors <- data.frame(
  cond = 1:N_cond,
  F1 = runif(N_cond, -1, 1),
  F2 = runif(N_cond, -1, 1),
  F3 = runif(N_cond, -1, 1)
)

# ベース波形と因子の形状関数
base_wave <- sin(2 * pi * time_grid) + 0.3 * sin(6 * pi * time_grid)
shape_F1 <- sin(4 * pi * time_grid)
shape_F2 <- (time_grid - 0.5)
shape_F3 <- (1 - (time_grid - 0.5)^2)

# 真の寄与度
beta_F1 <- 0.5
beta_F2 <- 0.4
beta_F3 <- 0.3

cat("=== 因子調整デモ ===\n\n")
cat("真の因子寄与:\n")
cat(sprintf("  F1: %.2f, F2: %.2f, F3: %.2f\n\n", beta_F1, beta_F2, beta_F3))

# 波形生成関数
generate_wave <- function(F1, F2, F3) {
  base_wave + F1 * beta_F1 * shape_F1 +
    F2 * beta_F2 * shape_F2 +
    F3 * beta_F3 * shape_F3
}

# 目標波形を設定（例：F1=0.8, F2=-0.6, F3=0.4）
target_factors <- c(F1 = 0.8, F2 = -0.6, F3 = 0.4)
target_wave <- generate_wave(target_factors[1], target_factors[2], target_factors[3])

cat("【ステップ1】目標波形の設定\n")
cat(sprintf("  目標因子: F1=%.2f, F2=%.2f, F3=%.2f\n\n",
            target_factors[1], target_factors[2], target_factors[3]))

# 初期値（全く異なる因子値）
initial_factors <- c(F1 = -0.3, F2 = 0.5, F3 = -0.7)
initial_wave <- generate_wave(initial_factors[1], initial_factors[2], initial_factors[3])

cat("【ステップ2】初期状態（目標とは異なる因子）\n")
cat(sprintf("  初期因子: F1=%.2f, F2=%.2f, F3=%.2f\n",
            initial_factors[1], initial_factors[2], initial_factors[3]))
cat(sprintf("  初期RMSE: %.6f\n\n", sqrt(mean((initial_wave - target_wave)^2))))

# 最適化による因子調整
cat("【ステップ3】最適化実行中...\n")

objective_function <- function(factor_values) {
  wave <- generate_wave(factor_values[1], factor_values[2], factor_values[3])
  mse <- mean((wave - target_wave)^2)
  return(mse)
}

opt_result <- optim(
  par = initial_factors,
  fn = objective_function,
  method = "L-BFGS-B",
  lower = c(-1, -1, -1),
  upper = c(1, 1, 1)
)

optimized_factors <- opt_result$par
optimized_wave <- generate_wave(optimized_factors[1], optimized_factors[2], optimized_factors[3])

cat("最適化完了！\n\n")
cat("【ステップ4】最適化結果\n")
cat(sprintf("  最適化因子: F1=%.2f, F2=%.2f, F3=%.2f\n",
            optimized_factors[1], optimized_factors[2], optimized_factors[3]))
cat(sprintf("  最適化RMSE: %.6f\n\n", sqrt(mean((optimized_wave - target_wave)^2))))

# 精度評価
cat("【ステップ5】因子推定の精度\n")
cat("  目標因子  vs  最適化因子:\n")
for (i in 1:3) {
  factor_name <- c("F1", "F2", "F3")[i]
  error <- abs(target_factors[i] - optimized_factors[i])
  error_pct <- error / abs(target_factors[i]) * 100
  cat(sprintf("    %s: %.3f  vs  %.3f  (誤差: %.3f, %.1f%%)\n",
              factor_name, target_factors[i], optimized_factors[i], error, error_pct))
}

# RMSE改善率
initial_rmse <- sqrt(mean((initial_wave - target_wave)^2))
final_rmse <- sqrt(mean((optimized_wave - target_wave)^2))
improvement <- (initial_rmse - final_rmse) / initial_rmse * 100

cat(sprintf("\nRMSE改善率: %.2f%%\n", improvement))

# 可視化
cat("\n【ステップ6】波形の可視化\n")

plot_data <- data.frame(
  time = rep(time_grid, 3),
  amplitude = c(target_wave, initial_wave, optimized_wave),
  type = factor(rep(c("目標波形", "初期波形", "最適化波形"), each = n_time),
                levels = c("目標波形", "初期波形", "最適化波形"))
)

p <- ggplot(plot_data, aes(x = time, y = amplitude, color = type, linetype = type)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("目標波形" = "black",
                                 "初期波形" = "red",
                                 "最適化波形" = "blue")) +
  scale_linetype_manual(values = c("目標波形" = "solid",
                                    "初期波形" = "dashed",
                                    "最適化波形" = "dotted")) +
  labs(
    title = "因子調整による波形の最適化",
    subtitle = sprintf("初期RMSE: %.4f → 最適化RMSE: %.4f (改善率: %.1f%%)",
                      initial_rmse, final_rmse, improvement),
    x = "時間",
    y = "振幅",
    color = "波形タイプ",
    linetype = "波形タイプ"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

print(p)
ggsave("/tmp/factor_optimization_demo.png", p, width = 10, height = 6, dpi = 300)

cat("\n=== デモ完了 ===\n")
cat("結論: 因子を調整することで、目標の時系列に合わせることができました！\n")
cat(sprintf("  - 因子推定精度: 平均誤差 %.1f%%\n", mean(abs(target_factors - optimized_factors) / abs(target_factors) * 100)))
cat(sprintf("  - 波形再現精度: RMSE %.6f\n", final_rmse))
cat("\nプロット保存先: /tmp/factor_optimization_demo.png\n")
