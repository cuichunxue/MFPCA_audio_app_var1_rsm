############################################################
# ダミーデータ生成 ＋ FPCA（一変量FPCA）＋重回帰分析
# ＋ 目標波形への因子条件探索
# 条件ごとに一次元の時系列データに対する「因子の寄与」を可視化し、
# 目標波形に近づくための因子条件を最適化する
############################################################

## 必要パッケージ -----------------------------------------
# install.packages("funData")
# install.packages("MFPCA")
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("tidyr")

library(funData)
library(MFPCA)
library(ggplot2)
library(dplyr)
library(tidyr)

set.seed(123)

############################################################
# 1. ダミーデータの設計
#   - 条件：N_cond 個
#   - 時間軸：0〜1を等間隔
#   - 因子：F1, F2, F3（連続因子）を用意
#   - 因子によって波形の形状が変わるようにシミュレーション
#   - 各条件で一次元の時系列データを生成
############################################################

N_cond   <- 60     # 実験条件数（=サンプル数）
n_time   <- 100    # 波形のサンプル数（時間点数）
time_grid <- seq(0, 1, length.out = n_time)

# 因子（ここでは3個）をダミーで作成
factors <- data.frame(
  cond = 1:N_cond,
  F1 = runif(N_cond, -1, 1),  # 例：「温度」のような連続因子
  F2 = runif(N_cond, -1, 1),  # 例：「圧力」
  F3 = runif(N_cond, -1, 1)   # 例：「流量」
)

# ベース波形（条件によらず共通）
base_wave <- sin(2 * pi * time_grid) + 0.3 * sin(6 * pi * time_grid)

# 因子が掛かる「形状（関数）」を定義
shape_F1 <- sin(4 * pi * time_grid)             # 高周波の変動
shape_F2 <- (time_grid - 0.5)                   # 緩やかな傾き
shape_F3 <- (1 - (time_grid - 0.5)^2)           # 中央付近のふくらみ

# 因子の真の寄与度（重み）
beta_F1 <- 0.5
beta_F2 <- 0.4
beta_F3 <- 0.3

cat("真の因子寄与（Beta係数）:\n")
cat(sprintf("  F1: %.2f\n", beta_F1))
cat(sprintf("  F2: %.2f\n", beta_F2))
cat(sprintf("  F3: %.2f\n", beta_F3))

# 波形データ：行列 [条件 × 時間]
Y_matrix <- matrix(NA_real_, nrow = N_cond, ncol = n_time)

for (i in 1:N_cond) {
  # 条件 i の真の関数
  signal_i <- base_wave +
    factors$F1[i] * beta_F1 * shape_F1 +
    factors$F2[i] * beta_F2 * shape_F2 +
    factors$F3[i] * beta_F3 * shape_F3

  # ノイズを加える
  noise_i <- rnorm(n_time, mean = 0, sd = 0.1)

  Y_matrix[i, ] <- signal_i + noise_i
}

############################################################
# 2. funData へ変換
#    条件ごとに一次元の関数データオブジェクトを作成
############################################################

# funData を作成（一次元の関数データ）
func_data <- funData(
  argvals = time_grid,
  X       = Y_matrix  # N_cond 行 × n_time 列
)

# 確認
cat("\n一次元関数データ構造:\n")
print(func_data)

# 例：1つ目の条件の波形を確認（コメントアウト解除で使用）
# plot(func_data, obs = 1)

############################################################
# 3. FPCA（一変量FPCA）を実行
#    - 主成分数 M = 4（必要に応じて増やす）
############################################################

M <- 4  # 取り出す主成分の数

cat("\nFPCA実行中...\n")
fpca_res <- UFPCA(
  type = "uFPCA",
  funDataObject = func_data,
  npc = M
)

# 結果の概要
cat("\nFPCA結果サマリー:\n")
print(summary(fpca_res))

# 主成分の寄与率を表示
cat("\n累積寄与率:\n")
print(cumsum(fpca_res$values) / sum(fpca_res$values))

# 主成分関数を可視化（コメントアウト解除で使用）
# plot(fpca_res)
# screeplot(fpca_res)

# 各条件ごとの FPCA スコア（N_cond x M）
scores <- fpca_res$scores
colnames(scores) <- paste0("PC", 1:M)

cat("\nFPCAスコア（先頭6行）:\n")
print(head(scores))

############################################################
# 4. FPCA スコアを目的変数とした重回帰分析
#    - 各主成分 PC1〜PCM を「波形の形状特徴」とみなす
#    - 因子 F1, F2, F3 から PCスコアを予測
#    - 標準化して寄与の大きさを比較可能にする
############################################################

# スコアと因子を結合
data_reg <- cbind(factors, scores)

# 標準化（平均0, 分散1）：寄与の大きさを比較しやすくするため
data_scaled <- data_reg %>%
  mutate(
    across(starts_with("F"), scale, .names = "z_{col}")
  ) %>%
  mutate(
    across(starts_with("PC"), scale, .names = "z_{col}")
  )

# 因子名とPC名のベクトル
factor_names <- c("z_F1", "z_F2", "z_F3")
pc_names     <- paste0("z_PC", 1:M)

# 各PCごとに重回帰を実行し、標準化係数を抽出
coef_list <- list()
r_squared_list <- list()

for (pc in pc_names) {
  # モデル式： z_PCk ~ z_F1 + z_F2 + z_F3
  form <- as.formula(
    paste(pc, "~", paste(factor_names, collapse = " + "))
  )
  
  fit <- lm(form, data = data_scaled)
  
  # 切片を除いた係数（標準化済み）
  coefs <- coef(fit)[factor_names]
  
  # R二乗値を保存
  r_squared_list[[pc]] <- summary(fit)$r.squared
  
  coef_list[[pc]] <- data.frame(
    PC     = pc,
    Factor = factor_names,
    Beta   = as.numeric(coefs),
    row.names = NULL
  )
}

coef_df <- bind_rows(coef_list)

# 表示用に名前を整える
coef_df <- coef_df %>%
  mutate(
    PC     = gsub("z_", "", PC),
    Factor = gsub("z_", "", Factor)
  )

cat("\n標準化回帰係数:\n")
print(coef_df)

# R二乗値の表示
cat("\n各主成分の回帰モデルのR二乗値:\n")
r_squared_df <- data.frame(
  PC = gsub("z_", "", names(r_squared_list)),
  R_squared = unlist(r_squared_list)
)
print(r_squared_df)

############################################################
# 5. 因子の寄与を可視化
#    (a) ヒートマップ：PC × 因子 の標準化係数
#    (b) バープロット：各因子の総合寄与
#    (c) 因子別の寄与プロファイル
############################################################

# (a) ヒートマップ（係数の大きさと符号を色で表現）
p_heat <- ggplot(coef_df, aes(x = PC, y = Factor, fill = Beta)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f", Beta)), size = 4, fontface = "bold") +
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white",
    midpoint = 0, name = "標準化係数",
    limits = c(-max(abs(coef_df$Beta)), max(abs(coef_df$Beta)))
  ) +
  labs(
    title = "FPCA スコアに対する因子の寄与（標準化回帰係数）",
    subtitle = "青：負の寄与、赤：正の寄与",
    x = "主成分 (PC)",
    y = "因子"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text = element_text(face = "bold")
  )

print(p_heat)

# (b) 因子ごとの総合寄与度（各PCでの Beta^2 の合計）
importance_df <- coef_df %>%
  group_by(Factor) %>%
  summarise(
    Importance = sum(Beta^2),
    .groups = "drop"
  ) %>%
  arrange(desc(Importance))

p_bar <- ggplot(importance_df, aes(x = reorder(Factor, Importance), y = Importance)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.3f", Importance)), 
            hjust = -0.1, size = 4, fontface = "bold") +
  coord_flip() +
  labs(
    title = "因子の総合寄与度（全主成分での Beta^2 合計）",
    subtitle = "値が大きいほど波形パターンへの影響が大きい",
    x = "因子",
    y = "寄与度（標準化係数^2 の合計）"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

print(p_bar)

# (c) 因子別の主成分への寄与プロファイル
p_profile <- ggplot(coef_df, aes(x = PC, y = Beta, group = Factor, color = Factor)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "因子別の主成分への寄与プロファイル",
    subtitle = "各因子が各主成分にどう影響するか",
    x = "主成分 (PC)",
    y = "標準化係数 (Beta)",
    color = "因子"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

print(p_profile)

############################################################
# 6. 追加の詳細分析
############################################################

# 各主成分の寄与率
cat("\n各主成分の説明する分散の割合:\n")
variance_explained <- fpca_res$values / sum(fpca_res$values)
variance_df <- data.frame(
  PC = paste0("PC", 1:M),
  Variance_Explained = variance_explained,
  Cumulative = cumsum(variance_explained)
)
print(variance_df)

# 主成分の寄与率をプロット
p_variance <- ggplot(variance_df, aes(x = PC, y = Variance_Explained)) +
  geom_col(fill = "coral", alpha = 0.8) +
  geom_line(aes(y = Cumulative, group = 1), color = "darkblue", linewidth = 1.2) +
  geom_point(aes(y = Cumulative), color = "darkblue", size = 3) +
  scale_y_continuous(
    name = "個別寄与率",
    sec.axis = sec_axis(~., name = "累積寄与率")
  ) +
  labs(
    title = "主成分の分散説明率",
    x = "主成分"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p_variance)

############################################################
# 7. まとめテーブルの出力
############################################################

cat("\n=== 分析結果サマリー ===\n")
cat("\n1. データ構造:\n")
cat(sprintf("  - 条件数: %d\n", N_cond))
cat(sprintf("  - 時間点数: %d\n", n_time))
cat(sprintf("  - 因子数: 3 (F1, F2, F3)\n"))

cat("\n2. FPCA結果:\n")
cat(sprintf("  - 主成分数: %d\n", M))
cat(sprintf("  - 累積寄与率 (PC1-PC%d): %.2f%%\n",
            M, variance_df$Cumulative[M] * 100))

cat("\n3. 因子の総合寄与度（降順）:\n")
print(importance_df)

cat("\n4. 各主成分の回帰適合度:\n")
print(r_squared_df)

############################################################
# 8. 目標波形への因子条件の逆推定・最適化
############################################################

## 8.1 目標波形の設定 ----------------------------------

# 【方法1】既存の条件から目標波形を選ぶ場合
target_cond_id <- 15  # 例：条件15の波形を目標とする
target_wave <- Y_matrix[target_cond_id, ]  # [時間]

cat("\n目標条件の因子値:\n")
print(factors[target_cond_id, ])

# 【方法2】理想的な波形を新規に定義する場合（コメントアウト解除で使用）
# target_wave <- sin(3 * pi * time_grid) + 0.5 * cos(5 * pi * time_grid)

## 8.2 目標波形のFPCAスコア算出 ----------------------

# 目標波形をfunData形式に変換
target_func_data <- funData(
  argvals = time_grid,
  X = matrix(target_wave, nrow = 1)  # 1行の行列として
)

# FPCAモデルを使って目標波形のスコアを計算
# （既存のFPCA固有関数への射影）
target_scores <- predict(fpca_res, newdata = target_func_data, scores = TRUE)

cat("\n目標波形のFPCAスコア:\n")
print(target_scores)

## 8.3 逆推定：スコアから因子条件を予測 ---------------

# 逆回帰モデル：因子 ~ FPCAスコア
# 各因子ごとに逆モデルを構築

inverse_models <- list()
factor_cols <- c("F1", "F2", "F3")

for (factor_name in factor_cols) {
  # 因子を標準化
  y_factor <- scale(data_reg[[factor_name]])
  
  # 説明変数：標準化されたPCスコア
  X_scores <- as.data.frame(data_scaled[, pc_names])
  
  # 逆回帰モデル
  inverse_formula <- as.formula(
    paste("y_factor ~", paste(pc_names, collapse = " + "))
  )
  
  inverse_fit <- lm(inverse_formula, data = cbind(y_factor = y_factor, X_scores))
  
  inverse_models[[factor_name]] <- list(
    model = inverse_fit,
    mean = attr(y_factor, "scaled:center"),
    sd = attr(y_factor, "scaled:scale")
  )
  
  cat(sprintf("\n%s の逆回帰モデル R^2: %.3f\n", 
              factor_name, summary(inverse_fit)$r.squared))
}

# 目標スコアから因子条件を予測
predicted_factors <- data.frame(cond = "predicted")

for (factor_name in factor_cols) {
  model_info <- inverse_models[[factor_name]]
  
  # 標準化されたスコアで予測
  target_scores_scaled <- scale(target_scores, 
                                 center = colMeans(scores),
                                 scale = apply(scores, 2, sd))
  
  newdata <- as.data.frame(t(target_scores_scaled))
  colnames(newdata) <- pc_names
  
  # 標準化スケールでの予測
  pred_std <- predict(model_info$model, newdata = newdata)
  
  # 元のスケールに戻す
  pred_original <- pred_std * model_info$sd + model_info$mean
  
  predicted_factors[[factor_name]] <- pred_original
}

cat("\n目標波形に近づくための推定因子条件:\n")
print(predicted_factors[, factor_cols])

cat("\n実際の目標条件（参考）:\n")
print(factors[target_cond_id, factor_cols])

## 8.4 最適化による因子条件の探索 ---------------------

# 目的関数：目標波形との二乗誤差
objective_function <- function(factor_values, target_wave,
                                base_wave,
                                beta_F1, beta_F2, beta_F3,
                                shape_F1, shape_F2, shape_F3) {
  F1 <- factor_values[1]
  F2 <- factor_values[2]
  F3 <- factor_values[3]

  # 予測波形を生成
  predicted_wave <- base_wave +
    F1 * beta_F1 * shape_F1 +
    F2 * beta_F2 * shape_F2 +
    F3 * beta_F3 * shape_F3

  # 二乗誤差（全時間点）
  mse <- mean((predicted_wave - target_wave)^2)
  return(mse)
}

# 初期値：逆推定の結果を使用
initial_factors <- as.numeric(predicted_factors[, factor_cols])

# 最適化実行（制約あり：因子の範囲を-1〜1に制限）
cat("\n最適化実行中...\n")
opt_result <- optim(
  par = initial_factors,
  fn = objective_function,
  target_wave = target_wave,
  base_wave = base_wave,
  beta_F1 = beta_F1,
  beta_F2 = beta_F2,
  beta_F3 = beta_F3,
  shape_F1 = shape_F1,
  shape_F2 = shape_F2,
  shape_F3 = shape_F3,
  method = "L-BFGS-B",
  lower = c(-1, -1, -1),
  upper = c(1, 1, 1)
)

optimized_factors <- data.frame(
  cond = "optimized",
  F1 = opt_result$par[1],
  F2 = opt_result$par[2],
  F3 = opt_result$par[3]
)

cat("\n最適化された因子条件:\n")
print(optimized_factors)

cat(sprintf("\n目的関数値（MSE）: %.6f\n", opt_result$value))

## 8.5 結果の比較可視化 -------------------------------

# 比較用の波形を生成
comparison_waves <- list(
  target = target_wave,
  predicted = numeric(n_time),
  optimized = numeric(n_time)
)

# 逆推定された因子での波形
comparison_waves$predicted <- base_wave +
  predicted_factors$F1 * beta_F1 * shape_F1 +
  predicted_factors$F2 * beta_F2 * shape_F2 +
  predicted_factors$F3 * beta_F3 * shape_F3

# 最適化された因子での波形
comparison_waves$optimized <- base_wave +
  optimized_factors$F1 * beta_F1 * shape_F1 +
  optimized_factors$F2 * beta_F2 * shape_F2 +
  optimized_factors$F3 * beta_F3 * shape_F3

# 波形の比較プロット
cat("\n波形比較プロット生成中...\n")
df_plot <- data.frame(
  time = rep(time_grid, 3),
  value = c(comparison_waves$target,
            comparison_waves$predicted,
            comparison_waves$optimized),
  type = rep(c("Target", "Predicted (Inverse)", "Optimized"),
             each = n_time)
)

p_comparison <- ggplot(df_plot, aes(x = time, y = value, color = type, linetype = type)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("Target" = "black",
                                 "Predicted (Inverse)" = "blue",
                                 "Optimized" = "red")) +
  scale_linetype_manual(values = c("Target" = "solid",
                                    "Predicted (Inverse)" = "dashed",
                                    "Optimized" = "dotted")) +
  labs(
    title = "目標波形と推定波形の比較",
    x = "時間",
    y = "振幅",
    color = "波形タイプ",
    linetype = "波形タイプ"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )

print(p_comparison)

## 8.6 誤差評価 ----------------------------------------

# RMSE計算
rmse_predicted <- sqrt(mean((comparison_waves$target - comparison_waves$predicted)^2))
rmse_optimized <- sqrt(mean((comparison_waves$target - comparison_waves$optimized)^2))

cat("\nRMSE:\n")
cat(sprintf("  逆推定: %.6f\n", rmse_predicted))
cat(sprintf("  最適化: %.6f\n", rmse_optimized))

# RMSEの可視化
rmse_comparison <- data.frame(
  Method = c("Predicted (Inverse)", "Optimized"),
  RMSE = c(rmse_predicted, rmse_optimized)
)

p_rmse <- ggplot(rmse_comparison, aes(x = Method, y = RMSE, fill = Method)) +
  geom_col(alpha = 0.8) +
  scale_fill_manual(values = c("Predicted (Inverse)" = "steelblue",
                                "Optimized" = "coral")) +
  labs(
    title = "予測誤差（RMSE）比較",
    subtitle = "逆推定 vs 最適化",
    x = "推定方法",
    y = "RMSE"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "none"
  )

print(p_rmse)

## 8.7 因子条件の総合比較 ------------------------------

factor_comparison <- rbind(
  data.frame(Method = "Target (Actual)", factors[target_cond_id, factor_cols]),
  data.frame(Method = "Predicted (Inverse)", predicted_factors[, factor_cols]),
  data.frame(Method = "Optimized", optimized_factors[, factor_cols])
)

cat("\n因子条件の総合比較:\n")
print(factor_comparison)

# 因子条件の可視化
factor_long <- factor_comparison %>%
  pivot_longer(cols = c(F1, F2, F3), 
               names_to = "Factor", 
               values_to = "Value")

p_factors <- ggplot(factor_long, aes(x = Factor, y = Value, fill = Method)) +
  geom_col(position = "dodge", alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = c("Target (Actual)" = "black",
                                "Predicted (Inverse)" = "steelblue",
                                "Optimized" = "coral")) +
  labs(
    title = "目標波形達成のための因子条件比較",
    x = "因子",
    y = "因子値",
    fill = "推定方法"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )

print(p_factors)

############################################################
# 9. 感度分析：因子を変化させたときの波形への影響
############################################################

## 9.1 各因子の感度分析 -------------------------------

# 基準となる因子条件（最適化結果を使用）
baseline_factors <- optimized_factors[, factor_cols]

# 各因子を±30%変化させたときの影響を評価
sensitivity_range <- seq(-0.3, 0.3, length.out = 21)  # -30%〜+30%

sensitivity_results <- list()

for (factor_name in factor_cols) {
  results_df <- data.frame()

  for (delta in sensitivity_range) {
    # 因子を変化させる
    test_factors <- baseline_factors
    test_factors[[factor_name]] <- test_factors[[factor_name]] * (1 + delta)

    # 波形を生成
    test_wave <- base_wave +
      test_factors$F1 * beta_F1 * shape_F1 +
      test_factors$F2 * beta_F2 * shape_F2 +
      test_factors$F3 * beta_F3 * shape_F3

    # 目標との誤差を計算
    rmse_total <- sqrt(mean((test_wave - target_wave)^2))

    results_df <- rbind(results_df, data.frame(
      Factor = factor_name,
      Delta_Percent = delta * 100,
      RMSE = rmse_total
    ))
  }

  sensitivity_results[[factor_name]] <- results_df
}

sensitivity_df <- bind_rows(sensitivity_results)

# 感度分析の可視化
p_sensitivity <- ggplot(sensitivity_df, 
                        aes(x = Delta_Percent, y = RMSE, color = Factor)) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "因子変動に対する波形誤差の感度分析",
    subtitle = "最適条件からの変化率 vs RMSE",
    x = "因子の変化率 (%)",
    y = "RMSE（目標波形との誤差）",
    color = "因子"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

print(p_sensitivity)

## 9.2 因子の許容範囲の推定 --------------------------

# 目標RMSEの閾値（例：最適値の1.1倍まで許容）
rmse_threshold <- opt_result$value * 1.1

tolerance_df <- sensitivity_df %>%
  filter(RMSE <= rmse_threshold) %>%
  group_by(Factor) %>%
  summarise(
    Min_Delta = min(Delta_Percent),
    Max_Delta = max(Delta_Percent),
    Tolerance_Range = Max_Delta - Min_Delta,
    .groups = "drop"
  )

cat("\n因子の許容変動範囲（RMSE閾値以内）:\n")
print(tolerance_df)

############################################################
# 10. 実験提案：推奨される因子条件
############################################################

cat("\n" , paste(rep("=", 60), collapse = ""), "\n")
cat("【実験提案】目標波形達成のための推奨因子条件\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat("1. 最適因子条件:\n")
for (fn in factor_cols) {
  cat(sprintf("   %s = %.4f\n", fn, optimized_factors[[fn]]))
}

cat(sprintf("\n2. 予測精度: RMSE = %.6f\n", opt_result$value))

cat("\n3. 因子の許容範囲:\n")
for (i in 1:nrow(tolerance_df)) {
  cat(sprintf("   %s: %.1f%% 〜 +%.1f%% (範囲: %.1f%%)\n",
              tolerance_df$Factor[i],
              tolerance_df$Min_Delta[i],
              tolerance_df$Max_Delta[i],
              tolerance_df$Tolerance_Range[i]))
}

cat("\n4. 推奨実験順序（感度の高い順）:\n")
sensitivity_order <- tolerance_df %>%
  arrange(Tolerance_Range) %>%
  pull(Factor)

for (i in 1:length(sensitivity_order)) {
  cat(sprintf("   %d. %s（高精度制御が必要）\n", i, sensitivity_order[i]))
}

cat("\n" , paste(rep("=", 60), collapse = ""), "\n\n")

############################################################
# コード終了
#
# 【実データへの適用手順】
#
# 1. データ読み込み部分を実データに置き換え:
#    - factors: 実験条件表（温度・圧力など）
#    - Y_matrix: 実測波形データ（条件×時間）
#
# 2. パラメータ調整:
#    - M: 主成分数（累積寄与率80-90%を目安）
#    - target_cond_id: 目標とする条件ID
#    - optim の lower/upper: 因子の実験可能範囲
#
# 3. 結果の解釈:
#    - ヒートマップ: 因子と波形パターン(PC)の関係
#    - 最適化結果: 目標達成のための因子設定値
#    - 感度分析: 各因子の制御精度要求
#
############################################################