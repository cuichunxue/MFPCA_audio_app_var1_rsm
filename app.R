############################################################
# FPCA 音響最適化 インタラクティブアプリ
# - 因子調整と目標波形の差異をリアルタイム確認（一次元時系列）
# - モデル選択（重回帰分析 vs ランダムフォレスト vs ガウス過程）
# - 最適化探索機能（L-BFGS-B vs ベイズ最適化）
############################################################

library(shiny)
library(shinydashboard)
library(fda)       # 関数型データ分析
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)  # インタラクティブプロット用
library(randomForest)  # ランダムフォレスト用
library(DiceKriging)  # ガウス過程用
library(rBayesianOptimization)  # ベイズ最適化用
library(fastshap)  # SHAP値計算用

set.seed(123)

############################################################
# データ生成とFPCA実行（グローバル）
############################################################

# データ生成パラメータ
N_cond   <- 60
n_time   <- 100
time_grid <- seq(0, 1, length.out = n_time)

# 因子データ生成
factors <- data.frame(
  cond = 1:N_cond,
  F1 = runif(N_cond, -1, 1),
  F2 = runif(N_cond, -1, 1),
  F3 = runif(N_cond, -1, 1)
)

# ベース波形
base_wave <- sin(2 * pi * time_grid) + 0.3 * sin(6 * pi * time_grid)

# 因子の形状関数
shape_F1 <- sin(4 * pi * time_grid)
shape_F2 <- (time_grid - 0.5)
shape_F3 <- (1 - (time_grid - 0.5)^2)

# 因子の真の寄与度（重み）
beta_F1 <- 0.5
beta_F2 <- 0.4
beta_F3 <- 0.3

# 波形データ生成：行列 [条件 × 時間]
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

# fdaパッケージでスムージングとFPCA
nbasis <- 25
norder <- 4
basis <- create.bspline.basis(rangeval = c(0, 1), nbasis = nbasis, norder = norder)

# GCVによる最適スムージングパラメータの選択
Y_t <- t(Y_matrix)
lambda_candidates <- 10^seq(-6, 2, by = 0.5)
gcv_values <- numeric(length(lambda_candidates))
for (i in seq_along(lambda_candidates)) {
  fdPar_temp <- fdPar(basis, Lfdobj = 2, lambda = lambda_candidates[i])
  smooth_temp <- smooth.basis(time_grid, Y_t, fdPar_temp)
  gcv_values[i] <- mean(smooth_temp$gcv)
}
optimal_lambda <- lambda_candidates[which.min(gcv_values)]

# 最適λでスムージング
fdPar_obj <- fdPar(basis, Lfdobj = 2, lambda = optimal_lambda)
smooth_result <- smooth.basis(time_grid, Y_t, fdPar_obj)
func_data <- smooth_result$fd

# FPCA実行
M <- 4
fpca_res <- pca.fd(func_data, nharm = M, centerfns = TRUE)

# PCスコア
all_pc_scores <- fpca_res$scores
colnames(all_pc_scores) <- paste0("PC", 1:M)

# 標準化回帰分析
data_reg <- cbind(factors, all_pc_scores)
data_scaled <- data_reg %>%
  mutate(across(starts_with("F"), scale, .names = "z_{col}")) %>%
  mutate(across(starts_with("PC"), scale, .names = "z_{col}"))

factor_names <- c("z_F1", "z_F2", "z_F3")
pc_names <- paste0("z_PC", 1:M)

# 回帰係数計算
coef_list <- list()
for (pc in pc_names) {
  form <- as.formula(paste(pc, "~", paste(factor_names, collapse = " + ")))
  fit <- lm(form, data = data_scaled)
  coefs <- coef(fit)[factor_names]
  coef_list[[pc]] <- data.frame(
    PC = pc,
    Factor = factor_names,
    Beta = as.numeric(coefs),
    row.names = NULL
  )
}
coef_df <- bind_rows(coef_list) %>%
  mutate(PC = gsub("z_", "", PC), Factor = gsub("z_", "", Factor))

# 分散説明率
eigenvalues <- fpca_res$values[1:M]
variance_explained <- eigenvalues / sum(eigenvalues)
variance_df <- data.frame(
  PC = paste0("PC", 1:M),
  Variance_Explained = variance_explained,
  Cumulative = cumsum(variance_explained)
)

############################################################
# モデル訓練（グローバル）
############################################################

# 訓練データとテストデータに分割（80/20）
set.seed(456)
train_idx <- sample(1:N_cond, size = floor(0.8 * N_cond))
test_idx <- setdiff(1:N_cond, train_idx)

train_data <- data_reg[train_idx, ]
test_data <- data_reg[test_idx, ]

# 重回帰モデル（因子 → PCスコア）
lm_models <- list()
for (pc in paste0("PC", 1:M)) {
  formula_str <- paste(pc, "~ F1 + F2 + F3")
  lm_models[[pc]] <- lm(as.formula(formula_str), data = train_data)
}

# ランダムフォレストモデル（因子 → PCスコア）
rf_models <- list()
for (pc in paste0("PC", 1:M)) {
  formula_str <- paste(pc, "~ F1 + F2 + F3")
  rf_models[[pc]] <- randomForest(as.formula(formula_str), data = train_data,
                                   ntree = 500, importance = TRUE)
}

# ガウス過程モデル（因子 → PCスコア）
gp_models <- list()
for (pc in paste0("PC", 1:M)) {
  # DiceKrigingはデータフレームではなく行列が必要
  X_train <- as.matrix(train_data[, c("F1", "F2", "F3")])
  y_train <- train_data[[pc]]

  # ガウス過程モデルの構築（デフォルトカーネル：Matern 5/2）
  gp_models[[pc]] <- km(
    formula = ~1,  # 定数平均関数
    design = X_train,
    response = y_train,
    covtype = "matern5_2",  # Maternカーネル
    control = list(trace = FALSE)
  )
}

# モデル評価関数
evaluate_model <- function(models, data, model_type) {
  results <- list()

  for (pc in paste0("PC", 1:M)) {
    actual <- data[[pc]]

    if (model_type == "lm") {
      predicted <- predict(models[[pc]], newdata = data)
    } else if (model_type == "rf") {
      predicted <- predict(models[[pc]], newdata = data)
    } else if (model_type == "gp") {
      # ガウス過程の予測
      X_test <- as.matrix(data[, c("F1", "F2", "F3")])
      pred_result <- predict(models[[pc]], newdata = X_test, type = "UK")
      predicted <- pred_result$mean
    }

    # メトリクス計算
    rmse <- sqrt(mean((actual - predicted)^2))
    mae <- mean(abs(actual - predicted))
    r_squared <- cor(actual, predicted)^2

    results[[pc]] <- data.frame(
      PC = pc,
      RMSE = rmse,
      MAE = mae,
      R_squared = r_squared
    )
  }

  bind_rows(results)
}

# 訓練データでの評価
lm_train_perf <- evaluate_model(lm_models, train_data, "lm")
rf_train_perf <- evaluate_model(rf_models, train_data, "rf")
gp_train_perf <- evaluate_model(gp_models, train_data, "gp")

# テストデータでの評価
lm_test_perf <- evaluate_model(lm_models, test_data, "lm")
rf_test_perf <- evaluate_model(rf_models, test_data, "rf")
gp_test_perf <- evaluate_model(gp_models, test_data, "gp")

############################################################
# ヘルパー関数
############################################################

# 因子から波形を生成
generate_waveform <- function(F1, F2, F3) {
  wave <- base_wave +
    F1 * beta_F1 * shape_F1 +
    F2 * beta_F2 * shape_F2 +
    F3 * beta_F3 * shape_F3
  return(wave)
}

# 波形からPCスコアを計算（最適化で使用）
calculate_pc_scores <- function(wave) {
  # 波形をfdオブジェクトに変換（同じ基底とスムージングパラメータを使用）
  target_smooth <- smooth.basis(time_grid, wave, fdPar_obj)
  target_fd <- target_smooth$fd

  # FPCAスコアを計算（既存の主成分への射影）
  # 目標波形を中心化
  target_centered <- target_fd - fpca_res$meanfd

  # 各主成分との内積を計算してスコアを求める
  target_scores <- numeric(M)
  for (k in 1:M) {
    target_scores[k] <- inprod(target_centered, fpca_res$harmonics[k])
  }

  return(target_scores)
}

# 差異メトリクス計算
calculate_metrics <- function(wave1, wave2) {
  mse <- mean((wave1 - wave2)^2)
  rmse <- sqrt(mse)

  # 相関係数
  correlation <- cor(wave1, wave2)

  list(
    mse = mse,
    rmse = rmse,
    correlation = correlation
  )
}

# 最適化関数（モデルベース）
optimize_factors_model <- function(target_pc_scores, model_type, models) {
  objective_function <- function(factor_values) {
    F1 <- factor_values[1]
    F2 <- factor_values[2]
    F3 <- factor_values[3]

    # 因子からPCスコアを予測
    newdata <- data.frame(F1 = F1, F2 = F2, F3 = F3)

    predicted_scores <- numeric(M)
    for (i in 1:M) {
      pc_name <- paste0("PC", i)
      if (model_type == "lm") {
        predicted_scores[i] <- predict(models[[pc_name]], newdata = newdata)
      } else if (model_type == "rf") {
        predicted_scores[i] <- predict(models[[pc_name]], newdata = newdata)
      } else if (model_type == "gp") {
        # ガウス過程の予測
        X_new <- matrix(c(F1, F2, F3), nrow = 1)
        pred_result <- predict(models[[pc_name]], newdata = X_new, type = "UK")
        predicted_scores[i] <- pred_result$mean
      }
    }

    # 目標PCスコアとの二乗誤差
    mse <- mean((predicted_scores - target_pc_scores)^2)
    return(mse)
  }

  # 初期値（ランダム）
  initial_factors <- c(0, 0, 0)

  opt_result <- optim(
    par = initial_factors,
    fn = objective_function,
    method = "L-BFGS-B",
    lower = c(-1, -1, -1),
    upper = c(1, 1, 1)
  )

  return(list(
    factors = opt_result$par,
    mse = opt_result$value,
    convergence = opt_result$convergence
  ))
}

# 最適化関数（直接波形ベース）
optimize_factors_direct <- function(target_wave) {
  objective_function <- function(factor_values) {
    F1 <- factor_values[1]
    F2 <- factor_values[2]
    F3 <- factor_values[3]

    predicted_wave <- generate_waveform(F1, F2, F3)
    mse <- mean((predicted_wave - target_wave)^2)
    return(mse)
  }

  # 初期値（ランダム）
  initial_factors <- c(0, 0, 0)

  opt_result <- optim(
    par = initial_factors,
    fn = objective_function,
    method = "L-BFGS-B",
    lower = c(-1, -1, -1),
    upper = c(1, 1, 1)
  )

  return(list(
    factors = opt_result$par,
    mse = opt_result$value,
    convergence = opt_result$convergence
  ))
}

# ベイズ最適化関数（直接波形ベース）
optimize_factors_bayesian_direct <- function(target_wave) {
  # ベイズ最適化では目的関数を最大化するため、負のMSEを返す
  objective_function <- function(F1, F2, F3) {
    predicted_wave <- generate_waveform(F1, F2, F3)
    mse <- mean((predicted_wave - target_wave)^2)
    return(list(Score = -mse))  # 負のMSEを返して最大化
  }

  # ベイズ最適化実行
  opt_result <- BayesianOptimization(
    FUN = objective_function,
    bounds = list(F1 = c(-1, 1), F2 = c(-1, 1), F3 = c(-1, 1)),
    init_points = 5,  # 初期ランダムサンプル数
    n_iter = 20,  # 最適化イテレーション数
    acq = "ucb",  # 獲得関数：Upper Confidence Bound
    kappa = 2.576,  # UCBパラメータ
    verbose = FALSE
  )

  return(list(
    factors = c(opt_result$Best_Par["F1"],
                opt_result$Best_Par["F2"],
                opt_result$Best_Par["F3"]),
    mse = -opt_result$Best_Value,  # 負を元に戻す
    convergence = 0  # ベイズ最適化は常に収束
  ))
}

# ベイズ最適化関数（モデルベース）
optimize_factors_bayesian_model <- function(target_pc_scores, model_type, models) {
  # ベイズ最適化では目的関数を最大化するため、負のMSEを返す
  objective_function <- function(F1, F2, F3) {
    newdata <- data.frame(F1 = F1, F2 = F2, F3 = F3)

    predicted_scores <- numeric(M)
    for (i in 1:M) {
      pc_name <- paste0("PC", i)
      if (model_type == "lm") {
        predicted_scores[i] <- predict(models[[pc_name]], newdata = newdata)
      } else if (model_type == "rf") {
        predicted_scores[i] <- predict(models[[pc_name]], newdata = newdata)
      } else if (model_type == "gp") {
        X_new <- matrix(c(F1, F2, F3), nrow = 1)
        pred_result <- predict(models[[pc_name]], newdata = X_new, type = "UK")
        predicted_scores[i] <- pred_result$mean
      }
    }

    mse <- mean((predicted_scores - target_pc_scores)^2)
    return(list(Score = -mse))  # 負のMSEを返して最大化
  }

  # ベイズ最適化実行
  opt_result <- BayesianOptimization(
    FUN = objective_function,
    bounds = list(F1 = c(-1, 1), F2 = c(-1, 1), F3 = c(-1, 1)),
    init_points = 5,  # 初期ランダムサンプル数
    n_iter = 20,  # 最適化イテレーション数
    acq = "ucb",  # 獲得関数：Upper Confidence Bound
    kappa = 2.576,  # UCBパラメータ
    verbose = FALSE
  )

  return(list(
    factors = c(opt_result$Best_Par["F1"],
                opt_result$Best_Par["F2"],
                opt_result$Best_Par["F3"]),
    mse = -opt_result$Best_Value,  # 負を元に戻す
    convergence = 0  # ベイズ最適化は常に収束
  ))
}

############################################################
# UI定義
############################################################

ui <- dashboardPage(
  dashboardHeader(title = "FPCA 音響最適化システム"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("概要とFPCA結果", tabName = "overview", icon = icon("chart-line")),
      menuItem("モデル選択と評価", tabName = "model_selection", icon = icon("brain")),
      menuItem("因子調整と最適化", tabName = "factor_adjust", icon = icon("sliders-h")),
      menuItem("感度分析", tabName = "sensitivity", icon = icon("chart-area"))
    )
  ),

  dashboardBody(
    tabItems(
      # タブ1: 概要とFPCA結果
      tabItem(
        tabName = "overview",
        fluidRow(
          box(
            title = "FPCA分散説明率", width = 6, status = "primary",
            plotOutput("variance_plot", height = 300)
          ),
          box(
            title = "因子寄与ヒートマップ", width = 6, status = "primary",
            plotOutput("coef_heatmap", height = 300)
          )
        ),
        fluidRow(
          box(
            title = "因子総合寄与度", width = 6, status = "info",
            plotOutput("importance_bar", height = 300)
          ),
          box(
            title = "データ情報", width = 6, status = "info",
            verbatimTextOutput("data_info")
          )
        )
      ),

      # タブ2: モデル選択と評価
      tabItem(
        tabName = "model_selection",
        fluidRow(
          box(
            title = "モデル選択", width = 4, status = "primary",
            radioButtons("model_type", "予測モデルを選択:",
                        choices = list(
                          "重回帰分析（線形）" = "lm",
                          "ランダムフォレスト（非線形）" = "rf",
                          "ガウス過程（確率的非線形）" = "gp"
                        ),
                        selected = "lm"),
            hr(),
            h4("モデルの特徴:"),
            conditionalPanel(
              condition = "input.model_type == 'lm'",
              tags$ul(
                tags$li("線形関係を仮定"),
                tags$li("解釈性が高い"),
                tags$li("計算が高速"),
                tags$li("外挿に注意")
              )
            ),
            conditionalPanel(
              condition = "input.model_type == 'rf'",
              tags$ul(
                tags$li("非線形関係を捉える"),
                tags$li("ロバスト性が高い"),
                tags$li("過学習のリスク低"),
                tags$li("変数重要度を提供")
              )
            ),
            conditionalPanel(
              condition = "input.model_type == 'gp'",
              tags$ul(
                tags$li("非線形関係を柔軟にモデル化"),
                tags$li("予測の不確実性を定量化"),
                tags$li("少ないデータで高精度"),
                tags$li("滑らかな予測曲面")
              )
            )
          ),
          box(
            title = "モデル性能比較", width = 8, status = "primary",
            tabsetPanel(
              tabPanel("テストデータ",
                      h4("テストデータでの予測精度"),
                      tableOutput("model_comparison_test")),
              tabPanel("訓練データ",
                      h4("訓練データでの予測精度"),
                      tableOutput("model_comparison_train"))
            )
          )
        ),
        fluidRow(
          box(
            title = "予測精度の可視化", width = 6, status = "info",
            plotOutput("model_performance_plot", height = 400)
          ),
          box(
            title = "変数重要度", width = 6, status = "info",
            plotOutput("variable_importance_plot", height = 400)
          )
        )
      ),

      # タブ3: 因子調整と最適化（統合版）
      tabItem(
        tabName = "factor_adjust",
        fluidRow(
          box(
            title = "因子調整", width = 4, status = "warning",
            sliderInput("f1", "F1 (因子1)",
                        min = -1, max = 1, value = 0, step = 0.01),
            sliderInput("f2", "F2 (因子2)",
                        min = -1, max = 1, value = 0, step = 0.01),
            sliderInput("f3", "F3 (因子3)",
                        min = -1, max = 1, value = 0, step = 0.01),
            hr(),
            selectInput("target_cond", "目標条件を選択",
                        choices = 1:N_cond, selected = 15),
            actionButton("reset_factors", "因子をリセット",
                        class = "btn-warning btn-block", icon = icon("undo"))
          ),
          box(
            title = "差異メトリクス", width = 8, status = "warning",
            fluidRow(
              valueBoxOutput("mse_box", width = 4),
              valueBoxOutput("rmse_box", width = 4),
              valueBoxOutput("cor_box", width = 4)
            )
          )
        ),
        fluidRow(
          box(
            title = "最適値探索", width = 4, status = "danger",
            radioButtons("opt_algorithm", "最適化アルゴリズム:",
                        choices = list(
                          "L-BFGS-B（勾配ベース）" = "lbfgs",
                          "ベイズ最適化（サンプルベース）" = "bayes"
                        ),
                        selected = "lbfgs"),
            hr(),
            radioButtons("opt_method", "目的関数:",
                        choices = list(
                          "直接法（波形ベース）" = "direct",
                          "モデルベース（PCスコア）" = "model"
                        ),
                        selected = "direct"),
            conditionalPanel(
              condition = "input.opt_method == 'model'",
              p(strong("使用モデル:")),
              textOutput("selected_model_text")
            ),
            conditionalPanel(
              condition = "input.opt_algorithm == 'bayes'",
              p(strong("ベイズ最適化設定:")),
              tags$ul(
                tags$li("初期サンプル: 5点"),
                tags$li("最適化イテレーション: 20回"),
                tags$li("獲得関数: UCB")
              )
            ),
            hr(),
            actionButton("run_optimization", "最適値探索を実行",
                        class = "btn-danger btn-lg btn-block",
                        icon = icon("rocket")),
            hr(),
            h4("最適化結果:"),
            verbatimTextOutput("opt_result"),
            hr(),
            actionButton("apply_opt_factors", "最適化結果を適用",
                        class = "btn-success btn-block",
                        icon = icon("check"))
          ),
          box(
            title = "波形比較", width = 8, status = "primary", solidHeader = TRUE,
            plotOutput("waveform_comparison", height = 600)
          )
        )
      ),

      # タブ4: 感度分析
      tabItem(
        tabName = "sensitivity",
        fluidRow(
          box(
            title = "因子感度分析", width = 12, status = "info", solidHeader = TRUE,
            plotOutput("sensitivity_plot", height = 500)
          )
        ),
        fluidRow(
          box(
            title = "因子許容範囲", width = 12, status = "info",
            tableOutput("tolerance_table")
          )
        )
      )
    )
  )
)

############################################################
# Server定義
############################################################

server <- function(input, output, session) {

  # リアクティブ値
  rv <- reactiveValues(
    opt_factors = NULL,
    opt_mse = NULL,
    opt_method_used = NULL
  )

  # 目標波形（リアクティブ）
  target_wave <- reactive({
    target_cond <- as.integer(input$target_cond)
    Y_matrix[target_cond, ]
  })

  # 目標のPCスコア（最適化で使用）
  target_pc <- reactive({
    calculate_pc_scores(target_wave())
  })

  # 現在の因子から波形を生成
  current_wave <- reactive({
    generate_waveform(input$f1, input$f2, input$f3)
  })

  # 差異メトリクス
  metrics <- reactive({
    calculate_metrics(current_wave(), target_wave())
  })

  # タブ1: 概要 ---------------------------------------------

  output$variance_plot <- renderPlot({
    ggplot(variance_df, aes(x = PC, y = Variance_Explained)) +
      geom_col(fill = "coral", alpha = 0.8) +
      geom_line(aes(y = Cumulative, group = 1), color = "darkblue", linewidth = 1.5) +
      geom_point(aes(y = Cumulative), color = "darkblue", size = 4) +
      geom_text(aes(label = sprintf("%.1f%%", Variance_Explained * 100)),
                vjust = -0.5, size = 4) +
      scale_y_continuous(
        name = "個別寄与率",
        sec.axis = sec_axis(~., name = "累積寄与率"),
        labels = scales::percent
      ) +
      labs(title = "主成分の分散説明率", x = "主成分") +
      theme_minimal(base_size = 14) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  })

  output$coef_heatmap <- renderPlot({
    ggplot(coef_df, aes(x = PC, y = Factor, fill = Beta)) +
      geom_tile(color = "white", linewidth = 1) +
      geom_text(aes(label = sprintf("%.2f", Beta)), size = 5, fontface = "bold") +
      scale_fill_gradient2(
        low = "blue", high = "red", mid = "white",
        midpoint = 0, name = "標準化係数",
        limits = c(-max(abs(coef_df$Beta)), max(abs(coef_df$Beta)))
      ) +
      labs(
        title = "因子の主成分への寄与",
        x = "主成分 (PC)", y = "因子"
      ) +
      theme_minimal(base_size = 14) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  })

  output$importance_bar <- renderPlot({
    importance_df <- coef_df %>%
      group_by(Factor) %>%
      summarise(Importance = sum(Beta^2), .groups = "drop") %>%
      arrange(desc(Importance))

    ggplot(importance_df, aes(x = reorder(Factor, Importance), y = Importance)) +
      geom_col(fill = "steelblue", alpha = 0.8) +
      geom_text(aes(label = sprintf("%.3f", Importance)),
                hjust = -0.1, size = 5, fontface = "bold") +
      coord_flip() +
      labs(
        title = "因子の総合寄与度",
        x = "因子", y = "寄与度（Beta² 合計）"
      ) +
      theme_minimal(base_size = 14) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  })

  output$data_info <- renderPrint({
    cat("=== データ構造 ===\n\n")
    cat(sprintf("条件数: %d\n", N_cond))
    cat(sprintf("時間点数: %d\n", n_time))
    cat(sprintf("因子数: 3 (F1, F2, F3)\n"))
    cat(sprintf("\n主成分数: %d\n", M))
    cat(sprintf("累積寄与率: %.2f%%\n", variance_df$Cumulative[M] * 100))
  })

  # タブ2: モデル選択と評価 ---------------------------------

  output$model_comparison_test <- renderTable({
    lm_perf <- lm_test_perf %>% mutate(Model = "重回帰分析")
    rf_perf <- rf_test_perf %>% mutate(Model = "ランダムフォレスト")
    gp_perf <- gp_test_perf %>% mutate(Model = "ガウス過程")

    comparison <- bind_rows(lm_perf, rf_perf, gp_perf) %>%
      select(Model, PC, R_squared, RMSE, MAE)

    comparison
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  output$model_comparison_train <- renderTable({
    lm_perf <- lm_train_perf %>% mutate(Model = "重回帰分析")
    rf_perf <- rf_train_perf %>% mutate(Model = "ランダムフォレスト")
    gp_perf <- gp_train_perf %>% mutate(Model = "ガウス過程")

    comparison <- bind_rows(lm_perf, rf_perf, gp_perf) %>%
      select(Model, PC, R_squared, RMSE, MAE)

    comparison
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  output$model_performance_plot <- renderPlot({
    lm_perf <- lm_test_perf %>% mutate(Model = "重回帰分析")
    rf_perf <- rf_test_perf %>% mutate(Model = "ランダムフォレスト")
    gp_perf <- gp_test_perf %>% mutate(Model = "ガウス過程")

    plot_data <- bind_rows(lm_perf, rf_perf, gp_perf)

    ggplot(plot_data, aes(x = PC, y = R_squared, fill = Model)) +
      geom_col(position = "dodge", alpha = 0.8) +
      geom_text(aes(label = sprintf("%.3f", R_squared)),
                position = position_dodge(width = 0.9),
                vjust = -0.5, size = 3.5) +
      scale_fill_manual(values = c("重回帰分析" = "steelblue",
                                     "ランダムフォレスト" = "darkgreen",
                                     "ガウス過程" = "purple")) +
      ylim(0, 1) +
      labs(
        title = "モデル予測精度比較（テストデータ）",
        subtitle = "R²スコア - 値が高いほど予測精度が高い",
        x = "主成分", y = "R² スコア"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom"
      )
  })

  output$variable_importance_plot <- renderPlot({
    importance_data <- data.frame()

    if (input$model_type == "lm") {
      # 重回帰分析：係数の絶対値を重要度とする
      for (pc in paste0("PC", 1:M)) {
        model_coefs <- coef(lm_models[[pc]])
        # 切片を除く
        factor_coefs <- model_coefs[c("F1", "F2", "F3")]
        imp_df <- data.frame(
          PC = pc,
          Variable = c("F1", "F2", "F3"),
          Importance = abs(factor_coefs)
        )
        importance_data <- bind_rows(importance_data, imp_df)
      }

      ggplot(importance_data, aes(x = Variable, y = Importance, fill = PC)) +
        geom_col(position = "dodge", alpha = 0.8) +
        labs(
          title = "重回帰分析：変数重要度",
          subtitle = "係数の絶対値（値が高いほど予測に重要）",
          x = "因子", y = "重要度（|係数|）"
        ) +
        theme_minimal(base_size = 14) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "bottom"
        )

    } else if (input$model_type == "rf") {
      # ランダムフォレスト：%IncMSEを使用
      for (pc in paste0("PC", 1:M)) {
        imp <- importance(rf_models[[pc]])
        imp_df <- data.frame(
          PC = pc,
          Variable = rownames(imp),
          Importance = imp[, "%IncMSE"]
        )
        importance_data <- bind_rows(importance_data, imp_df)
      }

      ggplot(importance_data, aes(x = Variable, y = Importance, fill = PC)) +
        geom_col(position = "dodge", alpha = 0.8) +
        labs(
          title = "ランダムフォレスト：変数重要度",
          subtitle = "値が高いほど予測に重要",
          x = "因子", y = "重要度（%IncMSE）"
        ) +
        theme_minimal(base_size = 14) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "bottom"
        )

    } else if (input$model_type == "gp") {
      # ガウス過程：SHAP値を使用
      for (pc in paste0("PC", 1:M)) {
        # GP予測関数のラッパー
        pfun <- function(object, newdata) {
          newdata_matrix <- as.matrix(newdata)
          pred_result <- predict(object, newdata = newdata_matrix, type = "UK")
          return(pred_result$mean)
        }

        # SHAP値を計算（訓練データのサブセットを使用）
        shap_values <- explain(
          gp_models[[pc]],
          X = train_data[, c("F1", "F2", "F3")],
          pred_wrapper = pfun,
          nsim = 10,  # シミュレーション回数
          adjust = TRUE
        )

        # 各変数のSHAP値の絶対値平均を重要度とする
        shap_importance <- colMeans(abs(shap_values))
        imp_df <- data.frame(
          PC = pc,
          Variable = names(shap_importance),
          Importance = as.numeric(shap_importance)
        )
        importance_data <- bind_rows(importance_data, imp_df)
      }

      ggplot(importance_data, aes(x = Variable, y = Importance, fill = PC)) +
        geom_col(position = "dodge", alpha = 0.8) +
        labs(
          title = "ガウス過程：変数重要度（SHAP値）",
          subtitle = "値が高いほど予測に重要",
          x = "因子", y = "重要度（|SHAP|平均）"
        ) +
        theme_minimal(base_size = 14) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "bottom"
        )
    }
  })

  output$selected_model_text <- renderText({
    if (input$model_type == "lm") {
      "重回帰分析（線形モデル）"
    } else if (input$model_type == "rf") {
      "ランダムフォレスト（非線形モデル）"
    } else {
      "ガウス過程（確率的非線形モデル）"
    }
  })

  # タブ3: 因子調整と最適化 ---------------------------------

  output$mse_box <- renderValueBox({
    valueBox(
      sprintf("%.6f", metrics()$mse),
      "MSE (平均二乗誤差)",
      icon = icon("calculator"),
      color = "red"
    )
  })

  output$rmse_box <- renderValueBox({
    valueBox(
      sprintf("%.6f", metrics()$rmse),
      "RMSE (二乗平均平方根誤差)",
      icon = icon("ruler"),
      color = "orange"
    )
  })

  output$cor_box <- renderValueBox({
    valueBox(
      sprintf("%.4f", metrics()$correlation),
      "相関係数",
      icon = icon("heart"),
      color = if(metrics()$correlation > 0.9) "green" else if(metrics()$correlation > 0.7) "yellow" else "red"
    )
  })

  # rmse_per_mic_plotは削除（一次元時系列のため不要）

  output$waveform_comparison <- renderPlot({
    # 現在の波形と目標波形、最適化結果がある場合は最適化波形も表示
    plot_data <- data.frame()

    # 目標波形
    df_target <- data.frame(
      time = time_grid,
      amplitude = target_wave(),
      type = "目標波形"
    )

    # 現在の波形
    df_current <- data.frame(
      time = time_grid,
      amplitude = current_wave(),
      type = "現在の波形"
    )

    plot_data <- rbind(plot_data, df_target, df_current)

    # 最適化結果がある場合
    if (!is.null(rv$opt_factors)) {
      opt_wave <- generate_waveform(rv$opt_factors[1], rv$opt_factors[2], rv$opt_factors[3])
      df_opt <- data.frame(
        time = time_grid,
        amplitude = opt_wave,
        type = "最適化波形"
      )
      plot_data <- rbind(plot_data, df_opt)
    }

    # プロット
    colors <- c("目標波形" = "black", "現在の波形" = "red", "最適化波形" = "green")
    linetypes <- c("目標波形" = "solid", "現在の波形" = "dashed", "最適化波形" = "dotdash")

    ggplot(plot_data, aes(x = time, y = amplitude, color = type, linetype = type)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = colors) +
      scale_linetype_manual(values = linetypes) +
      labs(
        title = "波形比較",
        x = "時間", y = "振幅",
        color = "波形タイプ", linetype = "波形タイプ"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom"
      )
  })

  observeEvent(input$reset_factors, {
    updateSliderInput(session, "f1", value = 0)
    updateSliderInput(session, "f2", value = 0)
    updateSliderInput(session, "f3", value = 0)
  })

  # 最適値探索
  observeEvent(input$run_optimization, {
    withProgress(message = '最適値探索実行中...', value = 0, {
      incProgress(0.3, detail = "目的関数を評価中")

      # アルゴリズムと目的関数の組み合わせで処理を分岐
      if (input$opt_algorithm == "lbfgs") {
        # L-BFGS-Bアルゴリズム
        if (input$opt_method == "direct") {
          # 直接法（波形ベース）
          opt_result <- optimize_factors_direct(target_wave())
          rv$opt_method_used <- "direct_lbfgs"
        } else {
          # モデルベース（PCスコア）
          if (input$model_type == "lm") {
            opt_result <- optimize_factors_model(target_pc(), "lm", lm_models)
          } else if (input$model_type == "rf") {
            opt_result <- optimize_factors_model(target_pc(), "rf", rf_models)
          } else {
            opt_result <- optimize_factors_model(target_pc(), "gp", gp_models)
          }
          rv$opt_method_used <- paste0("model_lbfgs_", input$model_type)
        }
      } else {
        # ベイズ最適化
        if (input$opt_method == "direct") {
          # 直接法（波形ベース）
          opt_result <- optimize_factors_bayesian_direct(target_wave())
          rv$opt_method_used <- "direct_bayes"
        } else {
          # モデルベース（PCスコア）
          if (input$model_type == "lm") {
            opt_result <- optimize_factors_bayesian_model(target_pc(), "lm", lm_models)
          } else if (input$model_type == "rf") {
            opt_result <- optimize_factors_bayesian_model(target_pc(), "rf", rf_models)
          } else {
            opt_result <- optimize_factors_bayesian_model(target_pc(), "gp", gp_models)
          }
          rv$opt_method_used <- paste0("model_bayes_", input$model_type)
        }
      }

      incProgress(0.6, detail = "最適解を計算中")

      rv$opt_factors <- opt_result$factors
      rv$opt_mse <- opt_result$mse

      incProgress(1, detail = "完了！")
    })

    showNotification("最適値探索が完了しました！", type = "message")
  })

  output$opt_result <- renderPrint({
    if (is.null(rv$opt_factors)) {
      cat("最適値探索を実行してください。\n")
    } else {
      cat("=== 最適化結果 ===\n\n")

      # 使用した手法を表示
      if (grepl("bayes", rv$opt_method_used)) {
        cat("アルゴリズム: ベイズ最適化\n")
      } else {
        cat("アルゴリズム: L-BFGS-B\n")
      }

      if (grepl("direct", rv$opt_method_used)) {
        cat("目的関数: 直接法（波形ベース）\n\n")
      } else {
        model_name <- sub(".*_(lm|rf|gp)", "\\1", rv$opt_method_used)
        model_full <- switch(model_name,
                            "lm" = "重回帰分析",
                            "rf" = "ランダムフォレスト",
                            "gp" = "ガウス過程")
        cat(sprintf("目的関数: モデルベース（%s）\n\n", model_full))
      }

      cat(sprintf("F1: %.4f\n", rv$opt_factors[1]))
      cat(sprintf("F2: %.4f\n", rv$opt_factors[2]))
      cat(sprintf("F3: %.4f\n", rv$opt_factors[3]))

      if (grepl("direct", rv$opt_method_used)) {
        cat(sprintf("\n波形MSE: %.6f\n", rv$opt_mse))
        cat(sprintf("波形RMSE: %.6f\n", sqrt(rv$opt_mse)))
      } else {
        cat(sprintf("\nPCスコアMSE: %.6f\n", rv$opt_mse))

        # 実際の波形誤差も計算
        opt_wave <- generate_waveform(rv$opt_factors[1], rv$opt_factors[2], rv$opt_factors[3])
        wave_mse <- mean((opt_wave - target_wave())^2)
        cat(sprintf("波形MSE: %.6f\n", wave_mse))
        cat(sprintf("波形RMSE: %.6f\n", sqrt(wave_mse)))
      }
    }
  })

  observeEvent(input$apply_opt_factors, {
    if (!is.null(rv$opt_factors)) {
      updateSliderInput(session, "f1", value = rv$opt_factors[1])
      updateSliderInput(session, "f2", value = rv$opt_factors[2])
      updateSliderInput(session, "f3", value = rv$opt_factors[3])
      showNotification("最適化結果を適用しました！", type = "message")
    } else {
      showNotification("まず最適値探索を実行してください。", type = "warning")
    }
  })

  # タブ4: 感度分析 ----------------------------------------

  output$sensitivity_plot <- renderPlot({
    # 基準因子（現在の設定）
    baseline_factors <- data.frame(F1 = input$f1, F2 = input$f2, F3 = input$f3)

    sensitivity_range <- seq(-0.3, 0.3, length.out = 21)
    sensitivity_results <- list()

    for (factor_name in c("F1", "F2", "F3")) {
      results_df <- data.frame()

      for (delta in sensitivity_range) {
        test_factors <- baseline_factors
        test_factors[[factor_name]] <- test_factors[[factor_name]] * (1 + delta)

        test_wave <- generate_waveform(test_factors$F1, test_factors$F2, test_factors$F3)
        rmse_total <- sqrt(mean((test_wave - target_wave())^2))

        results_df <- rbind(results_df, data.frame(
          Factor = factor_name,
          Delta_Percent = delta * 100,
          RMSE = rmse_total
        ))
      }

      sensitivity_results[[factor_name]] <- results_df
    }

    sensitivity_df <- bind_rows(sensitivity_results)

    ggplot(sensitivity_df, aes(x = Delta_Percent, y = RMSE, color = Factor)) +
      geom_line(linewidth = 1.5) +
      geom_point(size = 2) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 1) +
      labs(
        title = "因子変動に対する波形誤差の感度分析",
        subtitle = "現在の因子設定からの変化率 vs RMSE",
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
  })

  output$tolerance_table <- renderTable({
    # 簡易的な許容範囲テーブル
    data.frame(
      Factor = c("F1", "F2", "F3"),
      Current_Value = c(input$f1, input$f2, input$f3),
      Recommended_Range = c("±10%", "±10%", "±10%"),
      Sensitivity = c("中", "高", "低")
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
}

############################################################
# アプリ起動
############################################################

shinyApp(ui = ui, server = server)