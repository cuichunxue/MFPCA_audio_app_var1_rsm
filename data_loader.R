############################################################
# データ読み込みユーティリティ
# 実データを解析するためのデータローダー
############################################################

library(dplyr)

#' データ形式1: ワイド形式（因子と時系列が1行に結合）
#'
#' CSV形式例:
#' cond, F1, F2, F3, t_0.00, t_0.01, t_0.02, ..., t_1.00
#' 1, 0.5, -0.3, 0.2, 0.123, 0.145, 0.167, ..., 0.234
#' 2, -0.2, 0.6, -0.1, 0.234, 0.256, 0.278, ..., 0.345
#'
#' @param file_path CSVファイルのパス
#' @param factor_cols 因子列の名前（例: c("F1", "F2", "F3")）
#' @param time_prefix 時系列列のプレフィックス（例: "t_"）
#' @return list(factors = データフレーム, Y_matrix = 行列, time_grid = ベクトル)
load_data_wide <- function(file_path, factor_cols = c("F1", "F2", "F3"),
                           time_prefix = "t_") {
  cat("ワイド形式データを読み込み中...\n")
  cat(sprintf("  ファイル: %s\n", file_path))

  # データ読み込み
  data <- read.csv(file_path, stringsAsFactors = FALSE)

  # 因子データを抽出
  if (!"cond" %in% colnames(data)) {
    data$cond <- 1:nrow(data)
  }
  factors <- data[, c("cond", factor_cols), drop = FALSE]

  # 時系列列を特定
  time_cols <- grep(paste0("^", time_prefix), colnames(data), value = TRUE)

  if (length(time_cols) == 0) {
    stop("時系列列が見つかりません。time_prefix を確認してください。")
  }

  # 時系列データを行列に変換
  Y_matrix <- as.matrix(data[, time_cols])

  # 時間グリッドを抽出（列名から）
  time_grid <- as.numeric(sub(paste0("^", time_prefix), "", time_cols))

  cat(sprintf("  条件数: %d\n", nrow(factors)))
  cat(sprintf("  因子数: %d\n", length(factor_cols)))
  cat(sprintf("  時間点数: %d\n", length(time_grid)))
  cat(sprintf("  時間範囲: [%.3f, %.3f]\n", min(time_grid), max(time_grid)))

  return(list(
    factors = factors,
    Y_matrix = Y_matrix,
    time_grid = time_grid
  ))
}


#' データ形式2: 因子と時系列を別々のファイルで管理
#'
#' factors.csv:
#' cond, F1, F2, F3
#' 1, 0.5, -0.3, 0.2
#' 2, -0.2, 0.6, -0.1
#'
#' timeseries.csv:
#' cond, t_0.00, t_0.01, t_0.02, ..., t_1.00
#' 1, 0.123, 0.145, 0.167, ..., 0.234
#' 2, 0.234, 0.256, 0.278, ..., 0.345
#'
#' @param factor_file 因子データのCSVファイルパス
#' @param timeseries_file 時系列データのCSVファイルパス
#' @param time_prefix 時系列列のプレフィックス（例: "t_"）
#' @return list(factors = データフレーム, Y_matrix = 行列, time_grid = ベクトル)
load_data_separate <- function(factor_file, timeseries_file, time_prefix = "t_") {
  cat("別々のファイルからデータを読み込み中...\n")
  cat(sprintf("  因子ファイル: %s\n", factor_file))
  cat(sprintf("  時系列ファイル: %s\n", timeseries_file))

  # 因子データ読み込み
  factors <- read.csv(factor_file, stringsAsFactors = FALSE)

  # 時系列データ読み込み
  timeseries_data <- read.csv(timeseries_file, stringsAsFactors = FALSE)

  # 時系列列を特定
  time_cols <- grep(paste0("^", time_prefix), colnames(timeseries_data), value = TRUE)

  if (length(time_cols) == 0) {
    stop("時系列列が見つかりません。time_prefix を確認してください。")
  }

  # 時系列データを行列に変換
  Y_matrix <- as.matrix(timeseries_data[, time_cols])

  # 時間グリッドを抽出
  time_grid <- as.numeric(sub(paste0("^", time_prefix), "", time_cols))

  # 条件IDでマージ確認
  if ("cond" %in% colnames(factors) && "cond" %in% colnames(timeseries_data)) {
    if (!all(factors$cond == timeseries_data$cond)) {
      warning("因子データと時系列データの条件IDが一致しません。順序を確認してください。")
    }
  }

  cat(sprintf("  条件数: %d\n", nrow(factors)))
  cat(sprintf("  因子数: %d\n", ncol(factors) - 1))  # condを除く
  cat(sprintf("  時間点数: %d\n", length(time_grid)))
  cat(sprintf("  時間範囲: [%.3f, %.3f]\n", min(time_grid), max(time_grid)))

  return(list(
    factors = factors,
    Y_matrix = Y_matrix,
    time_grid = time_grid
  ))
}


#' データ形式3: 時系列のみの行列ファイル（因子なし、後で追加）
#'
#' timeseries.csv:
#' 0.123, 0.145, 0.167, ..., 0.234
#' 0.234, 0.256, 0.278, ..., 0.345
#'
#' @param file_path 時系列行列のCSVファイルパス
#' @param time_range 時間範囲（例: c(0, 1)）
#' @return list(Y_matrix = 行列, time_grid = ベクトル)
load_timeseries_only <- function(file_path, time_range = c(0, 1)) {
  cat("時系列行列データを読み込み中...\n")
  cat(sprintf("  ファイル: %s\n", file_path))

  # データ読み込み
  Y_matrix <- as.matrix(read.csv(file_path, header = FALSE))

  # 時間グリッドを生成
  n_time <- ncol(Y_matrix)
  time_grid <- seq(time_range[1], time_range[2], length.out = n_time)

  cat(sprintf("  条件数: %d\n", nrow(Y_matrix)))
  cat(sprintf("  時間点数: %d\n", n_time))
  cat(sprintf("  時間範囲: [%.3f, %.3f]\n", time_range[1], time_range[2]))

  return(list(
    Y_matrix = Y_matrix,
    time_grid = time_grid
  ))
}


#' サンプルデータを生成してCSV形式で保存
#'
#' @param output_dir 出力ディレクトリ
#' @param format データ形式（"wide", "separate", "matrix"）
create_sample_data <- function(output_dir = "sample_data", format = "wide") {
  cat("サンプルデータを生成中...\n")

  # ディレクトリ作成
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # サンプルデータ生成
  set.seed(123)
  N_cond <- 20
  n_time <- 100
  time_grid <- seq(0, 1, length.out = n_time)

  factors <- data.frame(
    cond = 1:N_cond,
    F1 = runif(N_cond, -1, 1),
    F2 = runif(N_cond, -1, 1),
    F3 = runif(N_cond, -1, 1)
  )

  # 波形生成
  base_wave <- sin(2 * pi * time_grid) + 0.3 * sin(6 * pi * time_grid)
  shape_F1 <- sin(4 * pi * time_grid)
  shape_F2 <- (time_grid - 0.5)
  shape_F3 <- (1 - (time_grid - 0.5)^2)

  Y_matrix <- matrix(NA_real_, nrow = N_cond, ncol = n_time)
  for (i in 1:N_cond) {
    signal_i <- base_wave +
      factors$F1[i] * 0.5 * shape_F1 +
      factors$F2[i] * 0.4 * shape_F2 +
      factors$F3[i] * 0.3 * shape_F3
    Y_matrix[i, ] <- signal_i + rnorm(n_time, 0, 0.1)
  }

  # 形式に応じて保存
  if (format == "wide") {
    # ワイド形式
    time_cols <- paste0("t_", sprintf("%.2f", time_grid))
    colnames(Y_matrix) <- time_cols
    data_wide <- cbind(factors, Y_matrix)

    output_file <- file.path(output_dir, "data_wide.csv")
    write.csv(data_wide, output_file, row.names = FALSE)
    cat(sprintf("  保存: %s\n", output_file))

  } else if (format == "separate") {
    # 分離形式
    factor_file <- file.path(output_dir, "factors.csv")
    write.csv(factors, factor_file, row.names = FALSE)
    cat(sprintf("  保存: %s\n", factor_file))

    time_cols <- paste0("t_", sprintf("%.2f", time_grid))
    colnames(Y_matrix) <- time_cols
    timeseries_data <- cbind(cond = factors$cond, Y_matrix)

    timeseries_file <- file.path(output_dir, "timeseries.csv")
    write.csv(timeseries_data, timeseries_file, row.names = FALSE)
    cat(sprintf("  保存: %s\n", timeseries_file))

  } else if (format == "matrix") {
    # 行列のみ
    output_file <- file.path(output_dir, "timeseries_matrix.csv")
    write.csv(Y_matrix, output_file, row.names = FALSE)
    cat(sprintf("  保存: %s\n", output_file))

    # 因子は別途保存
    factor_file <- file.path(output_dir, "factors.csv")
    write.csv(factors, factor_file, row.names = FALSE)
    cat(sprintf("  保存（因子）: %s\n", factor_file))
  }

  cat("\nサンプルデータの生成が完了しました！\n")
  cat(sprintf("出力ディレクトリ: %s\n", output_dir))

  # 使用例を表示
  cat("\n【使用例】\n")
  if (format == "wide") {
    cat(sprintf('data <- load_data_wide("%s")\n', output_file))
  } else if (format == "separate") {
    cat(sprintf('data <- load_data_separate("%s", "%s")\n', factor_file, timeseries_file))
  } else if (format == "matrix") {
    cat(sprintf('ts_data <- load_timeseries_only("%s")\n', output_file))
    cat(sprintf('factors <- read.csv("%s")\n', factor_file))
  }
}


############################################################
# 使用例
############################################################

if (FALSE) {
  # サンプルデータを生成
  create_sample_data(output_dir = "sample_data", format = "wide")
  create_sample_data(output_dir = "sample_data", format = "separate")

  # データを読み込み
  data1 <- load_data_wide("sample_data/data_wide.csv")
  data2 <- load_data_separate("sample_data/factors.csv", "sample_data/timeseries.csv")

  # MFPCA.R で使用
  factors <- data1$factors
  Y_matrix <- data1$Y_matrix
  time_grid <- data1$time_grid
}
