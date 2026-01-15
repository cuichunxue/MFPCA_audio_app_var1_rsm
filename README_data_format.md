# データフォーマットガイド

実データでFPCA解析を実行するためのデータ準備ガイドです。

## サポートされるデータ形式

3つのデータ形式に対応しています：

### 形式1: ワイド形式（推奨）

**1つのCSVファイルに因子と時系列を結合**

```csv
cond,F1,F2,F3,t_0.00,t_0.01,t_0.02,...,t_1.00
1,0.5,-0.3,0.2,0.123,0.145,0.167,...,0.234
2,-0.2,0.6,-0.1,0.234,0.256,0.278,...,0.345
3,0.8,-0.5,0.4,0.156,0.178,0.189,...,0.267
...
```

**列の説明:**
- `cond`: 条件ID（1, 2, 3, ...）
- `F1, F2, F3`: 因子の値（温度、圧力、流量など）
- `t_0.00, t_0.01, ...`: 時系列データ（列名が時刻を表す）

**読み込み方法:**
```r
source("data_loader.R")
data <- load_data_wide(
  file_path = "your_data.csv",
  factor_cols = c("F1", "F2", "F3"),
  time_prefix = "t_"
)
factors <- data$factors
Y_matrix <- data$Y_matrix
time_grid <- data$time_grid
```

---

### 形式2: 分離形式

**因子と時系列を別々のファイルで管理**

**factors.csv:**
```csv
cond,F1,F2,F3
1,0.5,-0.3,0.2
2,-0.2,0.6,-0.1
3,0.8,-0.5,0.4
...
```

**timeseries.csv:**
```csv
cond,t_0.00,t_0.01,t_0.02,...,t_1.00
1,0.123,0.145,0.167,...,0.234
2,0.234,0.256,0.278,...,0.345
3,0.156,0.178,0.189,...,0.267
...
```

**読み込み方法:**
```r
source("data_loader.R")
data <- load_data_separate(
  factor_file = "factors.csv",
  timeseries_file = "timeseries.csv",
  time_prefix = "t_"
)
factors <- data$factors
Y_matrix <- data$Y_matrix
time_grid <- data$time_grid
```

---

### 形式3: 時系列行列のみ

**時系列データのみの行列（因子は後で追加）**

**timeseries_matrix.csv:**
```csv
0.123,0.145,0.167,...,0.234
0.234,0.256,0.278,...,0.345
0.156,0.178,0.189,...,0.267
...
```

**読み込み方法:**
```r
source("data_loader.R")
ts_data <- load_timeseries_only(
  file_path = "timeseries_matrix.csv",
  time_range = c(0, 1)
)
Y_matrix <- ts_data$Y_matrix
time_grid <- ts_data$time_grid

# 因子は別途作成
factors <- data.frame(
  cond = 1:nrow(Y_matrix),
  Temperature = c(300, 310, 320, ...),
  Pressure = c(1.0, 1.2, 1.5, ...)
)
```

---

## データ準備のチェックリスト

### 1. 時系列データ
- [ ] 各行が1つの条件（実験条件、測定条件など）に対応
- [ ] 各列が1つの時間点に対応
- [ ] 欠損値がない（あればNAで埋める）
- [ ] 時間間隔は等間隔が望ましい

### 2. 因子データ
- [ ] 条件数が時系列データの行数と一致
- [ ] 因子は連続値（温度、圧力など）またはダミー変数
- [ ] 因子名は英数字（日本語も可だが英語推奨）

### 3. 列名
- [ ] 時系列列は共通のプレフィックスを使用（例: `t_`, `time_`）
- [ ] 列名に空白やカンマを含まない

---

## サンプルデータの生成

テスト用にサンプルデータを生成できます：

```r
source("data_loader.R")

# ワイド形式
create_sample_data(output_dir = "sample_data", format = "wide")

# 分離形式
create_sample_data(output_dir = "sample_data", format = "separate")

# 行列形式
create_sample_data(output_dir = "sample_data", format = "matrix")
```

生成されたサンプルデータで動作確認してから、実データに置き換えることを推奨します。

---

## 実データでの解析手順

### ステップ1: データファイルを準備

いずれかの形式でCSVファイルを準備します。Excelで作成する場合は「CSV UTF-8（コンマ区切り）」で保存してください。

### ステップ2: データローダーでデータを読み込み

```r
source("data_loader.R")

# 例: ワイド形式の場合
data <- load_data_wide(
  file_path = "your_data.csv",
  factor_cols = c("Temperature", "Pressure", "FlowRate"),
  time_prefix = "t_"
)

factors <- data$factors
Y_matrix <- data$Y_matrix
time_grid <- data$time_grid
```

### ステップ3: MFPCA.Rで解析実行

```r
# MFPCA.Rの冒頭のダミーデータ生成部分をコメントアウト
# （31行目〜75行目をスキップ）

# データローダーから読み込んだデータを使用
source("MFPCA.R")
```

または、`example_use_real_data.R`を参考にして独自のスクリプトを作成することもできます。

---

## データの例

### 音響データの例

**実験設定:**
- 6つの温度条件（290K, 300K, 310K, 320K, 330K, 340K）
- 各条件で音響信号を100時間点で測定

**データ形式（ワイド形式）:**
```csv
cond,Temperature,t_0.00,t_0.01,t_0.02,...,t_0.99
1,290,0.15,0.18,0.22,...,0.31
2,300,0.18,0.21,0.25,...,0.34
3,310,0.21,0.24,0.28,...,0.37
...
```

### センサーデータの例

**実験設定:**
- 温度、圧力、流量の3因子を変化
- 各条件でセンサー値を200時間点で記録

**データ形式（分離形式）:**

**factors.csv:**
```csv
cond,Temperature,Pressure,FlowRate
1,300,1.0,10.0
2,310,1.0,10.0
3,300,1.2,10.0
...
```

**timeseries.csv:**
```csv
cond,t_0.000,t_0.005,t_0.010,...,t_0.995
1,25.3,25.5,25.7,...,26.1
2,26.1,26.3,26.5,...,26.9
3,25.8,26.0,26.2,...,26.6
...
```

---

## トラブルシューティング

### エラー: "時系列列が見つかりません"

**原因:** `time_prefix` の設定が列名と一致していない

**解決策:**
```r
# 列名が "time_0.00", "time_0.01" の場合
data <- load_data_wide("your_data.csv", time_prefix = "time_")

# 列名が "T0", "T1", "T2" の場合
data <- load_data_wide("your_data.csv", time_prefix = "T")
```

### エラー: "因子データと時系列データの条件IDが一致しません"

**原因:** 因子ファイルと時系列ファイルの行の順序が異なる

**解決策:**
両方のファイルで `cond` 列を基準にソートしてから保存してください。

### エラー: 文字化け

**原因:** CSVファイルのエンコーディングが正しくない

**解決策:**
```r
# UTF-8で明示的に読み込み
data <- read.csv("your_data.csv", fileEncoding = "UTF-8")
```

---

## まとめ

1. **データをCSV形式で準備**（3つの形式から選択）
2. **`data_loader.R` でデータを読み込み**
3. **`MFPCA.R` または `example_use_real_data.R` で解析実行**

不明点があれば、`example_use_real_data.R` のコードを参考にしてください。
