## simProfiles_simDNAmixtures.R
## Description: Simulate single-source DNA profiles using simDNAmixtures

# 必要なパッケージを読み込み
library(simDNAmixtures)  # プロファイルシミュレーション本体
library(dplyr)           # データ整形
library(readr)           # CSV読み書き
library(tidyr)           # wide ↔ long 変換

# 1. 日本人GlobalFilerアレル頻度表（wide形式）を読み込む
allele_count_raw <- read_csv("data/Allele-count_GF_Japanese.csv", show_col_types = FALSE)

# 2. long形式に変換 → 各ローカス・アレルの相対頻度を計算
allele_freq <- allele_count_raw %>%
  pivot_longer(-Allele, names_to = "Marker", values_to = "Count") %>%
  group_by(Marker) %>%
  mutate(Frequency = Count / sum(Count, na.rm = TRUE)) %>%
  filter(!is.na(Frequency)) %>%
  select(Marker, Allele, Frequency) %>%
  arrange(Marker, Allele)

# 3. 中間ファイルとしてCSVに保存（他でも使えるように）
output_freq_path <- "data/kit/Allele-count_GF_Japanese_simDNA.csv"
write_csv(allele_freq, output_freq_path)

# 4. 再読み込みし、simDNAmixtures形式（Marker×Allele）に変換
freqs_raw <- read_csv(output_freq_path, show_col_types = FALSE)

# 5. AMEL（性染色体マーカー）を除外
freqs_filtered <- freqs_raw %>% filter(Marker != "AMEL")

# 6. list形式（Marker = named vector of frequencies）に変換
freqs <- freqs_filtered %>%
  group_by(Marker) %>%
  summarise(freq_vector = list(setNames(Frequency, Allele)), .groups = "drop") %>%
  deframe()

# 7. GlobalFiler設定（stutterなし）を読み込み、AMELを除外
kit_conf <- gf_configuration()
model_settings <- kit_conf$gamma_settings_no_stutter
model_settings$locus_names <- setdiff(model_settings$locus_names, "AMEL")  # 念のためAMELを除外

# 8. サンプリングパラメータ設定
sampling_parameters <- list(
  min_mu = 100, max_mu = 2000,              # DNA量（テンプレート）範囲（pg）
  min_cv = 0.1, max_cv = 0.3,               # CV（高さばらつき）
  degradation_shape1 = 10, degradation_shape2 = 1  # 分解度（ガンマ分布）
)

# 9. 単一寄与者のプロファイルを50件シミュレーション
sim_result <- sample_mixtures(
  n = 50,
  contributors = c("U1"),  # 寄与者ラベル
  freqs = freqs,           # アレル頻度リスト
  sampling_parameters = sampling_parameters,
  model_settings = model_settings,
  sample_model = sample_gamma_model         # ピーク高さにガンマ分布を使用
)

# 10. 出力ファイル名にタイムスタンプ（例：simulated_epg_202507251210.csv）
timestamp <- format(Sys.time(), "%Y%m%d%H%M")
output_csv_path <- sprintf("data/simulated_epg_%s.csv", timestamp)

# 11. SMASH形式（long形式）のプロファイルを保存
write_csv(sim_result$smash, output_csv_path)

# 12. 完了メッセージ表示
cat(sprintf("Simulated profiles saved to %s\n", output_csv_path))
