# 必要パッケージの読み込み
library(tidyverse)

# ファイル選択
file_path <- file.choose()
allele_raw <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)

# 欠損を0に置き換え
allele_raw[is.na(allele_raw)] <- 0

# "Allele" 列を文字列として保持
allele_raw$Allele <- as.character(allele_raw$Allele)

# Locusの順序を保存（ファイルの左から右の順）
locus_order <- colnames(allele_raw)[-1]  # "Allele"以外

# ワイド形式 → ロング形式へ変換
allele_long <- allele_raw %>%
  pivot_longer(
    cols = -Allele,
    names_to = "Locus",
    values_to = "Count"
  ) %>%
  mutate(
    Count = as.numeric(Count)
  ) %>%
  filter(Count > 0)  # 空白だったもの（0件）は削除

# 出現頻度（割合）の計算
freq_table <- allele_long %>%
  group_by(Locus) %>%
  mutate(Freq = Count / sum(Count)) %>%
  ungroup()

# 保存先（同じフォルダに保存）
save_dir <- dirname(file_path)
saveRDS(locus_order, file = file.path(save_dir, "data/locus_order.rds"))
saveRDS(freq_table, file = file.path(save_dir, "data/freq_table.rds"))

cat("✅ RDSファイルを保存しました。\n")
cat(" -", file.path(save_dir, "data/locus_order.rds"), "\n")
cat(" -", file.path(save_dir, "data/freq_table.rds"), "\n")
