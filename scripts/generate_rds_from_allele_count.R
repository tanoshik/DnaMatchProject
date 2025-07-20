# scripts/generate_rds_from_allele_count.R
# 使用目的: Allele-count_GF_Japanese.csv から locus_order.rds, freq_table.rds を生成
# ※ このファイルは source() するだけでは何も起こらない（関数定義のみ）

generate_rds_interactive <- function() {
  suppressPackageStartupMessages({
    library(tidyverse)
    library(tcltk)
  })
  
  # ファイル選択（ユーザーへの説明あり）
  file_path <- tcltk::tk_choose.files(
    caption = "Select allele count CSV file (e.g., Allele-count_GF_Japanese.csv)",
    multi = FALSE
  )
  
  if (length(file_path) == 0 || file_path == "") {
    stop("No file selected. Please select a valid allele count CSV file.")
  }
  
  # CSV読み込み
  allele_raw <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  allele_raw[is.na(allele_raw)] <- 0
  allele_raw$Allele <- as.character(allele_raw$Allele)
  
  # Locus順序
  locus_order <- colnames(allele_raw)[-1]
  
  # ワイド形式 → ロング形式変換
  allele_long <- allele_raw %>%
    pivot_longer(cols = -Allele, names_to = "Locus", values_to = "Count") %>%
    mutate(Count = as.numeric(Count)) %>%
    filter(Count > 0)
  
  # 頻度計算
  freq_table <- allele_long %>%
    group_by(Locus) %>%
    mutate(Freq = Count / sum(Count)) %>%
    ungroup()
  
  # 保存パス
  save_dir <- "data"
  dir.create(save_dir, showWarnings = FALSE)
  
  # RDS保存
  saveRDS(locus_order, file = file.path(save_dir, "locus_order.rds"))
  saveRDS(freq_table, file = file.path(save_dir, "freq_table.rds"))
  
  # 確認出力（ASCIIのみ）
  cat("RDS files successfully saved:\n")
  cat(" -", file.path(save_dir, "locus_order.rds"), "\n")
  cat(" -", file.path(save_dir, "freq_table.rds"), "\n")
}

# 非対話版 generate_rds（main.R から呼び出す用）
generate_rds <- function(file_path = "data/Allele-count_GF_Japanese.csv") {
  suppressPackageStartupMessages(library(tidyverse))
  allele_raw <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  allele_raw[is.na(allele_raw)] <- 0
  allele_raw$Allele <- as.character(allele_raw$Allele)
  
  locus_order <- colnames(allele_raw)[-1]
  
  allele_long <- allele_raw %>%
    pivot_longer(cols = -Allele, names_to = "Locus", values_to = "Count") %>%
    mutate(Count = as.numeric(Count)) %>%
    filter(Count > 0)
  
  freq_table <- allele_long %>%
    group_by(Locus) %>%
    mutate(Freq = Count / sum(Count)) %>%
    ungroup()
  
  dir.create("data", showWarnings = FALSE)
  saveRDS(locus_order, "data/locus_order.rds")
  saveRDS(freq_table, "data/freq_table.rds")
  cat("[OK] freq_table.rds and locus_order.rds have been generated.\n")
}
