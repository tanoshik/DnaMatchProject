# 欠損補完とホモ型変換処理
prepare_profile <- function(profile, homo_to_any = FALSE) {
  lapply(profile, function(alleles) {
    alleles <- as.character(alleles)
    if (length(alleles) != 2) alleles <- c("any", "any")
    alleles[is.na(alleles) | alleles == ""] <- "any"
    if (homo_to_any && alleles[1] == alleles[2] && alleles[1] != "any") {
      alleles[2] <- "any"
    }
    alleles
  })
}
# 検索クエリの準備
prepare_query <- function(query_file = "data/query_profile.csv",
                          locus_file = "data/locus_order.rds",
                          homo_to_any = TRUE) {
  # locus_order を読み込み
  if (!file.exists(locus_file)) {
    stop(paste("ローカス順ファイルが存在しません:", locus_file))
  }
  locus_order <- readRDS(locus_file)

  # query_profile 読み込み（補完付き）
  query_profile <- read_query_profile(query_file, locus_order, homo_to_any) # nolint: object_usage_linter, line_length_linter.

  # ログメッセージ（除外されたローカス）
  df_raw <- read.csv(query_file, stringsAsFactors = FALSE)
  used_loci <- intersect(df_raw$Locus, locus_order)
  unused_loci <- setdiff(df_raw$Locus, used_loci)
  missing_loci <- setdiff(locus_order, used_loci)

  if (length(unused_loci) > 0) {
    cat("[Warning] Excluded loci:", paste(unused_loci, collapse = ", "), "\n")
  }
  if (length(missing_loci) > 0) {
    cat("[Info] Imputed loci:", paste(missing_loci, collapse = ", "), "\n")
  }

  return(query_profile)
}
# データベースの準備（補完・整形・ローカス順に揃える）
prepare_database <- function(db_file = "data/database_profile.csv",
                             locus_file = "data/locus_order.rds",
                             homo_to_any = FALSE) {
  # ローカス順を読み込み
  if (!file.exists(locus_file)) {
    stop(paste("ローカス順ファイルが存在しません:", locus_file))
  }
  locus_order <- readRDS(locus_file)

  # データベース読み込み
  db_profiles <- read_db_profiles(db_file, locus_order, homo_to_any) # nolint: object_usage_linter, line_length_linter.

  # ログ：ローカスの過不足チェック（SampleID 1件目を対象に確認）
  df_raw <- read.csv(db_file, stringsAsFactors = FALSE)
  sample1 <- unique(df_raw$SampleID)[1]
  df1 <- df_raw[df_raw$SampleID == sample1, ]
  used_loci <- intersect(df1$Locus, locus_order)
  unused_loci <- setdiff(df1$Locus, used_loci)
  missing_loci <- setdiff(locus_order, used_loci)

  if (length(unused_loci) > 0) {
    cat("⚠️ 除外されたローカス（DB）:", paste(unused_loci, collapse = ", "), "\n")
  }
  if (length(missing_loci) > 0) {
    cat("🔧 補完されたローカス（DB）:", paste(missing_loci, collapse = ", "), "\n")
  }

  return(db_profiles)
}
