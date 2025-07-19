# ファイルのパス
file_path <- "data/database_profile_test.csv"
locus_order_path <- "data/locus_order.rds"

# ローカス順を読み込み
locus_order <- readRDS(locus_order_path)

# 修正版 read_db_profiles を再定義（必要ならここに貼る）
read_db_profiles <- function(file_path, locus_order, homo_to_any = FALSE) {
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  sample_ids <- unique(df$SampleID)
  
  profiles <- lapply(sample_ids, function(sid) {
    sub_df <- df[df$SampleID == sid, ]
    
    raw_profile <- setNames(
      lapply(locus_order, function(locus) {
        row <- sub_df[sub_df$Locus == locus, c("allele1", "allele2")]
        if (nrow(row) == 0) return(c("any", "any"))
        alleles <- unlist(row[1, ])
        alleles[is.na(alleles) | alleles == ""] <- "any"
        return(alleles)
      }),
      locus_order
    )
    
    prepare_profile(raw_profile, homo_to_any)
  })
  
  names(profiles) <- sample_ids
  profiles
}

# テスト実行
db_result <- read_db_profiles(file_path, locus_order, homo_to_any = FALSE)

# 内容を確認
str(db_result$Sample1)     # Sample1 のプロファイル構造
print(db_result$Sample1)   # 中身を表示（補完・整列の確認用）

