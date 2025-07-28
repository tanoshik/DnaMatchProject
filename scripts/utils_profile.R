# リスト形式への整形ラッパー関数（内部で prepare_profile_df() を呼び出し）
# 欠損補完（any）とローカス順の補完を行い、named list を返す
prepare_profile <- function(profile, homo_to_any = FALSE, locus_order = NULL) {
  df <- data.frame(
    Locus = names(profile),
    allele1 = sapply(profile, function(x) if (length(x) >= 1) x[1] else ""),
    allele2 = sapply(profile, function(x) if (length(x) >= 2) x[2] else ""),
    stringsAsFactors = FALSE
  )
  
  # データフレーム版整形処理を呼び出し（freq_tableは不要なので省略）
  df_prep <- prepare_profile_df(df, freq_table = list(), homo_to_any = homo_to_any)
  
  # ローカス順補完（locus_orderに準拠）
  if (!is.null(locus_order)) {
    missing_loci <- setdiff(locus_order, df_prep$Locus)
    if (length(missing_loci) > 0) {
      df_missing <- data.frame(
        Locus = missing_loci,
        allele1 = "any",
        allele2 = "any",
        stringsAsFactors = FALSE
      )
      df_prep <- rbind(df_prep, df_missing)
    }
    df_prep$Locus <- factor(df_prep$Locus, levels = locus_order)
    df_prep <- df_prep[order(df_prep$Locus), ]
  }
  
  result <- setNames(
    lapply(seq_len(nrow(df_prep)), function(i) {
      c(df_prep$allele1[i], df_prep$allele2[i])
    }),
    df_prep$Locus
  )
  return(result)
}

# データフレーム版 prepare_profile（Shiny GUIなど表形式入力用）
# すべての列は小文字ベース（allele1, allele2）で処理、整形は常に適用
prepare_profile_df <- function(profile_df, freq_table, homo_to_any = FALSE) {
  # 列名の正規化（小文字で統一）
  colnames(profile_df) <- tolower(colnames(profile_df))
  
  # 必須列チェック
  if (!all(c("locus", "allele1", "allele2") %in% colnames(profile_df))) {
    stop("列名に 'locus', 'allele1', 'allele2' が含まれていません")
  }
  
  # locus列のみ大文字保持（互換性のため）
  profile_df <- profile_df[, c("locus", "allele1", "allele2")]
  colnames(profile_df) <- c("Locus", "allele1", "allele2")
  
  profile_df$allele1 <- as.character(profile_df$allele1)
  profile_df$allele2 <- as.character(profile_df$allele2)
  
  for (i in seq_len(nrow(profile_df))) {
    locus <- profile_df$Locus[i]
    alleles <- c(profile_df$allele1[i], profile_df$allele2[i])
    
    # 欠損・空白 → any
    alleles[is.na(alleles) | alleles == ""] <- "any"
    
    # freq_table に存在しないアレル → any
    valid_alleles <- freq_table[[locus]]
    alleles[!alleles %in% valid_alleles] <- "any"
    
    # ホモ型 → any 変換（オプション）
    if (homo_to_any && alleles[1] == alleles[2] && alleles[1] != "any") {
      alleles <- c("any", "any")
    }
    
    # アレルペア整形（昇順＋any右詰め）
    alleles <- std_allele_pair(alleles)
    
    profile_df$allele1[i] <- alleles[1]
    profile_df$allele2[i] <- alleles[2]
  }
  
  return(profile_df)
}

# クエリ用プロファイルの準備関数
# ファイルから読み込み、locus_order に従って補完とログ出力を行う
prepare_query <- function(query_file = "data/query_profile.csv",
                          locus_file = "data/locus_order.rds",
                          homo_to_any = TRUE) {
  if (!file.exists(locus_file)) {
    stop(paste("ローカス順ファイルが存在しません:", locus_file))
  }
  locus_order <- readRDS(locus_file)
  query_profile <- read_query_profile(query_file, locus_order, homo_to_any)
  
  # ログ出力（除外ローカス、補完ローカス）
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

# database_profile.csv を SampleID → Locus のネストリストに変換
convert_db_to_list <- function(df, locus_order) {
  df$Locus <- factor(df$Locus, levels = locus_order)
  df <- df[order(df$SampleID, df$Locus), ]
  split_profiles <- split(df, df$SampleID)
  
  profiles <- lapply(split_profiles, function(sample_df) {
    profile <- setNames(
      lapply(split(sample_df[, c("allele1", "allele2")], sample_df$Locus), unlist),
      as.character(unique(sample_df$Locus))
    )
    
    # 欠損ローカスを any,any で補完
    missing_loci <- setdiff(locus_order, names(profile))
    for (locus in missing_loci) {
      profile[[locus]] <- c("any", "any")
    }
    
    profile[locus_order]
  })
  
  return(profiles)
}

# DBファイルの読み込み〜整形までを行う統合処理（raw と整形済の両方を返す）
prepare_database <- function(db_file = "data/database_profile.csv",
                             locus_file = "data/locus_order.rds",
                             homo_to_any = FALSE) {
  if (!file.exists(locus_file)) stop(paste("Locus order file not found:", locus_file))
  if (!file.exists(db_file)) stop(paste("Database file not found:", db_file))
  
  locus_order <- readRDS(locus_file)
  df_raw <- read.csv(db_file, stringsAsFactors = FALSE)
  df_raw <- df_raw[df_raw$Locus %in% locus_order, ]
  db_profiles_raw <- convert_db_to_list(df_raw, locus_order)
  db_profiles <- lapply(db_profiles_raw, prepare_profile, homo_to_any = homo_to_any, locus_order = locus_order)
  
  return(list(raw = db_profiles_raw, prepared = db_profiles))
}

# アレルペア整形関数：昇順＋"any"右詰め（空欄・NAも any として処理）
std_allele_pair <- function(alleles) {
  alleles <- as.character(alleles)
  if (length(alleles) != 2) return(c("any", "any"))
  
  # 空欄・NA は any に変換
  alleles[is.na(alleles) | alleles == ""] <- "any"
  
  # any 以外の部分だけ昇順ソート
  non_any <- alleles[alleles != "any"]
  sorted <- suppressWarnings({
    if (length(non_any) == 2 && all(grepl("^[0-9]+(\\.[0-9]+)?$", non_any))) {
      as.character(sort(as.numeric(non_any)))
    } else if (length(non_any) == 1) {
      non_any
    } else {
      character(0)
    }
  })
  
  result <- c(sorted, rep("any", 2 - length(sorted)))
  return(result)
}
