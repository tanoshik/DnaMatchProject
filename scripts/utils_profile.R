# アレルペア整形関数：昇順＋"any"右詰め（空欄・NAも any として処理）
std_allele_pair <- function(alleles) {
  alleles <- as.character(alleles)
  if (length(alleles) != 2) return(c("any", "any"))
  alleles[is.na(alleles) | alleles == ""] <- "any"
  
  non_any <- alleles[alleles != "any"]
  sorted <- suppressWarnings({
    if (length(non_any) == 2 && all(grepl("^[0-9]+(\\.[0-9]+)?$", non_any))) {
      as.character(sort(as.numeric(non_any)))
    } else {
      non_any
    }
  })
  result <- c(sorted, rep("any", 2 - length(sorted)))
  return(result)
}

# data.frame版 prepare_profile（GUI用）
# 空欄補完・ホモ型→any変換・std_allele_pair() による整形をすべて行う
prepare_profile_df <- function(df, homo_to_any = FALSE) {
  df$allele1 <- ifelse(is.na(df$allele1) | df$allele1 == "", "any", df$allele1)
  df$allele2 <- ifelse(is.na(df$allele2) | df$allele2 == "", "any", df$allele2)
  
  # ホモ型を any に変換
  if (homo_to_any) {
    df$allele2 <- ifelse(df$allele1 == df$allele2 & df$allele1 != "any", "any", df$allele2)
  }
  
  # アレルペアの整形（std_allele_pair）を適用
  std_pairs <- t(apply(df[, c("allele1", "allele2")], 1, std_allele_pair))
  df$allele1 <- std_pairs[, 1]
  df$allele2 <- std_pairs[, 2]
  
  return(df)
}

# リスト形式プロファイルの整形：prepare_profile_df を内部で使用
prepare_profile <- function(profile, homo_to_any = FALSE, locus_order = NULL) {
  df <- data.frame(
    Locus = names(profile),
    allele1 = sapply(profile, function(x) if (length(x) >= 1) x[1] else ""),
    allele2 = sapply(profile, function(x) if (length(x) >= 2) x[2] else ""),
    stringsAsFactors = FALSE
  )
  
  df_prep <- prepare_profile_df(df, homo_to_any)
  
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

# 検索クエリの準備
prepare_query <- function(query_file = "data/query_profile.csv",
                          locus_file = "data/locus_order.rds",
                          homo_to_any = TRUE) {
  if (!file.exists(locus_file)) {
    stop(paste("ローカス順ファイルが存在しません:", locus_file))
  }
  locus_order <- readRDS(locus_file)
  query_profile <- read_query_profile(query_file, locus_order, homo_to_any)
  
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

# Convert raw dataframe to nested list by SampleID and Locus
convert_db_to_list <- function(df, locus_order) {
  df$Locus <- factor(df$Locus, levels = locus_order)
  df <- df[order(df$SampleID, df$Locus), ]
  split_profiles <- split(df, df$SampleID)
  
  profiles <- lapply(split_profiles, function(sample_df) {
    profile <- setNames(
      lapply(split(sample_df[, c("allele1", "allele2")], sample_df$Locus), unlist),
      as.character(unique(sample_df$Locus))
    )
    missing_loci <- setdiff(locus_order, names(profile))
    for (locus in missing_loci) {
      profile[[locus]] <- c("any", "any")
    }
    profile[locus_order]
  })
  
  return(profiles)
}

# Prepare database: read, log, fill, convert
prepare_database <- function(db_file = "data/database_profile.csv",
                             locus_file = "data/locus_order.rds",
                             homo_to_any = FALSE) {
  if (!file.exists(locus_file)) stop(paste("Locus order file not found:", locus_file))
  if (!file.exists(db_file)) stop(paste("Database file not found:", db_file))
  
  locus_order <- readRDS(locus_file)
  df_raw <- read.csv(db_file, stringsAsFactors = FALSE)
  df_raw <- df_raw[df_raw$Locus %in% locus_order, ]
  db_profiles_raw <- convert_db_to_list(df_raw, locus_order)
  db_profiles <- lapply(db_profiles_raw, prepare_profile, homo_to_any = homo_to_any)
  
  return(list(raw = db_profiles_raw, prepared = db_profiles))
}
