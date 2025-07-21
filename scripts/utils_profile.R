# 欠損補完とホモ型変換処理（リスト形式：内部処理用）
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

# data.frame版 prepare_profile（GUI用）
prepare_profile_df <- function(df, homo_to_any = FALSE) {
  df$Allele1 <- ifelse(is.na(df$Allele1) | df$Allele1 == "", "any", df$Allele1)
  df$Allele2 <- ifelse(is.na(df$Allele2) | df$Allele2 == "", "any", df$Allele2)
  if (homo_to_any) {
    df$Allele2 <- ifelse(df$Allele1 == df$Allele2 & df$Allele1 != "any", "any", df$Allele2)
  }
  return(df)
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
