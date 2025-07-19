# query読み込み（補完付き）
read_query_profile <- function(file_path, locus_order, homo_to_any = FALSE) {
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  df <- df[df$Locus %in% locus_order, ]
  df$Locus <- as.character(df$Locus)
  
  profile <- list()
  for (locus in locus_order) {
    row <- df[df$Locus == locus, ]
    if (nrow(row) == 0) {
      alleles <- c("any", "any")
    } else {
      a1 <- ifelse(is.na(row$allele1) || row$allele1 == "", "any", row$allele1)
      a2 <- ifelse(is.na(row$allele2) || row$allele2 == "", "any", row$allele2)
      alleles <- c(a1, a2)
    }
    profile[[locus]] <- alleles
  }
  profile <- prepare_profile(profile, homo_to_any)
  names(profile) <- locus_order
  return(profile)
}

# db読み込み（補完付き）
read_db_profiles <- function(file_path, locus_order, homo_to_any = FALSE) {
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  sample_ids <- unique(df$SampleID)
  profiles <- lapply(sample_ids, function(sid) {
    sub_df <- df[df$SampleID == sid & df$Locus %in% locus_order, ]
    sub_df$Locus <- factor(sub_df$Locus, levels = locus_order)
    sub_df <- sub_df[order(sub_df$Locus), ]
    
    raw_profile <- setNames(
      lapply(split(sub_df[, c("allele1", "allele2")], sub_df$Locus), unlist),
      as.character(unique(sub_df$Locus))
    )
    
    missing_loci <- setdiff(locus_order, names(raw_profile))
    for (locus in missing_loci) {
      raw_profile[[locus]] <- c(NA, NA)
    }
    raw_profile <- raw_profile[locus_order]
    prepare_profile(raw_profile, homo_to_any)
  })
  names(profiles) <- sample_ids
  profiles
}

# DEBUG: query生データ vs prepared比較表示
print_raw_vs_prepared_query <- function(file_path, prepared, locus_order) {
  cat("=== Query Profile: Raw vs Prepared ===\n")
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  for (locus in locus_order) {
    row <- df[df$Locus == locus, ]
    raw <- if (nrow(row) == 0) c("NA", "NA") else c(row$allele1, row$allele2)
    prep <- prepared[[locus]]
    cat(sprintf("%-11s | query_raw: %-10s -> prepared: %-10s\n",
                locus,
                paste(raw, collapse = ","),
                paste(prep, collapse = ",")))
  }
}

# DEBUG: db生データ vs prepared比較表示
print_raw_vs_prepared_db <- function(file_path, prepared, locus_order) {
  cat("\n=== DB Profile: Raw vs Prepared ===\n")
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  df <- df[df$SampleID == unique(df$SampleID)[1], ]
  for (locus in locus_order) {
    row <- df[df$Locus == locus, ]
    raw <- if (nrow(row) == 0) c("NA", "NA") else c(row$allele1, row$allele2)
    prep <- prepared[[locus]]
    cat(sprintf("%-11s | db_raw:    %-10s -> prepared: %-10s\n",
                locus,
                paste(raw, collapse = ","),
                paste(prep, collapse = ",")))
  }
}
