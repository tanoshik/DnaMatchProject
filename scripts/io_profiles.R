# query読み込み（補完付き）
read_query_profile <- function(file_path, locus_order, homo_to_any = FALSE) {
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  df <- df[df$Locus %in% locus_order, ]
  df$Locus <- as.character(df$Locus)
  df <- df[match(df$Locus, locus_order, nomatch = 0) > 0, ]

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
  prepare_profile(profile, homo_to_any)
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
