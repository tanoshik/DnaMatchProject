# 出現頻度計算ユーティリティ

# ローカス単位の出現頻度計算（any対応）
# 出現頻度計算（long形式 freq_table対応）
calc_freq_loci <- function(locus, allele1, allele2, freq_table) {
  subset_table <- freq_table[freq_table$Locus == locus, ]
  total_count <- sum(subset_table$Count, na.rm = TRUE)
  if (total_count == 0) return(1)
  
  get_p <- function(allele) {
    count <- subset_table$Count[subset_table$Allele == allele]
    if (length(count) == 0 || is.na(count)) return(0)
    return(count / total_count)
  }
  
  if (allele1 == "any" && allele2 == "any") return(1)
  
  if (allele1 == "any" || allele2 == "any") {
    allele <- ifelse(allele1 != "any", allele1, allele2)
    p <- get_p(allele)
    return(2 * p - p^2)
  }
  
  p1 <- get_p(allele1)
  p2 <- get_p(allele2)
  if (p1 == 0 || p2 == 0) return(1)
  
  if (allele1 == allele2) {
    return(p1^2)
  } else {
    return(2 * p1 * p2)
  }
}

# データフレーム形式のプロファイルに対してFreq列を追加
calc_freq_loci_df <- function(profile_df, freq_table) {
  profile_df$Freq <- mapply(
    function(locus, a1, a2) {
      calc_freq_loci(locus, a1, a2, freq_table)
    },
    profile_df$Locus,
    profile_df$allele1,
    profile_df$allele2
  )
  return(profile_df)
}

# プロファイル全体の出現頻度（ローカスごとの積）
calc_total_freq <- function(profile_df, freq_table) {
  loci <- profile_df$Locus
  allele1 <- profile_df$allele1
  allele2 <- profile_df$allele2
  
  total <- 1
  
  for (i in seq_along(loci)) {
    freq <- calc_freq_loci(loci[i], allele1[i], allele2[i], freq_table)
    total <- total * freq
  }
  
  return(total)
}
