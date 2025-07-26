# ローカス単位のスコア
score_locus <- function(query, record) {
  if (all(query == "any") || all(record == "any")) return(2)
  q_real <- query[query != "any"]
  r_real <- record[record != "any"]
  matched <- rep(FALSE, length(r_real))
  score <- 0
  for (q in q_real) {
    for (i in seq_along(r_real)) {
      if (!matched[i] && q == r_real[i]) {
        matched[i] <- TRUE
        score <- score + 1
        break
      }
    }
  }
  min(score + sum(query == "any") + sum(record == "any"), 2)
}

# プロファイル単位のスコア＋詳細ログ
score_profile <- function(query_profile, record_profile) {
  # DEBUG: 構造確認
  cat(">> DEBUG: score_profile() called\n")
  cat(" - query_profile type: ", typeof(query_profile), "\n")
  cat(" - record_profile type: ", typeof(record_profile), "\n")
  cat(" - Names of query_profile:\n"); print(names(query_profile))
  cat(" - First query locus:\n"); print(query_profile[[1]])
  cat(" - First record locus:\n"); print(record_profile[[1]])
  
  loci <- names(query_profile)
  total <- 0
  per_locus_scores <- numeric(length(loci))
  names(per_locus_scores) <- loci
  
  log_df <- data.frame(
    Locus = character(),
    query_allele1 = character(),
    query_allele2 = character(),
    record_allele1 = character(),
    record_allele2 = character(),
    score = integer(),
    stringsAsFactors = FALSE
  )
  
  for (locus in loci) {
    q <- query_profile[[locus]]
    r <- record_profile[[locus]]
    
    # 🛠️ 安全のため長さ2を保証
    if (is.null(q) || length(q) != 2) q <- c("any", "any")
    if (is.null(r) || length(r) != 2) r <- c("any", "any")
    
    s <- score_locus(q, r)
    per_locus_scores[locus] <- s
    total <- total + s
    
    log_df <- rbind(log_df, data.frame(
      Locus = locus,
      query_allele1 = q[1], query_allele2 = q[2],
      record_allele1 = r[1], record_allele2 = r[2],
      score = s,
      stringsAsFactors = FALSE
    ))
  }
  
  list(
    scores = per_locus_scores,
    total = total,
    log = log_df
  )
}
