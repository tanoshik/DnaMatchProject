run_match <- function(query_profile, db_profiles, top_n = 10) {
  # スコア算出
  results <- lapply(names(db_profiles), function(sid) {
    rec <- db_profiles[[sid]]
    res <- score_profile(query_profile, rec)
    list(SampleID = sid, Score = res$total, Details = res$scores)
  })
  
  # スコア順に並べ替え
  score_df <- data.frame(
    SampleID = sapply(results, `[[`, "SampleID"),
    Score = sapply(results, `[[`, "Score"),
    stringsAsFactors = FALSE
  )
  score_df <- score_df[order(-score_df$Score), ]
  
  # ログ生成：Top N件
  top_ids <- head(score_df$SampleID, top_n)
  # ログ生成：Top N件
  log <- do.call(rbind, lapply(top_ids, function(sid) {
    rec <- db_profiles[[sid]]
    scores <- score_profile(query_profile, rec)
    
    loci <- locus_order  # 明示的に順序固定
    
    q1 <- sapply(loci, function(locus) as.character(query_profile[[locus]][1, 1]))
    q2 <- sapply(loci, function(locus) as.character(query_profile[[locus]][1, 2]))
    r1 <- sapply(loci, function(locus) rec[[locus]][1])
    r2 <- sapply(loci, function(locus) rec[[locus]][2])
    ss <- scores$scores[loci]  # ローカス順に再並び替え
    
    data.frame(
      SampleID = rep(sid, length(loci)),
      Locus = loci,
      Query_Allele1 = q1,
      Query_Allele2 = q2,
      DB_Allele1 = r1,
      DB_Allele2 = r2,
      Score = as.vector(ss),
      stringsAsFactors = FALSE
    )
  }))
  
  # 出力保存
  write.csv(score_df, "output/match_scores.csv", row.names = FALSE, quote = FALSE)
  write.csv(log, "output/match_log.csv", row.names = FALSE, quote = FALSE)
  
  # 戻り値
  list(score_df = score_df, log = log)
}
