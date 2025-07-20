run_match <- function(query_profile, db_profiles, top_n = 10) {
  results <- lapply(names(db_profiles), function(sid) {
    rec <- db_profiles[[sid]]
    res <- score_profile(query_profile, rec)
    list(SampleID = sid, Score = res$total, Details = res$scores)
  })
  
  # Topスコア順に並べ替え
  score_df <- data.frame(
    SampleID = sapply(results, `[[`, "SampleID"),
    Score = sapply(results, `[[`, "Score"),
    stringsAsFactors = FALSE
  )
  score_df <- score_df[order(-score_df$Score), ]
  
  # ログ生成：Top N件
  top_ids <- score_df$SampleID[seq_len(min(top_n, nrow(score_df)))]
  log <- do.call(rbind, lapply(top_ids, function(sid) {
    rec <- db_profiles[[sid]]
    scores <- score_profile(query_profile, rec)
    loci <- names(query_profile)
    data.frame(
      SampleID = sid,
      Locus = loci,
      Query = sapply(query_profile, paste, collapse = ","),
      Record = sapply(rec, paste, collapse = ","),
      Score = scores$scores,
      stringsAsFactors = FALSE
    )
  }))
  
  # CSV出力
  write.csv(score_df, "output/match_scores.csv", row.names = FALSE, quote = FALSE)
  write.csv(log, "output/match_log.csv", row.names = FALSE, quote = FALSE)
  
  # 戻り値
  list(score_df = score_df, log = log)
}
