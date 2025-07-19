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

# プロファイル単位のスコア
score_profile <- function(query_profile, record_profile) {
  cat("Length of query_profile:", length(query_profile), "\n")
  cat("Length of record_profile:", length(record_profile), "\n")
  cat("Names of query_profile:\n"); print(names(query_profile))
  cat("Names of record_profile:\n"); print(names(record_profile))
  scores <- mapply(score_locus, query_profile, record_profile)
  names(scores) <- names(query_profile)  # 明示的にローカス名を付与
  c(scores, total = sum(scores))         # ← 1つのnamed numeric vectorとして返す
}
