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
  scores <- mapply(score_locus, query_profile, record_profile)
  total <- sum(scores)
  list(scores = scores, total = total)
}
