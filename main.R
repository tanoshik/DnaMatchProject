source("scripts/utils_profile.R")
source("scripts/io_profiles.R")
source("scripts/scoring.R")
source("scripts/matcher.R")

query_file <- "data/query_profile.csv"
db_file <- "data/database_profile.csv"
homo_q <- TRUE
homo_r <- FALSE
top_n <- 10

# 読み込み
query_profile <- prepare_query(query_file, "data/locus_order.rds", homo_q)
db_profiles <- prepare_database(db_file, "data/locus_order.rds", homo_r)
locus_order <- readRDS("data/locus_order.rds")
db_profiles <- read_db_profiles(db_file, locus_order, homo_r)
results <- run_match(query_profile, db_profiles, top_n)
print(head(results$score_df, top_n))
