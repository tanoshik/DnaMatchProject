# ----------------------------
# Load core libraries and scripts
# ----------------------------
# source("scripts/init_project.R")
source("scripts/utils_profile.R")
source("scripts/io_profiles.R")
source("scripts/scoring.R")
source("scripts/matcher.R")

# ----------------------------
# Ensure required RDS files exist
# ----------------------------
if (!file.exists("data/freq_table.rds") || !file.exists("data/locus_order.rds")) {
  message("RDS files not found. Generating from Allele-count_GF_Japanese.csv ...")
  source("scripts/generate_rds_from_allele_count.R")
  generate_rds("data/Allele-count_GF_Japanese.csv")
}
# TODO: Support switching allele count sources (e.g., for non-GlobalFiler kits)

# ----------------------------
# Settings and profile loading
# ----------------------------
query_file <- "data/query_profile.csv"
db_file <- "data/database_profile.csv"
locus_order <- readRDS("data/locus_order.rds")
homo_q <- TRUE   # Convert homozygous query alleles to "any"
homo_r <- FALSE  # Do not convert homozygous db alleles
top_n <- 10      # Number of top hits to display

query_profile <- read_query_profile(query_file, locus_order, homo_to_any = homo_q)
db_profiles   <- read_db_profiles(db_file, locus_order, homo_to_any = homo_r)

# ----------------------------
# Run matching and show result
# ----------------------------
results <- run_match(query_profile, db_profiles, top_n)

cat("[INFO] Project initialized successfully.\n")
cat("[OK] Matching completed.\n")
cat(sprintf("Top hit: %s with score %d\n",
            results$score_df$SampleID[1],
            results$score_df$Score[1]))
