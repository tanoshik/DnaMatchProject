# Entry point for the DnaMatchProject

# Load libraries and scripts
# source("scripts/init_project.R")
source("scripts/utils_profile.R")
source("scripts/io_profiles.R")
source("scripts/scoring.R")
source("scripts/matcher.R")

# Settings
query_file <- "data/query_profile.csv"
db_file <- "data/database_profile.csv"
locus_order <- readRDS("data/locus_order.rds")
homo_q <- TRUE   # Convert homozygous query alleles to "any"
homo_r <- FALSE  # Do not convert homozygous db alleles
top_n <- 10      # Number of top hits to display

# Load and prepare profiles
query_profile <- read_query_profile(query_file, locus_order, homo_to_any = homo_q)
db_profiles <- read_db_profiles(db_file, locus_order, homo_to_any = homo_r)

# Run matching
results <- run_match(query_profile, db_profiles, top_n)

# Output results
cat("[INFO] Project initialized successfully.\n")
cat("[OK] Matching completed.\n")
cat(sprintf("Top hit: %s with score %d\n",
            results$score_df$SampleID[1],
            results$score_df$Score[1]))
