# Initialize project data and RDS files

# scripts/init_project.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(tcltk)
})

# generate_rds_interactive() は定義されているが、ここでは実行しない
source("scripts/generate_rds_from_allele_count.R")

check_file <- function(path) {
  if (file.exists(path)) {
    cat(sprintf("[OK] File exists: %s\n", path))
  } else {
    cat(sprintf("[WARN] File missing: %s\n", path))
  }
}

check_file("data/freq_table.rds")
check_file("data/locus_order.rds")

cat("[INFO] If RDS files are missing, run generate_rds_interactive() manually.\n")
