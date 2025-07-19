# DnaMatchProject

This project provides a lightweight DNA profile matching system based on STR loci.  
It enables profile comparison, scoring, and ranking based on allele overlap, using custom matching logic written in R.

The system is modular, testable, and structured for forensic or simulation purposes.

---

## Directory Structure

The project is organized as follows:

```
DnaMatchProject/
├── data/              # Input files: allele counts, profiles, locus order
├── output/            # Matching results
├── scripts/           # Core R scripts for processing
├── test/              # Test inputs and scripts (excluded from Git)
├── main.R             # Entry point to execute the matching system
└── DnaMatchProject.Rproj  # RStudio project file
```

---

## Setup

### Requirements

- R version: 4.5.0 or later
- RStudio (optional but recommended)
- Required R packages:
  - `readr`
  - `dplyr`

### Getting Started

```r
# Clone the repo and open in RStudio
# Then install packages:
install.packages(c("readr", "dplyr"))

# Optionally regenerate RDS files
source("scripts/generate_rds_from_allele_count.R")

# Execute main script
source("main.R")
```

Default settings are based on the GlobalFiler kit and allele frequencies from a Japanese population:

> Fujii K, Watahiki H, Mita Y, et al.  
> *Allele frequencies for 21 autosomal short tandem repeat loci obtained using GlobalFiler in a sample of 1501 individuals from the Japanese population*.  
> Legal Medicine (Tokyo). 2015 Sep;17(5):306–308.  
> doi: [10.1016/j.legalmed.2015.08.007](https://doi.org/10.1016/j.legalmed.2015.08.007)

---

## Usage

1. Prepare `data/database_profile.csv` and `data/query_profile.csv`
2. Configure matching behavior in `main.R`:

```r
homo_db <- TRUE  # treat homozygous in DB as [value, any]
homo_q  <- TRUE  # treat homozygous in query the same way
```

3. Run:

```r
source("main.R")
```

Results are saved in:

- `output/match_scores.csv`
- `output/match_log.csv`
