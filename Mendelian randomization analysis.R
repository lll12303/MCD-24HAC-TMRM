# ============================================================
# STEP 1: Load libraries
# ============================================================

library(TwoSampleMR)
library(data.table)
library(dplyr)
library(ggplot2)

# ============================================================
# STEP 2: File paths
# ============================================================

exposure_files <- list(
  exp1 = "data/exposure1.txt",
  exp2 = "data/exposure2.txt"
)

outcome_files <- list(
  out1 = "data/outcome1.gz",
  out2 = "data/outcome2.gz"
)

# ============================================================
# STEP 3: MR function
# ============================================================

run_mr_pair <- function(exposure_file, outcome_file,
                        exposure_name, outcome_name) {
  
  message("Running: ", exposure_name, " -> ", outcome_name)
  
  # ---------------------------
  # Exposure
  # ---------------------------
  exp <- tryCatch(fread(exposure_file), error = function(e) return(NULL))
  if (is.null(exp)) return(NULL)
  
  exp$PVAL <- 10^(-exp$LOG10P)
  
  exp <- exp %>%
    arrange(PVAL) %>%
    slice_head(n = 50000)   # safer cap
  
  if (!"A1FREQ" %in% names(exp)) {
    exp$A1FREQ <- 0.5
  }
  
  exp <- exp %>%
    transmute(
      CHR = as.character(CHROM),
      POS = GENPOS,
      effect_allele = ALLELE1,
      other_allele  = ALLELE0,
      beta = BETA,
      se   = SE,
      eaf  = A1FREQ,
      pval = PVAL
    )
  
  # ---------------------------
  # Outcome
  # ---------------------------
  out <- tryCatch(fread(outcome_file), error = function(e) return(NULL))
  if (is.null(out)) return(NULL)
  
  out <- out %>%
    transmute(
      CHR = as.character(`#chrom`),
      POS = pos,
      effect_allele = alt,
      other_allele  = ref,
      beta = beta,
      se   = sebeta,
      pval = pval
    )
  
  # ---------------------------
  # Merge
  # ---------------------------
  setDT(exp); setDT(out)
  setkey(exp, CHR, POS)
  setkey(out, CHR, POS)
  
  dat <- merge(exp, out, by = c("CHR", "POS"),
               suffixes = c("_exposure", "_outcome"))
  
  if (nrow(dat) < 10) return(NULL)
  
  dat$SNP <- paste0("chr", dat$CHR, "_", dat$POS)
  
  dat <- as.data.frame(dat)
  
  # ---------------------------
  # Format
  # ---------------------------
  exposure_dat <- format_data(
    dat,
    type = "exposure",
    snp_col = "SNP",
    beta_col = "beta_exposure",
    se_col = "se_exposure",
    pval_col = "pval_exposure",
    effect_allele_col = "effect_allele_exposure",
    other_allele_col = "other_allele_exposure",
    eaf_col = "eaf"
  )
  
  outcome_dat <- format_data(
    dat,
    type = "outcome",
    snp_col = "SNP",
    beta_col = "beta_outcome",
    se_col = "se_outcome",
    pval_col = "pval_outcome",
    effect_allele_col = "effect_allele_outcome",
    other_allele_col = "other_allele_outcome"
  )
  
  # ---------------------------
  # Harmonise + MR
  # ---------------------------
  harm <- harmonise_data(exposure_dat, outcome_dat)
  
  if (nrow(harm) < 10) return(NULL)
  
  res <- mr(harm)
  
  res <- res %>%
    mutate(
      exposure = exposure_name,
      outcome  = outcome_name,
      nsnp     = nrow(harm),
      OR       = exp(b),
      OR_lci   = exp(b - 1.96 * se),
      OR_uci   = exp(b + 1.96 * se)
    )
  
  return(res)
}

# ============================================================
# STEP 4: Batch run
# ============================================================

results <- list()

for (e in names(exposure_files)) {
  for (o in names(outcome_files)) {
    
    res <- tryCatch(
      run_mr_pair(
        exposure_files[[e]],
        outcome_files[[o]],
        e, o
      ),
      error = function(e) NULL
    )
    
    if (!is.null(res)) {
      results[[paste(e, o, sep = "_")]] <- res
    }
  }
}

mr_all <- bind_rows(results)

# ============================================================
# STEP 5: Format results
# ============================================================

plot_data <- mr_all %>%
  filter(method == "Inverse variance weighted") %>%
  mutate(
    sig = case_when(
      pval < 0.001 ~ "***",
      pval < 0.01  ~ "**",
      pval < 0.05  ~ "*",
      TRUE ~ ""
    )
  )