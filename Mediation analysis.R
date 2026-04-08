# ============================================================
# STEP 0: Load libraries and data
# ============================================================

library(tidyverse)
library(broom)
library(survival)
library(lavaan)

# ---------------------------
# Base data (covariates + exposure)
# ---------------------------
base_data <- read.csv("all_static_features.csv", check.names = FALSE) %>%
  select(
    `Participant ID`,
    age,
    Sex,
    `Body mass index (BMI)`,
    `Time difference`
  ) %>%
  inner_join(
    read.csv("exposure.csv", check.names = FALSE),
    by = "Participant ID"
  )

# ---------------------------
# Mediators (generic naming)
# ---------------------------
dataset1 <- read.csv("dataset1.csv", check.names = FALSE)
dataset2 <- read.csv("dataset2.csv", check.names = FALSE)
dataset3 <- read.csv("dataset3.csv", check.names = FALSE)

mediator_sets <- list(
  Dataset1 = dataset1,
  Dataset2 = dataset2,
  Dataset3 = dataset3
)

# ---------------------------
# Outcome
# ---------------------------
outcome_data <- read.csv("outcome.csv", check.names = FALSE)

# ---------------------------
# Survival data
# ---------------------------
disease_list <- list(
  Disease1 = read.csv("disease1_survival.csv"),
  Disease2 = read.csv("disease2_survival.csv"),
  Disease3 = read.csv("disease3_survival.csv"),
  Disease4 = read.csv("disease4_survival.csv")
)

# ============================================================
# STEP 1: Exposure -> Mediator
# ============================================================

exposures <- setdiff(names(base_data),
                     c("Participant ID", "age", "Sex",
                       "Body mass index (BMI)", "Time difference"))

covars <- c("age", "Sex", "Time difference")

run_step1 <- function(mediator_data, dataset_name) {
  
  full_data <- mediator_data %>%
    inner_join(base_data, by = "Participant ID")
  
  mediator_cols <- setdiff(names(mediator_data), "Participant ID")
  
  res_list <- list()
  
  for (exp_col in exposures) {
    for (m in mediator_cols) {
      
      df <- full_data %>%
        select(all_of(c(m, exp_col, covars))) %>%
        drop_na()
      
      if (nrow(df) < 300) next
      
      fit <- tryCatch(
        lm(as.formula(
          paste0("`", m, "` ~ scale(`", exp_col, "`) + ",
                 paste0("`", covars, "`", collapse = " + "))
        ), data = df),
        error = function(e) NULL
      )
      
      if (is.null(fit)) next
      
      res <- tidy(fit) %>%
        filter(term == paste0("scale(`", exp_col, "`)"))
      
      if (nrow(res) == 1) {
        res_list[[length(res_list) + 1]] <- tibble(
          dataset = dataset_name,
          exposure = exp_col,
          mediator = m,
          beta = res$estimate,
          p = res$p.value
        )
      }
    }
  }
  
  bind_rows(res_list) %>%
    mutate(FDR = p.adjust(p, "fdr"))
}

step1_all <- bind_rows(
  lapply(names(mediator_sets), function(n) {
    run_step1(mediator_sets[[n]], n)
  })
)

# ============================================================
# STEP 2: Mediator -> Disease (Cox)
# ============================================================

run_step2 <- function(mediator_data, survival_data, disease_name, dataset_name) {
  
  survival_data <- survival_data %>%
    rename(surv_time = time_years,
           surv_event = event)
  
  full_data <- mediator_data %>%
    inner_join(base_data, by = "Participant ID") %>%
    inner_join(survival_data, by = "Participant ID")
  
  mediator_cols <- setdiff(names(mediator_data), "Participant ID")
  
  res_list <- list()
  
  for (m in mediator_cols) {
    
    df <- full_data %>%
      select(all_of(c(m, covars, "surv_time", "surv_event"))) %>%
      drop_na()
    
    if (sum(df$surv_event) < 30) next
    
    fit <- tryCatch(
      coxph(as.formula(
        paste0("Surv(surv_time, surv_event) ~ scale(`", m, "`) + ",
               paste0("`", covars, "`", collapse = " + "))
      ), data = df),
      error = function(e) NULL
    )
    
    if (is.null(fit)) next
    
    res <- tidy(fit, exponentiate = TRUE, conf.int = TRUE)
    
    res_list[[length(res_list) + 1]] <- tibble(
      disease = disease_name,
      dataset = dataset_name,
      mediator = m,
      HR = res$estimate[1],
      p = res$p.value[1]
    )
  }
  
  bind_rows(res_list) %>%
    mutate(FDR = p.adjust(p, "fdr"))
}

step2_all <- bind_rows(
  lapply(names(disease_list), function(disease_name) {
    
    bind_rows(
      lapply(names(mediator_sets), function(dataset_name) {
        
        run_step2(
          mediator_sets[[dataset_name]],
          disease_list[[disease_name]],
          disease_name,
          dataset_name
        )
      })
    )
  })
)

# ============================================================
# STEP 3: Intersection
# ============================================================

get_intersection <- function(step1, step2, disease_name, dataset_name) {
  
  s1 <- step1 %>% filter(FDR < 0.05)
  s2 <- step2 %>% filter(FDR < 0.05)
  
  if (nrow(s1) == 0 || nrow(s2) == 0) return(NULL)
  
  inner_join(s1, s2, by = "mediator") %>%
    mutate(disease = disease_name, dataset = dataset_name)
}

all_intersections <- bind_rows(
  lapply(names(disease_list), function(disease_name) {
    
    bind_rows(
      lapply(names(mediator_sets), function(dataset_name) {
        
        s2 <- step2_all %>%
          filter(disease == disease_name,
                 dataset == dataset_name)
        
        get_intersection(step1_all, s2,
                         disease_name, dataset_name)
      })
    )
  })
)

# ============================================================
# STEP 4: SEM (Top mediators)
# ============================================================

master_data <- base_data %>%
  inner_join(outcome_data, by = "Participant ID")

get_top <- function(disease_name, dataset_name, exposure) {
  
  all_intersections %>%
    filter(
      disease == disease_name,
      dataset == dataset_name,
      exposure == exposure
    ) %>%
    arrange(FDR.y) %>%
    slice_head(n = 5) %>%
    pull(mediator)
}

run_sem <- function(data, mediators, exposure, disease_name) {
  
  if (length(mediators) == 0) return(NULL)
  
  names(data) <- make.names(names(data), unique = TRUE)
  
  mediators <- make.names(mediators)
  exposure  <- make.names(exposure)
  disease   <- make.names(disease_name)
  covs      <- make.names(covars)
  
  df <- data %>%
    select(all_of(c(mediators, exposure, disease, covs))) %>%
    drop_na()
  
  if (nrow(df) < 1000) return(NULL)
  
  df[[exposure]] <- scale(df[[exposure]])
  df[mediators]  <- lapply(df[mediators], scale)
  
  model <- paste(
    paste0(mediators, " ~ ", exposure, collapse = "\n"),
    paste0(disease, " ~ ",
           paste(mediators, collapse = " + "),
           " + ", exposure),
    sep = "\n"
  )
  
  fit <- tryCatch(
    sem(model, data = df,
        ordered = disease,
        estimator = "WLSMV"),
    error = function(e) NULL
  )
  
  if (is.null(fit)) return(NULL)
  
  parameterEstimates(fit, standardized = TRUE)
}

sem_results <- bind_rows(
  lapply(names(disease_list), function(disease_name) {
    
    bind_rows(
      lapply(names(mediator_sets), function(dataset_name) {
        
        bind_rows(
          lapply(exposures, function(exp) {
            
            meds <- get_top(disease_name, dataset_name, exp)
            
            if (length(meds) == 0) return(NULL)
            
            data_tmp <- master_data %>%
              inner_join(mediator_sets[[dataset_name]],
                         by = "Participant ID")
            
            run_sem(data_tmp, meds, exp, disease_name)
          })
        )
      })
    )
  })
)