# =========================================================
# 1. Libraries
# =========================================================
library(tidyverse)
library(rms)
library(survival)
library(patchwork)
library(numDeriv)

# =========================================================
# 2. Load and preprocess data
# =========================================================
data <- read.csv("SCR.csv")

dd <- datadist(data)
options(datadist = "dd")
# =========================================================
# 3. Generic RCS plotting function
# =========================================================
create_rcs_plot <- function(fit, var, label, color, mode = "intersection") {
  
  mode <- match.arg(mode, c("intersection", "minimum", "maximum", "inflection"))
  
  # Generate prediction grid automatically
  x_vals <- seq(
    quantile(data[[var]], 0.01, na.rm = TRUE),
    quantile(data[[var]], 0.99, na.rm = TRUE),
    length.out = 200
  )
  
  pred <- Predict(fit, x_vals, fun = exp)
  df <- as.data.frame(pred)
  
  colnames(df)[1] <- "x"
  
  # ---------------- Key point detection ----------------
  if (mode == "intersection") {
    idx <- which.min(abs(df$yhat - 1))
  } else if (mode == "minimum") {
    idx <- which.min(df$yhat)
  } else if (mode == "maximum") {
    idx <- which.max(df$yhat)
  } else if (mode == "inflection") {
    
    d1 <- grad(function(z){
      approx(df$x, df$yhat, xout=z)$y
    }, df$x)
    
    d2 <- grad(function(z){
      approx(df$x, d1, xout=z)$y
    }, df$x)
    
    idx <- which.min(abs(d2))
  }
  
  px <- df$x[idx]
  py <- df$yhat[idx]
  
  label_text <- sprintf("%.1f min\nHR=%.2f", px, py)
  
  # ---------------- Plot ----------------
  ggplot(df) +
    geom_ribbon(aes(x = x, ymin = lower, ymax = upper),
                fill = color, alpha = 0.25) +
    geom_line(aes(x = x, y = yhat),
              color = color, size = 1.3) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    
    geom_point(aes(x = px, y = py),
               color = "red", size = 3) +
    
    geom_text(aes(x = px, y = py, label = label_text),
              color = "red", vjust = -1.2, size = 5, fontface = "bold") +
    
    labs(
      x = paste0(label, " (minutes/day)"),
      y = "Hazard ratio"
    ) +
    
    theme_bw(base_size = 16) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "black")
    )
}

# =========================================================
# 4. Model fitting
# =========================================================
fit_mvpa <- cph(Surv(time_diff_days, disease_binary) ~ 
                  rcs(MVPA_min, 3) + age + Sex,
                data = data, x = TRUE, y = TRUE)

fit_light <- cph(Surv(time_diff_days, disease_binary) ~ 
                   rcs(Light_min, 3) + age + Sex,
                 data = data, x = TRUE, y = TRUE)

fit_sedentary <- cph(Surv(time_diff_days, disease_binary) ~ 
                       rcs(Sedentary_min, 3) + age + Sex,
                     data = data, x = TRUE, y = TRUE)

fit_sleep <- cph(Surv(time_diff_days, disease_binary) ~ 
                   rcs(Sleep_min, 3) + age + Sex,
                 data = data, x = TRUE, y = TRUE)

# =========================================================
# 5. Generate plots
# =========================================================
p_mvpa <- create_rcs_plot(fit_mvpa, "MVPA_min", "MVPA", "#E31A1C", "intersection")

p_light <- create_rcs_plot(fit_light, "Light_min", "LPA", "#1F78B4", "intersection")

p_sedentary <- create_rcs_plot(fit_sedentary, "Sedentary_min", "SB", "#33A02C", "intersection")

p_sleep <- create_rcs_plot(fit_sleep, "Sleep_min", "Sleep", "#FFD92F", "minimum")
