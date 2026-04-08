# ============================================================
# 0️⃣ Libraries
# ============================================================
library(data.table)
library(ggplot2)
library(patchwork)
library(cowplot)
library(scales)

# ============================================================
# 1️⃣ Read data
# ============================================================
ukb_pa <- fread("best_recommendation_data.csv",check.names = FALSE)

feat_cols <- paste0(0:191)
disease_list <- unique(ukb_pa$Disease)

dt <- as.data.table(ukb_pa)

stopifnot(all(feat_cols %in% names(dt)))

# ============================================================
# 2️⃣ Canonical mapping
# ============================================================
behavior_map <- data.table(index = 0:191)

behavior_map[, behavior := factor(
  index %% 4,
  levels = 0:3,
  labels = c("SB", "LPA", "MVPA", "Sleep")
)]

behavior_map[, day_type := factor(
  ifelse(index < 96, "Weekday", "Weekend"),
  levels = c("Weekday", "Weekend")
)]

behavior_map[, hour := factor(
  floor((index %% 96) / 4),
  levels = 0:23,
  ordered = TRUE
)]

setkey(behavior_map, index)

# ============================================================
# 3️⃣ Population profile
# ============================================================
compute_population_profile <- function(df_wide) {
  
  dt_wide <- as.data.table(df_wide)
  
  long <- melt(
    dt_wide,
    measure.vars = feat_cols,
    variable.name = "index",
    value.name = "value",
    variable.factor = FALSE
  )
  
  long[, index := as.integer(index)]
  
  long <- behavior_map[long, on = "index", nomatch = 0]
  
  prof <- long[
    ,
    .(value = mean(value, na.rm = TRUE)),
    by = .(day_type, hour, behavior)
  ]
  
  setorder(prof, day_type, hour, behavior)
  
  # sanity check
  chk <- prof[, .N, by = .(day_type, hour)]
  stopifnot(all(chk$N == 4))
  
  return(prof)
}

# ============================================================
# 4️⃣ Low-risk target profile
# ============================================================
compute_low_profile <- function(df, lo_q = 0.05) {
  
  lo_cut <- quantile(df$Predicted_Risk, lo_q, na.rm = TRUE)
  
  df_low <- df[Predicted_Risk <= lo_cut]
  
  long <- melt(
    df_low,
    measure.vars = feat_cols,
    variable.name = "index",
    value.name = "value"
  )
  
  long[, index := as.integer(index)]
  
  target <- long[, .(value = mean(value, na.rm = TRUE)), by = index]
  
  setorder(target, index)
  
  return(target$value)
}

# ============================================================
# 5️⃣ Risk weight
# ============================================================
risk_weight <- function(risk, lo, hi, gamma = 2) {
  w <- (risk - lo) / (hi - lo)
  w <- pmin(pmax(w, 0), 1)
  w^gamma
}

# ============================================================
# 6️⃣ Individual recommendation
# ============================================================
generate_individual_recommendation <- function(row, target_vec, lo, hi) {
  
  x <- as.numeric(row[, ..feat_cols])
  w <- risk_weight(as.numeric(row$Predicted_Risk), lo, hi)
  
  new_x <- x + w * (target_vec - x)
  new_x[new_x < 0] <- 0
  
  tmp <- data.table(index = 0:191, value = new_x)
  tmp <- behavior_map[tmp, on = "index"]
  
  # normalize hourly
  tmp[, value := value / sum(value), by = .(day_type, hour)]
  
  # constraints
  night_hours <- c(22,23,0,1,2,3,4,5)
  day_hours   <- 8:20
  
  tmp[as.integer(hour) %in% night_hours & behavior == "Sleep",
      value := pmax(value, 0.75)]
  
  tmp[as.integer(hour) %in% night_hours & behavior == "MVPA",
      value := pmin(value, 0.05)]
  
  tmp[as.integer(hour) %in% day_hours & behavior == "SB",
      value := pmin(value, 0.75)]
  
  tmp[, value := value / sum(value), by = .(day_type, hour)]
  
  setorder(tmp, index)
  
  return(tmp$value)
}

# ============================================================
# 7️⃣ Generate all recommendations
# ============================================================
generate_all_recommendations_long <- function(df) {
  
  df <- copy(df)
  
  target_vec <- compute_low_profile(df)
  
  lo <- quantile(df$Predicted_Risk, 0.05, na.rm = TRUE)
  hi <- quantile(df$Predicted_Risk, 0.995, na.rm = TRUE)
  
  rec_list <- lapply(seq_len(nrow(df)), function(i) {
    
    new_x <- generate_individual_recommendation(
      df[i], target_vec, lo, hi
    )
    
    out <- data.table(
      `Participant ID` = df[i, `Participant ID`],
      Disease = df[i, Disease],
      index = 0:191,
      value = new_x
    )
    
    behavior_map[out, on = "index"]
  })
  
  return(rbindlist(rec_list))
}

# ============================================================
# 8️⃣ Plot function
# ============================================================
behavior_colors <- c(
  SB = "#33A02C",
  LPA = "#1F78B4",
  MVPA = "#E31A1C",
  Sleep = "#FFD92F"
)

plot_behavior_weekday_weekend_R <- function(obs, rec, disease, save_dir) {
  
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  
  behaviors <- c("MVPA", "LPA", "SB", "Sleep")
  day_types <- c("Weekday", "Weekend")
  
  plot_list <- list()
  
  for (b in behaviors) {
    for (d in day_types) {
      
      avg <- obs[behavior == b & day_type == d]
      opt <- rec[behavior == b & day_type == d]
      
      df <- merge(
        avg[, .(hour, avg = value * 100)],
        opt[, .(hour, opt = value * 100)],
        by = "hour"
      )
      
      df[, diff := opt - avg]
      
      p <- ggplot(df, aes(x = as.numeric(hour))) +
        geom_ribbon(aes(ymin = 0, ymax = pmax(diff, 0)),
                    fill = behavior_colors[b], alpha = 0.25) +
        geom_ribbon(aes(ymin = pmin(diff, 0), ymax = 0),
                    fill = "grey60", alpha = 0.18) +
        geom_line(aes(y = avg), color = behavior_colors[b], linewidth = 1.4) +
        geom_line(aes(y = opt), color = behavior_colors[b],
                  linewidth = 1.4, linetype = "dashed") +
        theme_classic(base_size = 16) +
        labs(title = paste(disease, "-", b, "-", d))
      
      plot_list[[paste(b, d)]] <- p
    }
  }
  
  final_plot <- wrap_plots(plot_list, ncol = 2)
  
  ggsave(
    filename = file.path(save_dir, paste0(disease, "_plot.tiff")),
    plot = final_plot,
    dpi = 300,
    width = 12,
    height = 16
  )
}

# ============================================================
# 🔟 Run all diseases
# ============================================================
out_dir <- ".../Individual_recommendations"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (disease in disease_list) {
  
  message("Processing: ", disease)
  
  df_d <- dt[Disease == disease]
  if (nrow(df_d) == 0) next
  
  obs <- compute_population_profile(df_d)
  
  rec_long <- generate_all_recommendations_long(df_d)
  
  rec_wide <- dcast(
    rec_long,
    `Participant ID` ~ index,
    value.var = "value"
  )
  
  setnames(rec_wide, as.character(0:191), feat_cols)
  
  rec <- compute_population_profile(rec_wide)
  
  plot_behavior_weekday_weekend_R(
    obs, rec, disease, file.path(out_dir, disease)
  )
}