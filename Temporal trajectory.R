# =========================================================
# 1. Libraries
# =========================================================
library(data.table)
library(dplyr)
library(lubridate)
library(survival)
library(mstate)
library(igraph)
library(ggplot2)
library(scales)
library(ggVennDiagram)

# =========================================================
# 2. Load Data
# =========================================================
dat <- fread("analysis_master_dataset.csv")

date_cols <- c(
  "Start time of wear","obesity_date","t2d_date",
  "dlm_date","htn_date","Date of death","region_cutoff"
)

dat[, (date_cols) := lapply(.SD, as.Date), .SDcols = date_cols]

dat[, `:=`(
  Sex = factor(Sex),
  ethnicity_top = factor(ethnicity_top),
  smoking_cat = factor(smoking_cat)
)]

# =========================================================
# 3. Time-to-event construction
# =========================================================
dat[, `:=`(
  t_obesity = as.numeric(obesity_date - `Start time of wear`),
  t_t2d     = as.numeric(t2d_date - `Start time of wear`),
  t_htn     = as.numeric(htn_date - `Start time of wear`),
  t_dlm     = as.numeric(dlm_date - `Start time of wear`)
)]

dat[, t_end := pmin(
  as.numeric(`Date of death` - `Start time of wear`),
  as.numeric(region_cutoff - `Start time of wear`),
  na.rm = TRUE
)]

dt <- dat[, .(
  id = .I,
  t_obesity, t_t2d, t_htn, t_dlm, t_end,
  age, Sex, ethnicity_top, smoking_cat
)]

# =========================================================
# 4. Step 1: Cox + Bonferroni
# =========================================================
disease_list <- c("obesity","t2d","htn","dlm")
cox_results <- list()

for(d in disease_list){
  
  t_var <- paste0("t_", d)
  
  tmp <- dt[is.na(get(t_var)) | get(t_var) > 0]
  
  tmp[, event := as.numeric(!is.na(get(t_var)) & get(t_var) <= t_end)]
  tmp[, time  := ifelse(event==1, get(t_var), t_end)]
  
  fit <- coxph(Surv(time, event) ~ age + Sex + ethnicity_top + smoking_cat,
               data = tmp)
  
  cox_results[[d]] <- data.frame(
    disease = d,
    HR = exp(coef(fit)[1]),
    p = summary(fit)$coefficients[1,5]
  )
}

cox_df <- rbindlist(cox_results)
cox_df$p_adj <- p.adjust(cox_df$p, method="bonferroni")

valid_disease <- cox_df[HR > 1 & p_adj < 0.05]$disease

# =========================================================
# 5. Step 2–3: Temporal order + Logistic
# =========================================================
pairs <- expand.grid(D1=valid_disease, D2=valid_disease, stringsAsFactors=FALSE)
pairs <- pairs[pairs$D1 != pairs$D2,]

results <- list()
threshold <- 0.05 / nrow(pairs)

for(i in 1:nrow(pairs)){
  
  d1 <- pairs$D1[i]
  d2 <- pairs$D2[i]
  
  t1 <- paste0("t_", d1)
  t2 <- paste0("t_", d2)
  
  sub <- dt[
    !is.na(get(t1)) & !is.na(get(t2)) &
      get(t1) <= t_end & get(t2) <= t_end
  ]
  
  if(nrow(sub) < 20) next
  
  n_total <- nrow(sub)
  n_order <- nrow(sub[get(t1) < get(t2)])
  
  binom_p <- binom.test(n_order, n_total, p=0.5)$p.value
  if(binom_p >= threshold) next
  
  sub[, D1 := as.numeric(!is.na(get(t1)) & get(t1) <= t_end)]
  sub[, D2 := as.numeric(!is.na(get(t2)) & get(t2) <= t_end)]
  
  fit <- glm(D2 ~ D1 + age + Sex + ethnicity_top + smoking_cat,
             data = sub, family = binomial)
  
  coef_tab <- summary(fit)$coefficients
  if(!("D1" %in% rownames(coef_tab))) next
  
  beta <- coef_tab["D1",1]
  se   <- coef_tab["D1",2]
  p    <- coef_tab["D1",4]
  
  OR <- exp(beta)
  CI_low  <- exp(beta - 1.96 * se)
  CI_high <- exp(beta + 1.96 * se)
  
  if(OR > 1 & p < threshold){
    results[[length(results)+1]] <- data.frame(
      from=d1, to=d2,
      OR=OR, CI_low=CI_low, CI_high=CI_high,
      p=p, n=n_order
    )
  }
}

edge_df <- rbindlist(results)

# Keep strongest direction per pair
edge_df[, pair := paste(pmin(from,to), pmax(from,to), sep="_")]
edge_single <- edge_df[, .SD[which.max(OR)], by=pair]
edge_single[, pair := NULL]

# Save
write.csv(edge_single, "edge_df.csv", row.names = FALSE)

# =========================================================
# 6. Graph Visualization
# =========================================================
node_df <- data.frame(
  name = c("obesity","t2d","htn","dlm"),
  label = c("Obesity","Type 2 Diabetes","Hypertension","Dyslipidemia"),
  HR = cox_df$HR[match(c("obesity","t2d","htn","dlm"), cox_df$disease)]
)

g <- graph_from_data_frame(edge_single, vertices=node_df, directed=TRUE)

set.seed(123)
lay <- layout_with_fr(g)

V(g)$x <- lay[,1]
V(g)$y <- lay[,2]

node_radius <- 0.06

node_plot <- data.frame(
  name = V(g)$name,
  label = node_df$label,
  HR = V(g)$HR,
  x = V(g)$x,
  y = V(g)$y
)

edge_list <- as_data_frame(g, what="edges")

edge_list$x_from <- node_plot$x[match(edge_list$from, node_plot$name)]
edge_list$y_from <- node_plot$y[match(edge_list$from, node_plot$name)]
edge_list$x_to   <- node_plot$x[match(edge_list$to, node_plot$name)]
edge_list$y_to   <- node_plot$y[match(edge_list$to, node_plot$name)]

edge_list$dx <- edge_list$x_to - edge_list$x_from
edge_list$dy <- edge_list$y_to - edge_list$y_from
edge_list$dist <- sqrt(edge_list$dx^2 + edge_list$dy^2)

edge_list$ux <- edge_list$dx / edge_list$dist
edge_list$uy <- edge_list$dy / edge_list$dist

edge_list$x_start <- edge_list$x_from + edge_list$ux * node_radius
edge_list$y_start <- edge_list$y_from + edge_list$uy * node_radius

arrow_head_length <- 0.05
edge_list$x_end <- edge_list$x_to - edge_list$ux * (node_radius + arrow_head_length)
edge_list$y_end <- edge_list$y_to - edge_list$uy * (node_radius + arrow_head_length)

edge_list$x_mid <- edge_list$x_start + (edge_list$x_end - edge_list$x_start)/3
edge_list$y_mid <- edge_list$y_start + (edge_list$y_end - edge_list$y_start)/3

p <- ggplot() +
  geom_segment(data=edge_list,
               aes(x=x_start,y=y_start,xend=x_end,yend=y_end,color=OR),
               size=2.5,
               arrow=arrow(length=unit(5,"mm"), type="closed")) +
  geom_label(data=edge_list,
             aes(x=x_mid,y=y_mid,label=n),
             size=5) +
  geom_point(data=node_plot,
             aes(x=x,y=y,fill=HR),
             size=20, shape=21, color="black") +
  geom_text(data=node_plot,
            aes(x=x,y=y,label=label),
            size=5) +
  scale_color_gradient(low="#bdd7e7", high="#08306b") +
  scale_fill_gradient(low="#fee5d9", high="#a50f15") +
  theme_void()

ggsave("disease_trajectory.pdf", p, width=10, height=7)
