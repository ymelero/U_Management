# ---------------------------------------------------------------------------------------------- #
# Back-transform estimates (Estimate, SE, CI) from logit to proportion scale
# for the Beta submodels (C3, P, Herb, Viv)
# AND from log / logit scale for abundance and hurdle components
# Predict (marginal preditions) from SEM models 
# Author: [YM]
# ---------------------------------------------------------------------------------------------- #

library(tidyverse)
library(purrr)
library(brms)
detach("package:ggeffects")
detach("package:gridExtra")

fit_c2 <- readRDS("Results/Abundance_brms_SEM_Cluster4.rds")
fx_all <- as.data.frame(fixef(fit_c2)) %>% mutate(across(everything(), as.numeric))

# 1. Site variables
inv_logit <- function(x) plogis(as.numeric(x))
inv_exp   <- function(x) exp(as.numeric(x))

# --- Back-transform one row for Beta models (logit → proportion) ---
backtransform_row <- function(estimate, se, l95, u95) {
  tibble(
    Estimate_bt = inv_logit(estimate),
    SE_approx   = (inv_logit(estimate + se) - inv_logit(estimate - se)) / 2,
    CI95_low    = inv_logit(l95),
    CI95_high   = inv_logit(u95)
  )
}

# --- Apply to one mediator (C3, P, Herb, Viv) ---
extract_backtrans <- function(fx, mediator) {
  rows <- grep(paste0("^", mediator, "_"), rownames(fx), value = TRUE)
  
  map_dfr(rows, function(rn) {
    vals <- as.numeric(fx[rn, c("Estimate", "Est.Error", "Q2.5", "Q97.5")])
    trans <- backtransform_row(vals[1], vals[2], vals[3], vals[4])
    tibble(
      Mediator = mediator,
      Parameter = rn,
      Estimate = vals[1],
      SE = vals[2],
      L95 = vals[3],
      U95 = vals[4],
      Estimate_bt = trans$Estimate_bt,
      SE_approx = trans$SE_approx,
      CI95_low = trans$CI95_low,
      CI95_high = trans$CI95_high
    )
  })
}

# --- Run for all Beta mediators ---
mediators_beta <- c("C3beta", "Pbeta", "Herbbeta", "Vivbeta")
bt_beta <- map_dfr(mediators_beta, ~ extract_backtrans(fx_all, .x))


# 2. Back-transform abundance (log) and hurdle (logit) parts
backtransform_abundance_row <- function(estimate, se, l95, u95, type) {
  if (type == "log") {
    tibble(
      Estimate_bt = inv_exp(estimate),
      SE_approx   = (inv_exp(estimate + se) - inv_exp(estimate - se)) / 2,
      CI95_low    = inv_exp(l95),
      CI95_high   = inv_exp(u95)
    )
  } else if (type == "logit") {
    tibble(
      Estimate_bt = inv_logit(estimate),
      SE_approx   = (inv_logit(estimate + se) - inv_logit(estimate - se)) / 2,
      CI95_low    = inv_logit(l95),
      CI95_high   = inv_logit(u95)
    )
  }
}

extract_backtrans_abundance <- function(fx, mediator, type) {
  rows <- grep(paste0("^", mediator, "_"), rownames(fx), value = TRUE)
  
  map_dfr(rows, function(rn) {
    vals <- as.numeric(fx[rn, c("Estimate", "Est.Error", "Q2.5", "Q97.5")])
    trans <- backtransform_abundance_row(vals[1], vals[2], vals[3], vals[4], type)
    tibble(
      Mediator = mediator,
      Parameter = rn,
      Estimate = vals[1],
      SE = vals[2],
      L95 = vals[3],
      U95 = vals[4],
      Estimate_bt = trans$Estimate_bt,
      SE_approx = trans$SE_approx,
      CI95_low = trans$CI95_low,
      CI95_high = trans$CI95_high
    )
  })
}

# --- Run for both abundance parts ---
bt_abund <- extract_backtrans_abundance(fx_all, "abundance", "log")
bt_hu    <- extract_backtrans_abundance(fx_all, "hu_abundance", "logit")

# --- Merge all results into one unified table ---
bt_all_full <- bind_rows(bt_beta, bt_abund, bt_hu) %>%
  arrange(Mediator, Parameter)
write.table(bt_all_full, "Results/backtransdormed_C4.txt", sep=";", row.names = F)

# 3. Predict & Plot
library(ggeffects)
library(gridExtra)

# --- Function for main responses (Beta and abundance) ---
get_preds <- function(mediator) {
  ce <- conditional_effects(fit_c2, effects = "Type", resp = mediator)
  
  as.data.frame(ce[[1]]) %>%
    rename(predicted = estimate__, conf.low = lower__, conf.high = upper__, x = Type) %>%
    mutate(Mediator = mediator)
}

# --- Function for hurdle component (distributional parameter 'hu') ---
get_preds_hurdle <- function() {
  ce <- conditional_effects(fit_c2, effects = "Type", resp = "abundance", dpar = "hu")
  
  as.data.frame(ce[[1]]) %>%
    rename(predicted = estimate__, conf.low = lower__, conf.high = upper__, x = Type) %>%
    mutate(Mediator = "hu_abundance")
}

# --- Mediators present in the model ---
mediators_main <- c("C3beta", "Pbeta", "Herbbeta", "Vivbeta", "abundance")

# --- Run predictions for main responses ---
preds_main <- purrr::map_dfr(mediators_main, get_preds)

# --- Run predictions for hurdle (probability of zero) ---
preds_hu <- get_preds_hurdle()

# --- Combine all predictions ---
preds_all <- bind_rows(preds_main, preds_hu)

# --- Format for plotting ---
preds_all <- preds_all %>%
  mutate(
    Type = factor(x, levels = c("BIO", "PATRIMONIAL", "USOS"),
                  labels = c("Bio.", "Hist.", "Soc.")),
    Mediator = factor(Mediator,
                      levels = c("C3beta", "Pbeta", "Herbbeta", "Vivbeta",
                                 "abundance", "hu_abundance"),
                      labels = c("C3", "P", "Herb", "Viv", "Abundance", "Hurdle (P0)"))
  )

# --- Quick check ---
head(preds_all)

# --- Plot (optional) ---
ggplot(preds_all, aes(x = Type, y = predicted, color = Mediator)) +
  geom_point(position = position_dodge(width = 0.8), size = 2.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                position = position_dodge(width = 0.8), width = 0.25) +
  labs(x = "Space Type", y = "Proportion") +
  scale_color_brewer(palette = "Set2") +
  theme_classic(base_size = 12) +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(size = 11),
    axis.title = element_text(size = 12)
  )


# separado
preds_props <- preds_all %>%
  filter(Mediator %in% c("C3", "P", "Herb", "Viv"))

p1 <- ggplot(preds_props, aes(x = Type, y = predicted, color = Mediator)) +
  geom_point(position = position_dodge(width = 0.8), size = 2.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                position = position_dodge(width = 0.8), width = 0.25) +
  labs(x = "Space Type", y = "Proportions") +
  scale_color_brewer(palette = "Set2") +
  theme_classic(base_size = 12) +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(size = 11),
    axis.title = element_text(size = 12)
  )
p1

preds_hu_plot <- preds_all %>%
  filter(Mediator == "Hurdle (P0)")

p2 <- ggplot(preds_hu_plot, aes(x = Type, y = 1-predicted, color = Type)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = 1- conf.low, ymax = 1-conf.high), width = 0.25) +
  labs(x = "Space Type", y = "Presence") +
  scale_color_brewer(palette = "Set2") +
  theme_classic(base_size = 12) +
  theme(legend.title = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(size = 11),
    axis.title = element_text(size = 12)
  )
p2

preds_abund <- preds_all %>%
  filter(Mediator == "Abundance")

p3 <- ggplot(preds_abund, aes(x = Type, y = predicted, color = Type)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.25) +
  labs(x = "Space Type", y = "Abundance") +
  scale_color_brewer(palette = "Set2") +
  theme_classic(base_size = 12) +
  theme(legend.title = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(size = 11),
    axis.title = element_text(size = 12)
  ) 
p3 
grid.arrange(p1, p2,p3, ncol = 3, nrow = 1)

# Fig caption: Figure X. Mean and ±95% credible intervals from the Bayesian structural equation model for Cluster X. 
# of green urban space type (Biological, Historical, and Social) on: the proportion of vegetation composition (Proportions), presence probability, and butterfly abundance. 