# ---------------------------------------------------------------------------------------------- #
# Bayesian SEM with brms: Beta mediators + hurdle lognormal abundance
# Loop for Functional Cluster
# Author: [YM]
# ---------------------------------------------------------------------------------------------- #

library(tidyverse)
library(brms)
library(posterior)
library(bayesplot)

# -----------------------------
# Load and prepare data
# -----------------------------
# SINDEX_ZEROS contains counts and zeros, for the species x year not seen in a site when seen in another.
ubmsdata <- read.table("DATA/SINDEX_ZERO_MS.txt", header = TRUE, sep = "\t") %>%
  filter(YEAR > 2020) %>%
  mutate(
    SPECIES = case_when(
      SPECIES == "Glaucopsyche sp." ~ "Glaucopsyche melanops",
      SPECIES == "Melanargia sp."   ~ "Melanargia lachesis",
      SPECIES == "Satyridae"   ~ "Pyronia sp",
      
      TRUE                          ~ SPECIES
    )
  )

space    <- read.delim("DATA/space_data.txt", sep = "\t")
clusters <- read.csv2("DATA/Species_clusters.csv")

space <- space %>%
  mutate(
    C3.proportion   = C3ha / TotalA, # ornamental veg., hugh hyric resources, highly managed (sega y riego)
    P.proportion    = PUha / TotalA, # native veg. with spontaneous herbs and flowers, low management(one sega por año??)
    Herb.proportion = pmin(HERbha / TotalA, 1), # native herbs and flowers managed (water) but no siega (only once?? - IN Biod types)
    Viv.proportion  = VIVha / TotalA # ??
  )

ubmsdata <- ubmsdata %>%
  left_join(clusters, by = "SPECIES") %>%
  filter(!is.na(Cluster)) %>%
  left_join(space, by = "SITE_ID") %>%
  mutate(
    presence  = as.integer(SINDEX > 0),
    abundance = SINDEX,   # continuous proxy, >=0
    Cluster   = as.factor(Cluster)
  ) %>%
  rename(
    Area       = TotalA,
    Conn       = C500,
    C3_grass   = C3.proportion,
    P_grass    = P.proportion,
    Herb_grass = Herb.proportion,
    V_grass    = Viv.proportion
  )

na.ubmsdata<-ubmsdata %>% filter(is.na(C3_grass)) %>% distinct(SITE_ID)

# -----------------------------
# Standardize predictors + prepare mediators
# -----------------------------
ubmsdata <- ubmsdata %>%
  mutate(
    Type = factor(Type),
    # Smithson–Verkuilen transformation to keep proportions in ]0,1[
    C3_beta   = (C3_grass   * (n() - 1) + 0.5) / n(),
    P_beta    = (P_grass    * (n() - 1) + 0.5) / n(),
    Herb_beta = (Herb_grass * (n() - 1) + 0.5) / n(),
    Viv_beta  = (V_grass    * (n() - 1) + 0.5) / n(),
    # Standardize exogenous predictors
    z_Area = as.numeric(scale(Area)),
    z_Conn = as.numeric(scale(Conn)),
    # Logit-transform mediators for use as predictors in abundance model, then standardize
    zC3_logit   = as.numeric(scale(qlogis(C3_beta))),
    zP_logit    = as.numeric(scale(qlogis(P_beta))),
    zHerb_logit = as.numeric(scale(qlogis(Herb_beta))),
    zViv_logit  = as.numeric(scale(qlogis(Viv_beta)))
  )

# -----------------------------
# Function to fit SEM per cluster
# -----------------------------
fit_cluster_sem <- function(cluster_id, data) {
  
  message(">> Fitting SEM for cluster: ", cluster_id)
  
  df <- data %>% filter(Cluster == cluster_id)
  
  # --- Model formulas ---
  bf_C3   <- bf(C3_beta   ~ z_Area + Type, family = Beta(link = "logit"))
  bf_P    <- bf(P_beta    ~ z_Area + Type, family = Beta(link = "logit"))
  bf_Herb <- bf(Herb_beta ~ z_Area + Type, family = Beta(link = "logit"))
  bf_Viv  <- bf(Viv_beta  ~ z_Area + Type, family = Beta(link = "logit"))
  
  bf_abund <- bf(
    abundance ~ z_Area + z_Conn + Type + zC3_logit + zP_logit + zHerb_logit + zViv_logit + (1 | YEAR),
    hu       ~ z_Area + z_Conn + Type + zC3_logit + zP_logit + zHerb_logit + zViv_logit,
    family = hurdle_lognormal()
  )
  
  # --- Fit SEM ---
  fit <- brm(
    bf_C3 + bf_P + bf_Herb + bf_Viv + bf_abund,
    data = df,
    chains = 4, iter = 100000, warmup = 10000, cores = parallel::detectCores(),
    control = list(adapt_delta = 0.99, max_treedepth = 12),
    seed = 20251004
  )
  
  # --- Save ---
  saveRDS(fit, file = paste0("Results/Abundance_brms_SEM_Cluster", cluster_id, ".rds"))
  message(">> Model for cluster ", cluster_id, " saved successfully.")
  return(fit)
}

# -----------------------------
# Run for all clusters
# -----------------------------
cluster_levels <- levels(ubmsdata$Cluster)

fits <- purrr::map(cluster_levels, ~ fit_cluster_sem(.x, ubmsdata))
names(fits) <- cluster_levels

# ---------------------------------------------------------------------------------------------- #
# Post-processing Bayesian SEM (hurdle lognormal + beta mediators) per cluster
# Author: [YM]
# ---------------------------------------------------------------------------------------------- #

# ==== 1) Load models saved ====
fit_c4 <- readRDS("Results/Abundance_brms_SEM_Cluster4.rds")
fit_c2 <- readRDS("Results/Abundance_brms_SEM_Cluster2.rds")
fit_c3 <- readRDS("Results/Abundance_brms_SEM_Cluster3.rds")

# ==== 2) Basic diagnostics ====
# Always check Rhat < 1.01 and reasonable neff_ratio
summary(fit_c2)                          # detailed overview // check Rhat in summary for model fit Rhat ≈ 1.00 → chains have converged well.
summary(fit_c3)
# check R hat for each parameter SHOULD be < 1.05
range(neff_ratio(fit_c2), na.rm = TRUE)  # Effective sample size ratio
levels(fit_c2$data$Type)

# ==== 3) Posterior predictive checks ====
# Includes both zeros and positive abundances (hurdle model)
pp_check(fit_c2, resp = "abundance")

# ==== 4) Extract fixed effects ====
# Positive abundance part (conditional on presence)
abund_fx_c2 <- as.data.frame(fixef(fit_c2, resp = "abundance"))
abund_fx_c2[order(abs(abund_fx_c2$Estimate), decreasing = TRUE), ]

# Hurdle (zeros) part: probability of abundance = 0
hu_fx_c2 <- as.data.frame(fixef(fit_c2, resp = "abundance", dpar = "hu"))
hu_fx_c2[order(abs(hu_fx_c2$Estimate), decreasing = TRUE), ]

# Mediators (beta regressions)
c3_fx_c2   <- as.data.frame(fixef(fit_c2, resp = "C3beta"))
p_fx_c2    <- as.data.frame(fixef(fit_c2, resp = "Pbeta"))
herb_fx_c2 <- as.data.frame(fixef(fit_c2, resp = "Herbbeta"))
viv_fx_c2  <- as.data.frame(fixef(fit_c2, resp = "Vivbeta"))

# ==== 5) Predicted probability of zeros ====
# Estimated hurdle probability per observation
p_zero_c2 <- fitted(fit_c2, resp = "abundance", dpar = "hu")
head(p_zero_c2)

# You can merge this with the original data (subset by cluster)
# ubmsdata_c2 <- subset(ubmsdata, Cluster == levels(ubmsdata$Cluster)[2])
# ubmsdata_c2$p_zero_est <- p_zero_c2[, "Estimate"]

# ==== 6) Marginal effects for plotting ====
# Abundance conditional on presence
ce_abund_c2 <- conditional_effects(fit_c2, resp = "abundance")
# Probability of zero (absence process)
ce_hu_c2    <- conditional_effects(fit_c2, resp = "abundance", dpar = "hu")
# plot(ce_abund_c2); plot(ce_hu_c2)

# ==== 7) Helper function to extract all effects in one go ====
extract_all_effects <- function(fit, tag) {
  out <- list(
    abundance = as.data.frame(fixef(fit, resp = "abundance")),
    hu        = as.data.frame(fixef(fit, resp = "abundance", dpar = "hu")),
    C3beta    = as.data.frame(fixef(fit, resp = "C3beta")),
    Pbeta     = as.data.frame(fixef(fit, resp = "Pbeta")),
    Herbbeta  = as.data.frame(fixef(fit, resp = "Herbbeta")),
    Vivbeta   = as.data.frame(fixef(fit, resp = "Vivbeta"))
  )
  # optional: save as CSV for later inspection
  dir.create("Results/Tables", showWarnings = FALSE)
  purrr::iwalk(out, ~ write.csv(.x, file = file.path("Results/Tables", paste0(tag, "_", .y, ".csv")), row.names = TRUE))
  out
}

# Example run:
fx_c2 <- extract_all_effects(fit_c2, "Cluster2")
