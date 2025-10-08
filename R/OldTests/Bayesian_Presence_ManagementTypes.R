# ---------------------------------------------------------------------------------------------- #
# Bayesian SEM with brms: standardized predictors
# Author: [YM]
# ---------------------------------------------------------------------------------------------- #

library(tidyverse)
library(brms)
library(posterior)
library(bayesplot)

# -----------------------------
# Load and prepare data
# -----------------------------
ubmsdata <- read.table("DATA/SINDEX_ZEROS.txt", header = TRUE, sep = "\t") %>%
  filter(YEAR > 2020) %>%
  mutate(
    SPECIES = case_when(
      SPECIES == "Glaucopsyche sp." ~ "Glaucopsyche melanops",
      SPECIES == "Melanargia sp."   ~ "Melanargia lachesis",
      TRUE                          ~ SPECIES
    )
  )

space    <- read.delim("DATA/space_data.txt", sep="\t")
clusters <- read.csv2("DATA/Species_clusters.csv")

space <- space %>%
  mutate(
    C3.proportion   = C3ha / TotalA,
    P.proportion    = PUha / TotalA,
    Herb.proportion = pmin(HERbha / TotalA, 1),
    Viv.proportion  = VIVha / TotalA
  )

ubmsdata <- ubmsdata %>%
  left_join(clusters, by = "SPECIES") %>%
  filter(!is.na(Cluster)) %>%
  left_join(space, by = "SITE_ID") %>%
  mutate(
    presence = as.integer(SINDEX > 0),
    Cluster  = as.factor(Cluster)
  ) %>%
  rename(
    Area       = TotalA,
    Conn       = C500,
    C3_grass   = C3.proportion,
    P_grass    = P.proportion,
    Herb_grass = Herb.proportion,
    V_grass    = Viv.proportion
  )

# -----------------------------
# Standardize continuous vars
# -----------------------------
eps <- 1e-6
ubmsdata <- ubmsdata %>%
  mutate(
    Type = factor(Type),
    # Beta-corrected proportions
    C3_beta   = (C3_grass   * (n() - 1) + 0.5) / n(),
    P_beta    = (P_grass    * (n() - 1) + 0.5) / n(),
    Herb_beta = (Herb_grass * (n() - 1) + 0.5) / n(),
    Viv_beta  = (V_grass    * (n() - 1) + 0.5) / n()
  ) %>%
  # Z-scores for comparability
  mutate(across(c(Area, Conn, C3_beta, P_beta, Herb_beta, Viv_beta),
                ~ as.numeric(scale(.)),
                .names = "z_{.col}"))

# -----------------------------
# SEM specification in brms
# -----------------------------
bf_C3   <- bf(z_C3_beta   ~ z_Area + Type, family = gaussian())
bf_P    <- bf(z_P_beta    ~ z_Area + Type, family = gaussian())
bf_Herb <- bf(z_Herb_beta ~ z_Area + Type, family = gaussian())
bf_Viv  <- bf(z_Viv_beta  ~ z_Area + Type, family = gaussian())

bf_pres <- bf(
  presence ~ z_Area + z_Conn + Type + z_C3_beta + z_P_beta + z_Herb_beta + z_Viv_beta + (1 | YEAR),
  family = bernoulli(link = "logit")
)

# -----------------------------
# Fit the model
# -----------------------------
fit_sem <- brm(
  bf_C3 + bf_P + bf_Herb + bf_Viv + bf_pres,
  data = ubmsdata,
  chains = 4, iter = 3000, warmup = 1000, cores = parallel::detectCores(),
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  seed = 20251001
)

# -----------------------------
# Diagnostics
# -----------------------------
summary(fit_sem)

pp_check(fit_sem, resp = "presence")
pp_check(fit_sem, resp = "z_C3_beta")
pp_check(fit_sem, resp = "z_P_beta")
pp_check(fit_sem, resp = "z_Herb_beta")
pp_check(fit_sem, resp = "z_Viv_beta")

# -----------------------------
# Effects: direct, indirect, total (for Area -> presence)
# -----------------------------
draws <- as_draws_df(fit_sem)

coef_area_C3   <- "b_z_C3_beta_z_Area"
coef_pres_C3   <- "b_presence_z_C3_beta"

coef_area_P    <- "b_z_P_beta_z_Area"
coef_pres_P    <- "b_presence_z_P_beta"

coef_area_Herb <- "b_z_Herb_beta_z_Area"
coef_pres_Herb <- "b_presence_z_Herb_beta"

coef_area_Viv  <- "b_z_Viv_beta_z_Area"
coef_pres_Viv  <- "b_presence_z_Viv_beta"

coef_pres_area <- "b_presence_z_Area"

# Indirect effects
ind_C3   <- draws[[coef_area_C3]]   * draws[[coef_pres_C3]]
ind_P    <- draws[[coef_area_P]]    * draws[[coef_pres_P]]
ind_Herb <- draws[[coef_area_Herb]] * draws[[coef_pres_Herb]]
ind_Viv  <- draws[[coef_area_Viv]]  * draws[[coef_pres_Viv]]

# Direct and total effects of Area
dir_Area <- draws[[coef_pres_area]]
tot_Area <- dir_Area + ind_C3 + ind_P + ind_Herb + ind_Viv

# Summary helper
post_summary <- function(x) {
  c(median = median(x), 
    l95 = quantile(x, 0.025), 
    u95 = quantile(x, 0.975))
}

eff_table <- rbind(
  Area_ind_via_C3   = post_summary(ind_C3),
  Area_ind_via_P    = post_summary(ind_P),
  Area_ind_via_Herb = post_summary(ind_Herb),
  Area_ind_via_Viv  = post_summary(ind_Viv),
  Area_direct       = post_summary(dir_Area),
  Area_total        = post_summary(tot_Area)
) %>% as.data.frame()

print(round(eff_table, 3))

# -----------------------------
# Save model
# -----------------------------
saveRDS(fit_sem, file = "Results/Presence_brms_sem_standardized.rds")
