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
ubmsdata <- read.table("DATA/SINDEX_ZEROS.txt", header = TRUE, sep = "\t") %>%
  filter(YEAR > 2020) %>%
  mutate(
    SPECIES = case_when(
      SPECIES == "Glaucopsyche sp." ~ "Glaucopsyche melanops",
      SPECIES == "Melanargia sp."   ~ "Melanargia lachesis",
      TRUE                          ~ SPECIES
    )
  )

space    <- read.delim("DATA/space_data.txt", sep = "\t")
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

# -----------------------------
# Standardize predictors + prepare mediators
# -----------------------------
ubmsdata <- ubmsdata %>%
  mutate(
    Type = factor(Type),
    # Smithson–Verkuilen transformation to keep proportions in (0,1)
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
  
  return(fit)
}

# -----------------------------
# Run for all clusters
# -----------------------------
cluster_levels <- levels(ubmsdata$Cluster)

fits <- purrr::map(cluster_levels, ~ fit_cluster_sem(.x, ubmsdata))
names(fits) <- cluster_levels
