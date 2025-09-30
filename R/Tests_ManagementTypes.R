# ---------------------------------------------------------------------------------------------- #
# Script to test Urban Management Strategies: Presence and Abundances versus management types 
# Author: [YM]
# Inputs: 
#   - Data from uBMS license agreement 
# Outputs:
#   - XX
#  ---------------------------------------------------------------------------------------------- #

library(tidyverse)
library(ggplot2)
library(glmmTMB)
library(corrplot) 
library(DHARMa)
library(car)
library(piecewiseSEM)

library(psych)
library(ggeffects)

# -----------------------------
# Load and prepare data
# -----------------------------
ubmsdata <- read.table("DATA/SINDEX_ZEROS.txt", header = TRUE, sep = "\t") %>%
  filter(YEAR > 2020)%>%
  mutate(
    SPECIES = case_when(
      SPECIES == "Glaucopsyche sp." ~ "Glaucopsyche melanops",
      SPECIES == "Melanargia sp."   ~ "Melanargia lachesis",
      TRUE                          ~ SPECIES
    )
  )
# Glaucopsyche sp. => Glaucopsyche melanops
# Melanargia sp. => Melanargia lachesis
space <- read.delim("DATA/space_data.txt", sep="\t")
clusters <- read.csv2("DATA/Species_clusters.csv")

# "Since we are interested in the effect of each strategy, without accounting for the Area occupied by each one 
# -----------------------------
# Calculate proportion of area of each strategy type. 
# -----------------------------
space <- space %>%
  mutate(
    C3.proportion = C3ha / TotalA,
    P.proportion = PUha / TotalA,
    Herb.proportion = pmin(HERbha / TotalA, 1), #bc there is an error
    Viv.proportion = VIVha / TotalA,
  )

ubmsdata <- ubmsdata %>%
  left_join(clusters, by = "SPECIES") %>%
  filter(!is.na(Cluster)) %>%
  left_join(space, by = "SITE_ID") %>%
  mutate(
    presence = as.numeric(SINDEX > 0),
    Cluster = as.factor(Cluster)
  )

# -----------------------------
# Correlation matrix of landscape variables
# -----------------------------
# If total area is included in the SEM, then mathematically the area of vegetation should be in proportions
# then H would be type of vegetation (not directly the area of which, but also incuded (nested))
p.var <- space[, c(4:9,24:27)]
M <- corr.test(p.var)
corrplot(M$r, method = "number", type = "upper", sig.level = c(0.001, 0.01, 0.05), 
         insig = "label_sig", pch.cex = 0.9, pch.col = "grey20", order = "original")
corrplot(M$r, p.mat = M$p, type = "upper", sig.level = c(0.001, 0.01, 0.05), 
         insig = "label_sig", pch.cex = 0.9, pch.col = "grey", order = "original")


ubmsdata <- ubmsdata %>%
  rename(
    Area       = TotalA,
    Conn       = C500,
    C3_grass   = C3.proportion,
    P_grass    = P.proportion,
    Herb_grass = Herb.proportion,
    V_grass    = Viv.proportion
  )

hist(ubmsdata$C3_grass)
hist(ubmsdata$P_grass)
hist(ubmsdata$Herb_grass)
hist(ubmsdata$V_grass)
# -----------------------------
# SEM_B: garden type with covariance with TotalA

# 1. Presence
# Models for the mediators (vegetation proportions)
 # Vegetation proportions were modeled with fixed-effects GLMs, as random intercepts did not improve fit and caused convergence problems.
 # Also no site or year as random, since the proportion of vegetation is part of each site
 # models m_V1 - m_V4: residuals do not follow a Gaussian distribution, so we applied beta family (logit transformation was not Gaussian either)
 # We modeled vegetation proportions using a Beta family with a logit link, because these variables are bounded between 0 and 1 and showed deviations from normality under Gaussian models. The Beta distribution is specifically designed for continuous proportions, allowing for flexible shapes (skewed, symmetric, U-shaped) while keeping estimates within the valid range.
 
ubmsdata <- ubmsdata %>% # beta family is  ]0,1[]
  mutate(
    C3_beta   = (C3_grass   * (n() - 1) + 0.5) / n(),
    P_beta    = (P_grass    * (n() - 1) + 0.5) / n(),
    Herb_beta = (Herb_grass * (n() - 1) + 0.5) / n(),
    Viv_beta  = (V_grass  * (n() - 1) + 0.5) / n()
  )

m_V1 <- glmmTMB(C3_beta   ~ Area + Type, data = ubmsdata, family = beta_family())
m_V2 <- glmmTMB(P_beta    ~ Area + Type, data = ubmsdata, family = beta_family())
m_V3 <- glmmTMB(Herb_beta ~ Area + Type, data = ubmsdata, family = beta_family())
m_V4 <- glmmTMB(Viv_beta    ~ Area + Type, data = ubmsdata, family = beta_family())

models <- list(m_V1 = m_V1, m_V2 = m_V2, m_V3 = m_V3, m_V4 = m_V4)

# Summary 
for (i in names(models)) {
  cat("\n========================\n")
  cat("Summary for", i, "\n")
  print(summary(models[[i]]))
  cat("========================\n\n")
}
# Diagnostics
pdf("BetaModels_DHARMa_diagnostics.pdf", width = 7, height = 7)
for (i in names(models)) {
  res <- simulateResiduals(models[[i]], n = 1000)
  plot(res, main = paste("Residual diagnostics", i))
}
dev.off()

# Model for the response (presence/absence) 
m_resp <- glmmTMB(presence ~ Area + Conn + Type + C3_beta + P_beta + Viv_beta + Herb_beta + (1|YEAR), # interactions not significant and lower AIC
                  data = ubmsdata,
                  family = binomial(link = "logit"))

AIC(m_resp); summary(m_resp)

# Build the SEM_B with covariance between Area and Type
sem_B <- psem(
  m_V1, m_V2, m_V3, m_V4, m_resp,
  Area %~~% Type
)


# SEM outputs
summary(sem_B, conserve = TRUE)        # coefficients and global Fisher's C test
anova(sem_B)          # d-sep tests of implied independencies
rsquared(sem_B)       # marginal and conditional R² for each submodel
coefs(sem_B)          # standardized path coefficients



# -----------------------------