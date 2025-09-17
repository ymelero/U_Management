# ---------------------------------------------------------------------------------------------- #
# Script to test Urban Management Strategies: Presence and Abundances versus Garden Types 
# Author: [YM, MR]
# Inputs: 
#   - Data from uBMS license agreement 
# Outputs:
#   - XX
## ---------------------------------------------------------------------------------------------- #

library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggbreak)
library(glmmTMB)
library(MuMIn)
library(corrplot) 
library(psych)
library(performance)
library(tidyverse)
library(yhat)

# -----------------------------
# Load and prepare data
# -----------------------------
#MAYBE SHOUD I USED THE NEW CALCULATED ONES WITH MARTA...!
ubmsdata <- read.table("DATA/SINDEX_ZEROS.txt", header = TRUE, sep = "\t") # we should use from 2020?
ubmsdata <- ubmsdata %>%
  filter(YEAR > 2019)
space <- read.table("DATA/space_data.txt", header = TRUE, sep = "\t", fill = TRUE)[, -2]
clusters <- read.csv2("DATA/Species_clusters.csv")

# Standardize variable formats
space$Type <- factor(space$Type, levels = c("BIO", "PATRIMONIAL", "USOS"), 
                     labels = c("Biodiversity", "Historical", "Social"))

ubmsdata <- ubmsdata %>%
  left_join(clusters, by = "SPECIES") %>%
  filter(!is.na(Cluster)) %>%
  left_join(space, by = "SITE_ID") %>%
  filter(!is.na(Type)) %>%
  mutate(
    presence = as.numeric(SINDEX > 0),
    Cluster = as.factor(Cluster)
  )

# -----------------------------
# Correlation matrix of landscape variables
# -----------------------------
p.var <- space[, 5:14]
M <- corr.test(p.var)
corrplot(M$r, method = "number", type = "upper", sig.level = c(0.001, 0.01, 0.05), 
         insig = "label_sig", pch.cex = 0.9, pch.col = "grey20", order = "original")
corrplot(M$r, p.mat = M$p, type = "upper", sig.level = c(0.001, 0.01, 0.05), 
         insig = "label_sig", pch.cex = 0.9, pch.col = "grey", order = "original")

# -----------------------------
# PCA on landscape variables (with mean imputation for NAs)
# -----------------------------
pca_data <- p.var %>% mutate(across(everything(), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))
pca <- prcomp(pca_data, scale. = TRUE)
biplot(pca)

# -----------------------------
# Explore area and connectivivy differences by park type
# -----------------------------
ggplot(ubmsdata, aes(x = Type, y = TotalA)) +
  geom_boxplot() +
  labs(x = "Park Type", y = "Total Area") +
  theme_minimal()

summary(lm(TotalA ~ Type, data = ubmsdata))

ggplot(ubmsdata, aes(x = Type, y = C500)) +
  geom_boxplot() +
  labs(x = "Park Type", y = "C500") +
  theme_minimal()

summary(lm(C500 ~ Type, data = ubmsdata))

# -----------------------------
# Q1: Effect of park type and cluster on butterfly presence
# -----------------------------
#Include Area and Connectivity since they vary across types
model_bin <- glmmTMB(presence ~ Type * Cluster + TotalA + C500 + (1|SITE_ID) + (1|SPECIES),
                   data = ubmsdata, family = binomial)
summary(model_bin)

model_bin2 <- glmmTMB(presence ~ Type * Cluster + C500 + (1|SITE_ID) + (1|SPECIES),
                    data = ubmsdata, family = binomial)
summary(model_bin2)

model_bin3 <- glmmTMB(presence ~ Type * Cluster + (1|SITE_ID) + (1|SPECIES),
                    data = ubmsdata, family = binomial)
summary(model_bin3)

AIC(model_bin); AIC(model_bin2); AIC(model_bin3)

# R2 for models
r.squaredGLMM(model_bin)
r.squaredGLMM(model_bin2)
r.squaredGLMM(model_bin3)

# -----------------------------
# Models by Cluster - MAYBE BETTER ## for each cluster needs to be reduced? And add the Connectivity
# -----------------------------
#model_list <- ubmsdata %>%
#  split(.$Cluster) %>%
# lapply(function(df) {
#   glmmTMB(presence ~ Type * TotalA + (1 | SITE_ID) + (1 | SPECIES),
#         data = df, family = binomial)
# })

### PB is that Type and TotalA are very correlated. Pontential solution:
# model totalA and type and extract residuals: 
ubmsdata$resid_TotalA <- resid(lm(TotalA ~ Type, data = ubmsdata)) 
ubmsdata$resid_C500 <- resid(lm(C500 ~ Type, data = ubmsdata)) 

# o con ambas: 
lm_multi <- lm(cbind(TotalA, C500) ~ Type, data = ubmsdata)
resid_multi <- residuals(lm_multi)

ubmsdata$resid_TotalA <- resid_multi[,1]
ubmsdata$resid_C500   <- resid_multi[,2]

#model the residuals for our models
model_list <- ubmsdata %>%
  split(.$Cluster) %>%
  lapply(function(df) {
    glmmTMB(presence ~ Type+resid_TotalA +resid_C500+ (1 | SITE_ID) + (1 | SPECIES),
            data = df, family = binomial)
  })

summary(model_list[[1]])
summary(model_list[[2]])
summary(model_list[[3]])



# -----------------------------
# Plot proportion of presence by Type and Cluster
# -----------------------------
datos_resumen <- ubmsdata %>%
  group_by(Type, Cluster) %>%
  summarise(
    n = n(),
    positivos = sum(presence, na.rm = TRUE),
    prop = positivos / n,
    se = sqrt((prop * (1 - prop)) / n),
    .groups = "drop"
  ) %>%
  mutate(
    ic_inf = pmax(0, prop - 1.96 * se),
    ic_sup = pmin(1, prop + 1.96 * se)
  )

ggplot(datos_resumen, aes(x = Type, y = prop, color = Cluster, group = Cluster)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ic_inf, ymax = ic_sup), width = 0.2) +
  geom_line() +
  labs(x = "Park Type", y = "Proportion of presence", color = "Cluster") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "grey"),
    legend.position = "none"
  )

# -----------------------------
# Q2: Abundance model (only presence cases) - To do same as above once decided (MR, YM)
# -----------------------------
ubms_pos <- ubmsdata %>% filter(SINDEX > 0)

model_abund <- glmmTMB(
  SINDEX ~ Type * Cluster + TotalA + (1 | SITE_ID) + (1 | SPECIES),
  family = Gamma(link = "log"),
  data = ubms_pos
)

summary(model_abund)

# -----------------------------
# Plot abundance by space type and cluster
# -----------------------------
ggplot(ubms_pos, aes(x = Type, y = SINDEX, fill = Type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = Type), width = 0.2, size = 1.5, shape= 1, alpha = 0.6, show.legend = FALSE) +
  scale_y_break(c(100, 200), scales = 0.2) +
  scale_fill_manual(values = c("lightgreen", "lightyellow", "grey")) +
  scale_color_manual(values = c("lightgreen", "lightyellow", "grey")) +
  facet_wrap(~ Cluster) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "grey"),
    legend.position = "none"
  ) +
  labs(x = "Park Type", y = "Abundance")

