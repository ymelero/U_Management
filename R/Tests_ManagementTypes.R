# ---------------------------------------------------------------------------------------------- #
# Script to test Urban Management Strategies: Presence and Abundances versus management types 
# Author: [YM, MR]
# Inputs: 
#   - Data from uBMS license agreement 
# Outputs:
#   - XX
#  ---------------------------------------------------------------------------------------------- #

library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggbreak)
library(glmmTMB)
library(MuMIn)
library(corrplot) 
library(psych)
library(ggeffects)

# -----------------------------
# Load and prepare data
# -----------------------------
ubmsdata <- read.table("DATA/SINDEX_ZEROS.txt", header = TRUE, sep = "\t")
space <- read.table("DATA/space_data.txt", header = TRUE, sep = "\t", fill = TRUE)[, -2]
clusters <- read.csv2("DATA/Species_clusters.csv")

# "Since we are interested in the effect of each strategy, without accounting for the Area occupied by each one 
# -----------------------------
# Calculate proportion of area of each strategy type. 
# -----------------------------
space <- space %>%
  mutate(
    C3.proportion = C3ha / TotalA,
    P.proportion = PUha / TotalA,
    Herb.proportion = HERbha / TotalA,
    Viv.proportion = VIVha / TotalA,
  )

ubmsdata <- ubmsdata %>%
  left_join(clusters, by = "SPECIES") %>%
  filter(!is.na(Cluster)) %>%
  left_join(space, by = "SITE_ID") %>%
  filter(!is.na(C3.proportion)) %>%
  mutate(
    presence = as.numeric(SINDEX > 0),
    Cluster = as.factor(Cluster)
  )

# -----------------------------
# Correlation matrix of landscape variables
# -----------------------------
p.var <- space[, c(4:10,25:28)]
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
# Q3: Effect of management type and cluster on butterfly presence
# -----------------------------

# -----------------------------
# Models by Cluster - EASIER TO INTERPRET. MR ok? Accounts for the % of each management type, independent of TotalA.
# -----------------------------
model_list <- ubmsdata %>%
  split(.$Cluster) %>%
  lapply(function(df) {
    glmmTMB(presence ~ C3.proportion*TotalA + P.proportion*TotalA+
              Herb.proportion*TotalA + Viv.proportion*TotalA + C500+
              (1 | SITE_ID) + (1 | SPECIES),
            data = df, family = binomial)
  })

summary(model_list[[1]])
summary(model_list[[2]])
summary(model_list[[3]])
library(performance) # MUMIN does not work. Why?
r2(model_list[[1]])
r2(model_list[[2]]) 
r2(model_list[[1]]) # high sp variabiity

# -----------------------------
# Plot proportion of presence by Management Type and Cluster. !!! PLOT TO BE UPDATED
# -----------------------------
pred_list <- list(
  C3 = ggpredict(model_list[[1]], terms = "C3.proportion [0:0.4 by=0.001]"),
  P   = ggpredict(model_list[[1]], terms = "P.proportion [0:0.15 by=0.01]"),
  Herb = ggpredict(model_list[[1]], terms = "Herb.proportion [0:0.15 by=0.01]"),
  Viv  = ggpredict(model_list[[1]], terms = "Viv.proportion [0:0.4 by=0.01]")
)

pred_list2 <- list(
  C3 = ggpredict(model_list[[2]], terms = "C3.proportion [0:0.4 by=0.001]"),
  P   = ggpredict(model_list[[2]], terms = "P.proportion [0:0.15 by=0.01]"),
  Herb = ggpredict(model_list[[2]], terms = "Herb.proportion [0:0.15 by=0.01]"),
  Viv  = ggpredict(model_list[[2]], terms = "Viv.proportion [0:0.4 by=0.01]")
)

pred_list3 <- list(
  C3 = ggpredict(model_list[[3]], terms = "C3.proportion [0:0.4 by=0.001]"),
  P   = ggpredict(model_list[[3]], terms = "P.proportion [0:0.15 by=0.01]"),
  Herb = ggpredict(model_list[[3]], terms = "Herb.proportion [0:0.15 by=0.01]"),
  Viv  = ggpredict(model_list[[3]], terms = "Viv.proportion [0:0.4 by=0.01]")
)

pred_list_to_df <- function(pred_list, cluster_name) {
  bind_rows(
    lapply(names(pred_list), function(var) {
      df <- pred_list[[var]]
      df$Variable <- var
      df$Cluster <- cluster_name
      df
    })
  )
}

df1 <- pred_list_to_df(pred_list, "2")
df2 <- pred_list_to_df(pred_list2, "3")
df3 <- pred_list_to_df(pred_list3, "4")

predictions_all <- data.frame(bind_rows(df1, df2, df3))
predictions_all$Cluster <- as.factor(predictions_all$Cluster)

library(ggpubr)
# TO DO: Mejroar la estética de todos los plots!
vars <- c("C3.proportion", "P.proportion", "Herb.proportion", "Viv.proportion")
names(vars) <- c("C3", "P", "Herb", "Viv")

cluster_colors <- c("2" = "#56B4E9", "3" = "#009E73", "4" = "#E69F00")  # azul, verde, naranja pastel

plots_by_cluster <- list()

for (i in seq_along(levels(ubmsdata$Cluster))) {
  cl <- levels(ubmsdata$Cluster)[i]
  p_list <- list()
  
  for (j in seq_along(names(vars))) {
    v <- names(vars)[j]
    
    pred <- predictions_all[predictions_all$Cluster == cl & predictions_all$Variable == v, ]
    data_cl <- ubmsdata[ubmsdata$Cluster == cl, ]

    data_cl$presence[data_cl$presence > 1] <- 1
    data_cl$presence[data_cl$presence < 0] <- 0
    
    p <- ggplot() +
      geom_jitter(data = data_cl, aes_string(x = vars[v],
                                             y = "jitter(presence, amount = 0.08)"),
                  color = "lightgrey", alpha = 0.1, size = 0.6, width = 0.02) +
      geom_ribbon(data = pred, aes(x = x, ymin = conf.low, ymax = conf.high),
                  fill = cluster_colors[cl], alpha = 0.2) +
      geom_line(data = pred, aes(x = x, y = predicted),
                color = cluster_colors[cl], linewidth = 1) +
      labs(x = if (i == length(levels(ubmsdata$Cluster))) paste(v, "proportion") else NULL,
           y = if (j == 1) "Presence" else NULL) +
      coord_cartesian(ylim = c(0, 1)) +
      theme_minimal(base_size = 12) +
      theme(panel.grid = element_blank(),
            axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10),
            axis.text.x = if (i == length(levels(ubmsdata$Cluster))) element_text() else element_blank(),
            axis.text.y = if (j == 1) element_text() else element_blank(),
            plot.title = element_blank())
    
    p_list[[v]] <- p
  }
  
  plots_by_cluster[[cl]] <- ggarrange(plotlist = p_list, ncol = 4)
}

ggarrange(plotlist = plots_by_cluster, ncol = 1)

# -----------------------------
# Q4: Abundance model (only presence cases) - To do same as above once decided (MR, YM)
# -----------------------------
ubms_pos <- ubmsdata %>% filter(SINDEX > 0)

model_N <-ubms_pos %>%
  split(.$Cluster) %>%
  lapply(function(df) {
    glmmTMB(SINDEX ~ C3.proportion + P.proportion+
              Herb.proportion + Viv.proportion + C200 + TotalA +# !! MAYBE ADD AREA FOR ALL MODELS?? APENAS CAMBIA
              (1 | SITE_ID) + (1 | SPECIES),
            data = df, family =  Gamma(link = "log"))
  })

summary(model_N[[1]])
summary(model_N[[2]])
summary(model_N[[3]])
r2(model_list[[1]])
r2(model_list[[2]]) 
r2(model_list[[1]]) # high sp variabiity

# -----------------------------
# Plot abundance by variable and cluster
# -----------------------------
