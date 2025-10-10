# ---------------------------------------------------------------------------------------------- #
# Raw / Proportional Data Exploratory Plots
# Author: [YM]
# ---------------------------------------------------------------------------------------------- #


library(tidyverse)
library(ggplot2)
library(gridExtra)

space    <- read.delim("DATA/space_data.txt", sep = "\t")
space <- space %>%
  mutate(
    C3.proportion   = C3ha / TotalA, # ornamental veg., hugh hyric resources, highly managed (sega y riego)
    P.proportion    = PUha / TotalA, # native veg. with spontaneous herbs and flowers, low management(one sega por año??)
    Herb.proportion = pmin(HERbha / TotalA, 1), # native herbs and flowers managed (water) but no siega (only once?? - IN Biod types)
    Viv.proportion  = VIVha / TotalA # ??
  )

# --- 1. Space composition versus Type
# Reshape to long format for easy plotting
space_long <- space %>%
  pivot_longer(cols = c(C3.proportion, P.proportion, Herb.proportion, Viv.proportion),
               names_to = "Group", values_to = "Proportion")

# Clean group names
space_long$Group <- factor(space_long$Group,
                           levels = c("C3.proportion", "P.proportion", "Herb.proportion", "Viv.proportion"),
                           labels = c("C3", "P", "Herb", "Viv"))

# --- Boxplot (overall) ---
ggplot(space_long, aes(x = Group, y = Proportion, fill = Group)) +
  geom_boxplot(alpha = 0.6, width = 0.7, outlier.shape = 21) +
  scale_fill_brewer(palette = "Set2") + ylim(0, 0.5) +
  theme_classic(base_size = 12) +
  labs(x = "Functional group", y = "Proportion"
       ) +
  theme(legend.position = "none",
        plot.title = element_text(size = 13, face = "bold"))

# --- Boxplot by Type ---
p1<-ggplot(space_long, aes(x = Type, y = Proportion, fill = Group)) +
  geom_boxplot(alpha = 0.6, width = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_brewer(palette = "Set2") +
  theme_classic(base_size = 12) + #ylim(0, 0.5) +
  labs(x = "", y = "Proportion") +
  theme(legend.title = element_blank(),
        legend.position = "top") + theme(axis.text.x = element_blank(),
           axis.ticks.x = element_blank())

# --- remove C3 for better visulisation
space_long_red <- space_long %>% filter(Group !="C3")
cols <- RColorBrewer::brewer.pal(4, "Set2")[-1]
p2<-ggplot(space_long_red, aes(x = Type, y = Proportion, fill = Group)) +
  geom_boxplot(alpha = 0.6, width = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = cols) +  
  theme_classic(base_size = 12) + ylim(0, 0.15) +
  labs(x = "Space type", y = "Proportion") +
  theme(legend.position = "none")

grid.arrange(p1, p2, ncol = 1, nrow = 2)

# --- 2. Area & Connectivity versus Type
p_area <- ggplot(space, aes(x = Type, y = TotalA, fill = Type)) +
  geom_boxplot(alpha = 0.6, width = 0.7) +
  scale_fill_brewer(palette = "Set2") +
  theme_classic(base_size = 12) +
  labs(x = "", y = "Area") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())

p_con <- ggplot(space, aes(x = Type, y = C500, fill = Type)) +
  geom_boxplot(alpha = 0.6, width = 0.7) +
  scale_fill_brewer(palette = "Set2") +
  theme_classic(base_size = 12) +
  labs(x = "Space type", y = "Connectivity") +
  theme(legend.position = "none")

grid.arrange(p_area, p_con, ncol = 1, nrow = 2)

# --- 3. Plot variables correlations
library(corrplot) 
library(psych)
p.var <- space[, c(4:9,24:27)]
p.var <- p.var %>% rename("Open Area" = OpenA,
                          "Closed Area" = ClosedA,
                          "Total Area" = TotalA,
                          "Conn. 200m" = C200,
                          "Conn. 500m" = C500,
                          "Conn. 1km" = C1k,
                          "C3 Proportion" = C3.proportion,
                          "P Proportion" = P.proportion,
                          "Herb. Proportion" = Herb.proportion,
                          "Viv. Proportion" = Viv.proportion
                          )
M <- corr.test(p.var)
corrplot(M$r, method = "number", type = "upper", sig.level = c(0.001, 0.01, 0.05), 
         insig = "label_sig", pch.cex = 0.9, pch.col = "grey20", order = "original")
corrplot(M$r, p.mat = M$p, type = "upper", sig.level = c(0.001, 0.01, 0.05), 
         insig = "label_sig", pch.cex = 0.9, pch.col = "grey", order = "original", tl.col = "black")

# --- 4. Plot raw presence and abundance per cluster
# read ubms and cluster data from 1_Model_Bayesian....
p_abundance <- ggplot(ubmsdata, aes(x = Cluster, y = SINDEX, fill = Cluster)) +
  geom_boxplot(alpha = 0.6, width = 0.7) +
  scale_fill_brewer(palette = "Set2") +
  theme_classic(base_size = 12) +
  labs(x = "", y = "Presence / Abundance") +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())

p_abundance_zoom <- ggplot(ubmsdata, aes(x = Cluster, y = SINDEX, fill = Cluster)) +
  geom_boxplot(alpha = 0.6, width = 0.7) +
  scale_fill_brewer(palette = "Set2") +
  theme_classic(base_size = 12) + ylim(10, 100) +
  labs(x = "", y = "") +
  theme(legend.position = "none") 

grid.arrange(p_abundance, p_abundance_zoom, ncol = 1, nrow = 2)
