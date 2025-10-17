# ---------------------------------------------------------------------------------------------- #
# Rao’s Quadratic Entropy (RaoQ) per SITE × YEAR
# Author: [YM]
# ---------------------------------------------------------------------------------------------- #


library(tidyverse)
library(FD)

#--- 0. Load and prepare data
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

list(unique(ubmsdata$SPECIES))

traits <- read.csv2("DATA/Traits.csv")[,1:9]
traits <- traits %>%
  rename("SPECIES" = Scientific.Name)

#ubmsdata <- ubmsdata %>%
#  left_join(traits, by = "SPECIES") 
#na.ubmsdata<-ubmsdata %>% filter(is.na(HSI))


# --- 1. Abundance matrix (SITE × SPECIES) 
abun <- ubmsdata %>%
  group_by(SITE_ID, YEAR, SPECIES) %>%
  summarise(abun = sum(SINDEX, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = SPECIES, values_from = abun, values_fill = 0) %>%
  unite("site_year", SITE_ID, YEAR, sep = "_") %>%
  arrange(site_year) %>%                # opcional, solo para orden
  as.data.frame() %>%                   # <- convertir a data.frame
  column_to_rownames("site_year") %>%   # <- ahora sí quedan como rownames
  as.matrix()

abun <- abun[rowSums(abun) > 0, , drop = FALSE]

# --- 2. Trait Matrix (SPECIES × TRAITS) 
traits.m <- traits %>%
  filter(SPECIES %in% colnames(abun)) %>%
  distinct(SPECIES, .keep_all = TRUE) %>%
  as.data.frame() %>%
  column_to_rownames("SPECIES")

# order objects
common_sp <- intersect(rownames(traits.m), colnames(abun))
common_sp <- sort(common_sp)
traits.m <- traits.m[common_sp, , drop = FALSE]
abun   <- abun[, common_sp, drop = FALSE]

# clean DBs
abun.clean <- abun[rowSums(abun) > 0, , drop = FALSE]
abun.clean <- abun.clean[, colSums(abun.clean) > 0, drop = FALSE]
traits.m <- traits.m[rownames(traits.m) %in% colnames(abun.clean), , drop = FALSE]
traits.m  <- traits.m[,-1]
traits.m$SSI <- as.numeric(traits.m$SSI); traits.m$HSI <- as.factor(traits.m$HSI)
traits.m$Voltinism <- as.factor(traits.m$Voltinism); traits.m$Overwintering <- as.factor(traits.m$Overwintering)
traits.m$STI <- as.numeric(traits.m$STI); traits.m$Mobility <- as.factor(traits.m$Mobility)
traits.m$TAO <- as.numeric(traits.m$TAO)

# --- 3. Functional Distance (Gower: categorical and continous) & Rao's Q
stopifnot(identical(rownames(traits.m), colnames(abun.clean)))
gowerD <- gowdis(traits.m)   

rao <- dbFD(
  x = traits.m,
  a = abun.clean,
  calc.FRic = FALSE,
  calc.FDiv = FALSE,
  calc.CWM  = FALSE,
  stand.x   = TRUE,
  messages  = FALSE
)

fd.df <- tibble(
  site_year = names(rao$RaoQ),
  nbsp  = as.numeric(rao$nbsp),
  sing  = as.numeric(rao$sing.sp),
  FEve  = as.numeric(rao$FEve),
  FDis  = as.numeric(rao$FDis),
  RaoQ  = as.numeric(rao$RaoQ)
) %>%
  separate(site_year, into = c("SITE_ID", "YEAR"), sep = "_") %>%
  mutate(YEAR = as.numeric(YEAR)) %>%
  relocate(SITE_ID, YEAR)
