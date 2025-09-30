# ---------------------------------------------------------------------------------------------- #
# SEM Diagram with brms results: only significant paths (95% CI does not cross 0)
# Author: YM
# ---------------------------------------------------------------------------------------------- #

library(DiagrammeR)
library(glue)
library(dplyr)

# 1. Extract coefficients for each submodel
coef_C3   <- fixef(fit_sem, resp = "z_C3_beta")
coef_P    <- fixef(fit_sem, resp = "z_P_beta")
coef_Herb <- fixef(fit_sem, resp = "z_Herb_beta")
coef_Viv  <- fixef(fit_sem, resp = "z_Viv_beta")
coef_pres <- fixef(fit_sem, resp = "presence")

# 2. Collect effects of interest, including CI
edges <- tibble(
  from  = c("Area","Area","Area","Area","Tipo","Tipo","Tipo","Tipo",
            "V1","V2","V3","V4","Tipo","Area","Conn"),
  to    = c("V1","V2","V3","V4","V1","V2","V3","V4",
            "Resp","Resp","Resp","Resp","Resp","Resp","Resp"),
  est   = c(
    coef_C3["z_Area","Estimate"],
    coef_P["z_Area","Estimate"],
    coef_Herb["z_Area","Estimate"],
    coef_Viv["z_Area","Estimate"],
    coef_C3[grep("Type", rownames(coef_C3))[1],"Estimate"],
    coef_P[grep("Type", rownames(coef_P))[1],"Estimate"],
    coef_Herb[grep("Type", rownames(coef_Herb))[1],"Estimate"],
    coef_Viv[grep("Type", rownames(coef_Viv))[1],"Estimate"],
    coef_pres["z_C3_beta","Estimate"],
    coef_pres["z_P_beta","Estimate"],
    coef_pres["z_Herb_beta","Estimate"],
    coef_pres["z_Viv_beta","Estimate"],
    coef_pres[grep("Type", rownames(coef_pres))[1],"Estimate"],
    coef_pres["z_Area","Estimate"],
    coef_pres["z_Conn","Estimate"]
  ),
  l95   = c(
    coef_C3["z_Area","Q2.5"],
    coef_P["z_Area","Q2.5"],
    coef_Herb["z_Area","Q2.5"],
    coef_Viv["z_Area","Q2.5"],
    coef_C3[grep("Type", rownames(coef_C3))[1],"Q2.5"],
    coef_P[grep("Type", rownames(coef_P))[1],"Q2.5"],
    coef_Herb[grep("Type", rownames(coef_Herb))[1],"Q2.5"],
    coef_Viv[grep("Type", rownames(coef_Viv))[1],"Q2.5"],
    coef_pres["z_C3_beta","Q2.5"],
    coef_pres["z_P_beta","Q2.5"],
    coef_pres["z_Herb_beta","Q2.5"],
    coef_pres["z_Viv_beta","Q2.5"],
    coef_pres[grep("Type", rownames(coef_pres))[1],"Q2.5"],
    coef_pres["z_Area","Q2.5"],
    coef_pres["z_Conn","Q2.5"]
  ),
  u95   = c(
    coef_C3["z_Area","Q97.5"],
    coef_P["z_Area","Q97.5"],
    coef_Herb["z_Area","Q97.5"],
    coef_Viv["z_Area","Q97.5"],
    coef_C3[grep("Type", rownames(coef_C3))[1],"Q97.5"],
    coef_P[grep("Type", rownames(coef_P))[1],"Q97.5"],
    coef_Herb[grep("Type", rownames(coef_Herb))[1],"Q97.5"],
    coef_Viv[grep("Type", rownames(coef_Viv))[1],"Q97.5"],
    coef_pres["z_C3_beta","Q97.5"],
    coef_pres["z_P_beta","Q97.5"],
    coef_pres["z_Herb_beta","Q97.5"],
    coef_pres["z_Viv_beta","Q97.5"],
    coef_pres[grep("Type", rownames(coef_pres))[1],"Q97.5"],
    coef_pres["z_Area","Q97.5"],
    coef_pres["z_Conn","Q97.5"]
  )
)

# 3. Keep only edges where CI excludes zero
edges <- edges %>%
  filter(!(l95 < 0 & u95 > 0)) %>%
  mutate(
    width = pmax(0.5, abs(est)*3),
    color = ifelse(est < 0, "red3", "grey30"),
    label = paste0(round(est,2), " [", round(l95,2), ", ", round(u95,2), "]")
  )

# 4. Build DOT code
edge_lines <- glue("{from} -> {to} [label = '{label}', color = '{color}', penwidth = {width}]",
                   from = edges$from, to = edges$to,
                   label = edges$label, color = edges$color, width = edges$width)

# Add covariance Area <-> Tipo (always drawn, dashed)
edge_lines <- c(edge_lines,
                "Area -> Tipo [dir=both, arrowhead=none, arrowtail=none, style=dashed, color=grey70]")

# 5. Render the diagram
grViz(glue("
digraph SEM_B {
  graph [layout = dot, rankdir = TB]

  # Node style
  node [fontname = 'Calibri', fontsize = 10, color = grey50, fontcolor = black,
        style = rounded]
  # Edge style
  edge [fontname = 'Calibri', fontsize = 9, arrowsize = 0.6]

  # Nodes
  Area [label = 'Area total', shape = box]
  Conn [label = 'Connectivity (C500)', shape = box]
  Tipo [label = 'Garden type', shape = box]

  V1 [label = 'C3-grass']
  V2 [label = 'P-grass']
  V3 [label = 'Herb-grass']
  V4 [label = 'V-grass']

  Resp [label = 'Presence', shape = box]

  # Horizontal ranks
  {rank = same; Area; Conn; Tipo}
  {rank = same; V1; V2; V3; V4}
  {rank = same; Resp}

  # Edges
  {paste(edge_lines, collapse = '\n  ')}
}
"))
