# Resulting SEMs per cluster
# Author: [YM]
# ---------------------------------------------------------------------------------------------- #

library(DiagrammeR)

#--- 1. SEM C2 
c2 <- grViz("
digraph SEM_B {
  graph [layout = dot, rankdir = TB]

  node [fontname = 'Calibri', fontsize = 10, color = grey50, penwidth = 0.5,
        fontcolor = black, style = rounded]
  edge [color = grey70, penwidth = 0.5, arrowsize = 0.5]

  # Nodos
  Area [label = 'Area total', shape = box,style = 'rounded,filled', colorscheme=gnbu3, fillcolor=1]
  Conn [label = 'Connectivity (C500)', shape = box, style = 'rounded,filled', colorscheme=gnbu3, fillcolor=1]
  Tipo [label = 'Garden type', shape = box, style = 'rounded,filled', colorscheme=gnbu3, fillcolor=1]
  graph [nodesep = 0.2]  # hacer que nodos del mismo nivel estén más juntos

{ rank = same; Area; Conn; Tipo }

Area -> Conn [style = invis, constraint = false, weight = 50]

  V1 [label = 'C3-grass', shape = ellipse, style = filled, fillcolor = peachpuff]
  V2 [label = 'P-grass',  shape = ellipse, style = filled, colorscheme=gnbu3, fillcolor=1]
  V3 [label = 'Herb-grass', shape = ellipse]
  V4 [label = 'V-grass',  shape = ellipse]

  Resp [label = 'Presence / Abundance', shape = box]

  # Rangos
  {rank = same; Area; Conn; Tipo}
  {rank = same; V1; V2; V3; V4}
  {rank = same; Resp}

  # Exógenos -> vegetación
  Area -> V1 [color=none, style=invis, arrowsize=0]   # C3-grass
  Area -> V2
  Area -> V3 [color=none, style=invis, arrowsize=0]   # Herb-grass
  Area -> V4

  Tipo -> V1
  Tipo -> V2
  Tipo -> V3
  Tipo -> V4

  # Vegetación -> respuesta
  V1 -> Resp #[color=lightcoral, penwidth=1.2]
  V2 -> Resp #[color=palegreen3, penwidth=1.2]
  V3 -> Resp [color=none, style=invis, arrowsize=0]   # Herb-grass
  V4 -> Resp [color=none, style=invis, arrowsize=0]   # V-grass

  # Directos exógenos -> respuesta
  Tipo -> Resp
  Area -> Resp [color=none, style=invis, arrowsize=0] # eliminada
  Conn -> Resp

  # Covarianza entre Area y Tipo
  Area -> Tipo [dir=both, arrowhead=none, arrowtail=none,
                style=dashed, color=grey70]
}
")


library(DiagrammeRsvg)
library(rsvg)
svg_code <- export_svg(c2)
rsvg_pdf(charToRaw(svg_code), file = "Results/Figures/SEM_results_C2.pdf")

# Bind with boxplot
#install.packages("pdftools")
library(magick)
library(grid)
library(gridExtra)

# Read PDFs
img1 <- image_read_pdf("Results/Figures/SEM_results_C2.pdf", density = 300)
img2 <- image_read_pdf("Results/Figures/Fig_Predicted_Type_C2.pdf", density = 300)

# Convert to rasterGrob 
c2_g1 <- rasterGrob(as.raster(img1), interpolate = TRUE)
c2_g2 <- rasterGrob(as.raster(img2), interpolate = TRUE)

# Combine
grid.arrange(c2_g1, c2_g2, nrow = 2)
grid.text("(a)", x = 0.12, y = 0.97, gp = gpar(fontsize = 10))
grid.text("(b)", x = 0.12, y = 0.47, gp = gpar(fontsize = 10))


#--- 2. SEM C3
# Resulting SEMs per cluster
# Author: [YM]
# ---------------------------------------------------------------------------------------------- #

library(DiagrammeR)

c3 <- grViz("
digraph SEM_B {
  graph [layout = dot, rankdir = TB]

  node [fontname = 'Calibri', fontsize = 10, color = grey50, penwidth = 0.5,
        fontcolor = black, style = rounded]
  edge [color = grey70, penwidth = 0.5, arrowsize = 0.5]

  # Nodos
  Area [label = 'Area total', shape = box,style = 'rounded,filled', colorscheme=gnbu3, fillcolor=1]
  Conn [label = 'Connectivity (C500)', shape = box, style = 'rounded,filled', colorscheme=gnbu3, fillcolor=1]
  Tipo [label = 'Garden type', shape = box, style = 'rounded,filled', colorscheme=gnbu3, fillcolor=1]
  graph [nodesep = 0.2]  # hacer que nodos del mismo nivel estén más juntos

{ rank = same; Area; Conn; Tipo }

Area -> Conn [style = invis, constraint = false, weight = 50]

  V1 [label = 'C3-grass', shape = ellipse, style = filled, fillcolor = peachpuff]
  V2 [label = 'P-grass',  shape = ellipse]
  V3 [label = 'Herb-grass', shape = ellipse]
  V4 [label = 'V-grass',  shape = ellipse,style = filled, fillcolor = peachpuff]

  Resp [label = 'Presence / Abundance', shape = box]

  # Rangos
  {rank = same; Area; Conn; Tipo}
  {rank = same; V1; V2; V3; V4}
  {rank = same; Resp}

  # Exógenos -> vegetación
  Area -> V1 [color=none, style=invis, arrowsize=0]   # C3-grass
  Area -> V2
  Area -> V3 [color=none, style=invis, arrowsize=0]   # Herb-grass
  Area -> V4 [color=none, style=invis, arrowsize=0]   # V-grass

  Tipo -> V1
  Tipo -> V2
  Tipo -> V3
  Tipo -> V4

  # Vegetación -> respuesta
  V1 -> Resp #[color=lightcoral, penwidth=1.2]
  V2 -> Resp [color=none, style=invis, arrowsize=0]
  V3 -> Resp [color=none, style=invis, arrowsize=0]   # Herb-grass
  V4 -> Resp # [color=none, style=invis, arrowsize=0]   # V-grass

  # Directos exógenos -> respuesta
  Tipo -> Resp
  Area -> Resp # [color=none, style=invis, arrowsize=0] # eliminada
  Conn -> Resp

  # Covarianza entre Area y Tipo
  Area -> Tipo [dir=both, arrowhead=none, arrowtail=none,
                style=dashed, color=grey70]
}
")

c3


library(DiagrammeRsvg)
library(rsvg)
svg_code <- export_svg(c3)
rsvg_pdf(charToRaw(svg_code), file = "Results/Figures/SEM_results_C3.pdf")

# Bind with boxplot
img1 <- image_read_pdf("Results/Figures/SEM_results_C3.pdf", density = 300)
img2 <- image_read_pdf("Results/Figures/Fig_Predicted_Type_C3.pdf", density = 300)

c3_g1 <- rasterGrob(as.raster(img1), interpolate = TRUE)
c3_g2 <- rasterGrob(as.raster(img2), interpolate = TRUE)

grid.arrange(c3_g1, c3_g2, nrow = 2)
grid.text("(a)", x = 0.12, y = 0.97, gp = gpar(fontsize = 10))
grid.text("(b)", x = 0.12, y = 0.47, gp = gpar(fontsize = 10))

#--- 3. Cluster C4
# Resulting SEMs per cluster
# Author: [YM]
# ---------------------------------------------------------------------------------------------- #

library(DiagrammeR)

c4 <- grViz("
digraph SEM_B {
  graph [layout = dot, rankdir = TB]

  node [fontname = 'Calibri', fontsize = 10, color = grey50, penwidth = 0.5,
        fontcolor = black, style = rounded]
  edge [color = grey70, penwidth = 0.5, arrowsize = 0.5]

  # Nodos
  Area [label = 'Area total', shape = box,style = 'rounded,filled', colorscheme=gnbu3, fillcolor=1]
  Conn [label = 'Connectivity (C500)',
      shape = box,
      style = 'rounded,filled',
      fillcolor = white]
  Tipo [label = 'Garden type', shape = box, style = 'rounded,filled', colorscheme=gnbu3, fillcolor=1]
  graph [nodesep = 0.2]  # hacer que nodos del mismo nivel estén más juntos

{ rank = same; Area; Conn; Tipo }

Area -> Conn [style = invis, constraint = false, weight = 50]

  V1 [label = 'C3-grass', shape = ellipse, style = filled, fillcolor = peachpuff]
  V2 [label = 'P-grass',  shape = ellipse]
  V3 [label = 'Herb-grass', shape = ellipse]
  V4 [label = 'V-grass',  shape = ellipse]

  Resp [label = 'Presence / Abundance', shape = box]

  # Rangos
  {rank = same; Area; Conn; Tipo}
  {rank = same; V1; V2; V3; V4}
  {rank = same; Resp}

  # Exógenos -> vegetación
  Area -> V1 [color=none, style=invis, arrowsize=0]   # C3-grass
  Area -> V2
  Area -> V3 [color=none, style=invis, arrowsize=0]   # Herb-grass
  Area -> V4 [color=none, style=invis, arrowsize=0]   # V-grass

  Tipo -> V1
  Tipo -> V2
  Tipo -> V3
  Tipo -> V4

  # Vegetación -> respuesta
  V1 -> Resp #[color=lightcoral, penwidth=1.2]
  V2 -> Resp [color=none, style=invis, arrowsize=0]
  V3 -> Resp [color=none, style=invis, arrowsize=0]   # Herb-grass
  V4 -> Resp [color=none, style=invis, arrowsize=0]   # V-grass

  # Directos exógenos -> respuesta
  Tipo -> Resp
  Area -> Resp # [color=none, style=invis, arrowsize=0] # eliminada
  Conn -> Resp [color=none, style=invis, arrowsize=0]

  # Covarianza entre Area y Tipo
  Area -> Tipo [dir=both, arrowhead=none, arrowtail=none,
                style=dashed, color=grey70]
}
")

svg_code <- export_svg(c4)
rsvg_pdf(charToRaw(svg_code), file = "Results/Figures/SEM_results_C4.pdf")

# Bind with boxplot
img1 <- image_read_pdf("Results/Figures/SEM_results_C4.pdf", density = 300)
img2 <- image_read_pdf("Results/Figures/Fig_Predicted_Type_C4.pdf", density = 300)

# Convertir a rasterGrob (para que grid lo entienda)
c4_g1 <- rasterGrob(as.raster(img1), interpolate = TRUE)
c4_g2 <- rasterGrob(as.raster(img2), interpolate = TRUE)

# Combinar los dos con grid.arrange
grid.arrange(c4_g1, c4_g2, nrow = 2)
grid.text("(a)", x = 0.12, y = 0.97, gp = gpar(fontsize = 10))
grid.text("(b)", x = 0.12, y = 0.47, gp = gpar(fontsize = 10))




