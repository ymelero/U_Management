library(DiagrammeR)

grViz("
digraph sem {
  graph [layout = dot, rankdir = TB]

  node [fontname = 'Calibri', fontsize = 10, color = grey50,penwidth = 0.5, fontcolor = black,
        style = rounded]

  edge [color = grey70, penwidth = 0.5, arrowsize = 0.5]

  # Nodos
  Area [label = 'Area total', shape = box]
  Conn [label = 'Connectivity (C500)', shape = box]
  Tipo [label = 'Garden type', shape = box] # anidado pq parte de la variación en el tipo de jardín está explicada por el área.

  V1 [label = 'C3-grass', shape = ellipse]
  V2 [label = 'P-grass',  shape = ellipse]
  V3 [label = 'Herb-grass', shape = ellipse]
  V4 [label = 'V-grass',  shape = ellipse]

  Resp [label = 'Presence / Abundance', shape = box]

  # Rangos (alineación horizontal)
  {rank = same; Area; Conn}
  {rank = same; Tipo}
  {rank = same; V1; V2; V3; V4}
  {rank = same; Resp}

  # Flechas Area/Conn -> Garden type
  Area -> Tipo
  Conn -> Tipo

  # Flechas exógenas -> vegetación
  Area -> V1
  Area -> V2
  Area -> V3
  Area -> V4

  #Conn -> V1
  #Conn -> V2
  #Conn -> V3
  #Conn -> V4

  Tipo -> V1
  Tipo -> V2
  Tipo -> V3
  Tipo -> V4

  # Vegetación -> respuesta
  V1 -> Resp
  V2 -> Resp
  V3 -> Resp
  V4 -> Resp

  # Directos exógenos -> respuesta
  Tipo -> Resp
  Area -> Resp
  Conn -> Resp
}
")


## SEM_B with area and garden type with co-variance: 
grViz("
digraph SEM_B {
  graph [layout = dot, rankdir = TB]

  node [fontname = 'Calibri', fontsize = 10, color = grey50, arrowsize = 0.5, fontcolor = black,
        style = rounded]
  edge [color = grey70, penwidth = 0.5, arrowsize = 0.5]

  # Nodos
  Area [label = 'Area total', shape = box]
  Conn [label = 'Connectivity (C500)', shape = box]
  Tipo [label = 'Garden type', shape = box]

  V1 [label = 'C3-grass', shape = ellipse]
  V2 [label = 'P-grass',  shape = ellipse]
  V3 [label = 'Herb-grass', shape = ellipse]
  V4 [label = 'V-grass',  shape = ellipse]

  Resp [label = 'Presence / Abundance', shape = box]

  # Rangos
  {rank = same; Area; Conn; Tipo}
  {rank = same; V1; V2; V3; V4}
  {rank = same; Resp}

  # Exógenos -> vegetación
  Area -> V1
  Area -> V2
  Area -> V3
  Area -> V4

  Tipo -> V1
  Tipo -> V2
  Tipo -> V3
  Tipo -> V4

  # Vegetación -> respuesta
  V1 -> Resp
  V2 -> Resp
  V3 -> Resp
  V4 -> Resp

  # Directos exógenos -> respuesta
  Tipo -> Resp
  Area -> Resp
  Conn -> Resp

  # Covarianza entre Area y Tipo
  Area -> Tipo [dir=both, arrowhead=none, arrowtail=none, style=dashed, color=grey70]
}
")
