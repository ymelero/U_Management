# ---------------------------------------------
# NDVI Barcelona - ggplot version
# ---------------------------------------------

library(sf)
library(dplyr)
library(raster)      # or terra
library(viridis)
library(ggplot2)
library(readr)

# 1. Load ICGC municipal boundaries (line version)
gdf <- st_read("DATA/mapa-municipal-estat-v1r0-20250701.shp")

# 2. Filter only Barcelona municipality (INE code = 080193)
bcn_lines <- gdf %>%
  filter(CODIMUNI1 == "080193" | CODIMUNI2 == "080193")

# 3. Merge and polygonize lines
bcn_poly <- st_union(bcn_lines) |> st_polygonize()

# 4. Convert to valid MULTIPOLYGON and add attribute
bcn_multipoly <- st_collection_extract(bcn_poly, "POLYGON")
bcn_sf <- st_sf(NOMMUNI = "Barcelona", geometry = bcn_multipoly)

# 5. Load NDVI raster
ndvi <- raster("DATA/NDVI_2019.tif")

# 6. Make sure polygon and raster use the same CRS
bcn_sf <- st_transform(bcn_sf, crs(ndvi))

# 7. Clip raster with Barcelona polygon
ndvi_mask <- mask(crop(ndvi, bcn_sf), bcn_sf)

# 8. Aggregate raster to reduce resolution (avoid memory issues)
ndvi_lowres <- aggregate(ndvi_mask, fact = 10, fun = mean)

# 9. Convert raster to dataframe for ggplot
ndvi_df <- as.data.frame(ndvi_lowres, xy = TRUE)
colnames(ndvi_df)[3] <- "NDVI"

# 1. Read CSV with coordinates (assume columns are lon, lat)
sites <- read_delim("DATA/ubms_sites.csv", delim = ";") %>% filter (!is.na(transect_id)) %>% 
  filter(!grepl("_A", transect_id))

# 2. Convert to sf object (EPSG:4326 = lon/lat WGS84)
head(sites)
sites_sf <- st_as_sf(sites, coords = c("transect_longitude", "transect_latitude"), crs = 4326)

# 3. Transform to same CRS as raster and polygon (EPSG:25831)
sites_sf <- st_transform(sites_sf, crs(bcn_sf))
# Keep only sites inside Barcelona polygon
sites_sf <- st_intersection(sites_sf, bcn_sf)

# 10. Plot with ggplot2
ggplot() +
  geom_raster(data = ndvi_df, aes(x = x, y = y, fill = NDVI)) +
  scale_fill_gradientn(
    colours = c("peru", "peru",   "lightgreen","green3", "darkgreen"),
    values = scales::rescale(c(-0.5, 0, 0.2, 0.6, 1)),
    na.value = "transparent"
  ) +
  geom_sf(data = bcn_sf, fill = NA, color = "black", size = 0.3) +
  geom_sf(data = sites_sf, color = "red", size = 2, shape = 21, fill = "yellow") +
  theme_minimal() +
  labs(fill = "NDVI")+
  theme(
    axis.title = element_blank(),   # remove axis labels
    axis.text  = element_blank(),   # remove axis tick labels
    axis.ticks = element_blank(),   # remove ticks
    panel.grid = element_blank()    # remove grid
  )

# or
ggplot() +
  geom_raster(data = ndvi_df, aes(x = x, y = y, fill = NDVI)) +
  scale_fill_viridis_c(option = "C", na.value = "transparent") +  # viridis!
  geom_sf(data = bcn_sf, fill = NA, color = "black", size = 0.3) +
  geom_sf(data = sites_sf, color = "purple4", size = 2, shape = 21, fill = "purple4") +
  theme_minimal() +
  labs(fill = "NDVI")+
  theme(
    axis.title = element_blank(),   # remove axis labels
    axis.text  = element_blank(),   # remove axis tick labels
    axis.ticks = element_blank(),   # remove ticks
    panel.grid = element_blank()    # remove grid
  )

# or
ggplot() +
  geom_raster(data = ndvi_df, aes(x = x, y = y, fill = NDVI)) +
  scale_fill_viridis_c(option = "C", na.value = "transparent") +  # viridis!
  geom_sf(data = bcn_sf, fill = NA, color = "black", size = 0.3) +
  geom_sf(data = sites_sf, color = "grey", size = 3, shape = 21, fill = "yellow2") +
  theme_minimal() +
  labs(fill = "NDVI")+
  theme(
    axis.title = element_blank(),   # remove axis labels
    axis.text  = element_blank(),   # remove axis tick labels
    axis.ticks = element_blank(),   # remove ticks
    panel.grid = element_blank()    # remove grid
  )
