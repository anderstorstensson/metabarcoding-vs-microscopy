library(tidyverse)
library(ggOceanMaps)
library(ggOceanMapsData)
library(sf)
library(viridis)
library(inlmisc)

# Define file paths
station_file <- 'data/everything_from_shark/stations_to_extract.txt'
shark_salinity_file <- "data/shapefiles/sharkweb_data - 2023-05-22T130134.610.txt" # Data between 1980-2023
basins_file <- 'data/shapefiles/sharkweb_shapefiles/Havsomr_SVAR_2016_3b_CP1252.shp'

# Read files
stations <- read.table(station_file, sep = '\t', header = TRUE)
shark_salinity <- read.table(shark_salinity_file, sep = '\t', header = TRUE)
basins <- st_read(basins_file)

# Set CRS
basins <- st_set_crs(basins, 3006)

# Aggregate basins
all_basins <- basins %>%
  group_by(BASIN_NR) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup()

# Change CRS
all_basins <- st_transform(all_basins, 4326)

# Define major basins
basin_vec <- c("I", rep("II", 2), rep("III", 2), "X", rep("III", 9), "IV", "V")

# Remove basins that are not defined as Baltic proper
basin_vec[c(10, 14, 15, 4, 5, 11)] <- c("XX", "XXX", "XXXX", "XXXXX", "XXXXX", "XXXXXX")

# Add sea basin name
all_basins$sea_basin <- basin_vec

# Aggregate major basins
major_basins <- all_basins %>%
  group_by(sea_basin) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup()

# get borders
borders = st_intersection(major_basins) %>%
  filter(n.overlaps == 2)

# Change name of B7
stations$station_plots <- gsub("B7", "B7 / B3", stations$station_plots)

# Create basemap
baltic_sea_map <- basemap(
  limits = c(min(stations$lon) - 1, max(stations$lon) + 1, min(stations$lat) - 1, max(stations$lat) + 1),
  land.col = "#eeeac4",
  land.border.col = "black",
  rotate = TRUE,
  bathymetry = FALSE
)

# Wrangle salinity data
smhi_salinity <- shark_salinity %>%
  filter(Provtagningsdjup..m. <= 10) %>%
  dplyr::select(
    Stationsnamn,
    År,
    Månad,
    Provets.longitud..DD.,
    Provets.latitud..DD.,
    Provtagningsdjup..m.,
    Mätvärde
  ) %>%
  filter(!is.na(Mätvärde)) %>%
  dplyr::group_by(Provets.longitud..DD., Provets.latitud..DD.) %>%
  dplyr::summarise(
    salinity = mean(Mätvärde),
    lat = mean(Provets.latitud..DD.),
    lon = mean(Provets.longitud..DD.)
  )

# Convert points to sf
points_sf <- st_as_sf(smhi_salinity, coords = c("lon", "lat"), crs = st_crs(all_basins))

# Join salinity data with ICES area
shark_areas_salinity <- st_join(points_sf, all_basins)

salinity_by_area <- smhi_salinity %>% 
  left_join(shark_areas_salinity) %>%
  dplyr::group_by(BASIN_NR) %>%
  summarise(salinity = mean(salinity)) %>%
  filter(!is.na(BASIN_NR)) 

# write.table(salinity_by_area, "data/shapefiles/salinty_by_area.txt", sep = "\t", row.names = FALSE)

# Join salinity data with ICES area
shark_areas <- all_basins %>% 
  left_join(salinity_by_area)

# Copy value within the Åland sea as we do not have any data from eastern Åland sea
shark_areas$salinity[5] <- shark_areas$salinity[4]

# Add sea basins to the map
map <- baltic_sea_map + 
  geom_sf(data = st_difference(st_buffer(shark_areas, dist = 1000, nQuadSegs = 1)), aes(fill = salinity), colour = NA, show.legend = TRUE) +
  geom_sf(data = major_basins, fill = NA, colour = NA) +
  geom_sf(data = st_cast(borders, "LINESTRING"), colour = "black", lty = "solid", linewidth = 0.1, lineend = "round")+
  scale_size_identity()

# Make the graticules
lims <- attributes(map)$limits 
graticule <- sf::st_graticule(
  c(lims[1], lims[3], lims[2], lims[4]), 
  crs = attributes(map)$proj,
  lon = seq(-180, 180, 15),
  lat = seq(-90, 90, 5)
)

# Define dimensions
lon <- c(9.5, 11.3, 19.1, 19.0, 22.6)
lat <- c(58.1, 56.9, 56.3, 61.6, 64.5)
cord <- data.frame(lon, lat)
cord_sf <- st_as_sf(cord, coords = c("lon", "lat")) 
cord_sf_4326 <- st_set_crs(cord_sf, 4326)
cord_sf_4326$sea_basin <- c("V", "IV", "III", "II", "I")

# Define a biased color scheme
cols <- GetColors(256, scheme = "viridis", bias = 2, alpha = 0.8) 

# Reorder map and add labels
reordered_map <- reorder_layers(map) +
  geom_sf(data = graticule, color = "transparent", size = LS(0.5)) +
  coord_sf(
    xlim = lims[1:2],
    ylim = lims[3:4],
    crs = attributes(map)$proj
  ) +
  geom_spatial_text_repel(
    data = stations,
    aes(x = lon, y = lat, label = station_plots),
    size = 5,
    fontface = "bold",
    color = "white",
    bg.color = "black",
    bg.r = .1
  ) +
  geom_label(
    data = cord_sf_4326,
    aes(label = sea_basin, geometry = geometry),
    stat = "sf_coordinates",
    alpha = 0.7
  ) +
  # scale_fill_viridis(
  #   option = "viridis",
  #   name = "Salinity\n(psu)",
  #   na.value = "white",
  #   alpha = 0.5
  # ) +
  scale_fill_gradientn(colours = alpha(cols, 0.5),
                       name = "Salinity\n(psu)",
                       na.value = "white",
                       # limits = c(0, 30),
                       breaks = c(5, 10, 15, 20, 25)
  ) +
  theme(
    panel.grid.major = element_line(colour = "transparent"),
    panel.grid.minor = element_line(colour = "transparent")
  ) +
  xlab("Longitude") +
  ylab("Latitude") +
  geom_spatial_point(
    data = stations,
    aes(x = lon, y = lat),
    pch = 21,
    size = 2,
    fill = "white",
    colour = "black"
  ) +
  # labs(caption = "I: Bothnian bay, II: Bothnian sea, III: Baltic proper, IV: Kattegat, V: Skagerrak") +
  theme(plot.title = element_text(hjust = 0.5))

# Save map object
save(reordered_map, file = "plots/microscopy_vs_barcoding/map.RData")

# Save map
ggsave(
  plot = reordered_map,
  path = "plots/microscopy_vs_barcoding/",
  filename = "fig_1_map2.png",
  device = "png",
  dpi = 300,
  width = 7,
  height = 7,
  bg = "white"
)

ggsave(
  plot = reordered_map,
  path = "plots/microscopy_vs_barcoding/pdf",
  filename = "fig_1_map.pdf",
  device = "pdf",
  width = 7,
  height = 7,
  bg = "white"
)

ggsave(
  plot = reordered_map,
  path = "plots/microscopy_vs_barcoding/eps",
  filename = "fig_1_map.eps",
  device = "eps",
  width = 7,
  height = 7,
  bg = "white"
)
