# Make a map of the 7 Inside South Coast Chum Conservation Units

library(bcmaps)
library(ggplot2)
library(sf)
library(ggspatial)
library(rgdal)
# library(rnaturalearth)
# library(rnaturalearthdata)
# library(ggmap)
#remotes::install_github("ropensci/rnaturalearthhires")
library(leaflet)
library(PNWColors)

# Download chum CU data shapefile from:
# https://open.canada.ca/data/en/dataset/f86c0867-d38d-4072-bd08-57cbbcbafa46
chum_cu <- st_read("DataIn/chum_CU_boundary_shapefile/Chum Salmon CU Boundary_En.shp")
# Dowdload DFO Fishery Management areas shapefile from:
# https://catalogue.data.gov.bc.ca/dataset/dfo-statistical-areas-boundaries#edc-pow
fma <- st_read("DataIn/dfo_fishery_mgmt_areas_shapefile/DFO_STAT_A_polygon.shp")
areas <- c(12:19,28,29) # make vector of areas used in chum analysis
fma <- fma[!fma$MNGMNTR==0, ] # remove small areas (?)
fma <- fma[fma$MNGMNTR %in% areas, ] # keep only areas used in chum analysis 
# Freshwater Atlas Large Rivers (polygons)  https://catalogue.data.gov.bc.ca/dataset/freshwater-atlas-rivers
riv <- st_read("DataIn/rivers/FWRVRSPL_polygon.shp")
# Freshwater Atlas Lakes https://catalogue.data.gov.bc.ca/dataset/freshwater-atlas-lakes#edc-pow
lakes <- st_read("DataIn/lakes/FWLKSPL_polygon.shp")
lakes <- lakes[lakes$AREA_HA >100, ] # select lakes larger than 100 ha
# get south coast chum CUs 
sc_cus <- c("Southern Coastal Streams", "Northeast Vancouver Island", "Upper Knight", "Loughborough", "Bute Inlet", "Georgia Strait", "Howe Sound-Burrard Inlet")
# Get just the 5 CUs without CU-level infilling
sc_cus_no_CU_infilled <- c("Southern Coastal Streams", "Northeast Vancouver Island", "Loughborough","Georgia Strait", "Howe Sound-Burrard Inlet")
#scc <- chum_cu[chum_cu$CU_name %in% sc_cus_no_CU_infilled, ] # get sf object of just south coast chum CUs 
scc <- chum_cu[chum_cu$CU_name %in% sc_cus, ] # get sf object of just south coast chum CUs 

# get bc_maps layers
bc <- bc_bound_hres()
nb <- bc_neighbours()

# get bounds 
bounds <- as.numeric(st_bbox(scc))

# Get palette
pal <- pnw_palette("Cascades", length(scc$CU_name), type="continuous")

# function to get fishery management areas with multiple polygons
drop_multi_poly <- function(y) {
  a <- as.data.frame(y[ y$MNGMNTR %in% names(which(table(y$MNGMNTR) >1)), ]) # get multi-polygon area rows
  small_areas <- aggregate(formula = AREA_SQM ~ MNGMNTR, data=a, FUN="min") # get smallest polygons in sets
  y[!y$AREA_SQM %in% small_areas$AREA_SQM, ] # select not small polygons
}

# plot with ggplot
png("Figures/fig_chum_CU_map.png", width=8, height=7, units="in", res=300)
ggplot(scc) +
  geom_sf(data=fma,colour="coral", size=1, fill=NA) +
  geom_sf(data=bc, fill="antiquewhite", size=0) +
  geom_sf(data=nb[nb$name=="Washington",], fill="antiquewhite", size=0) +
  geom_sf(data=scc, aes(fill=CU_name), size=0) + 
  geom_sf(data=lakes, fill="cornflowerblue", colour="cornflowerblue", size=0.001) +
  geom_sf(data=riv, fill="cornflowerblue", colour="cornflowerblue", size=0.001) +
  geom_sf_label(data=scc, aes(label = CU_name, colour=CU_name),  size=3, fontface="bold") +
  geom_sf_label(data=drop_multi_poly(fma), aes(label = MNGMNTR), size=3,label.size=0.1, colour="coral",alpha=0.7, fontface="bold") +
  annotate( geom="text", label = "BRITISH COLUMBIA", x = -121.8, y = 49.5, color = "grey22", size = 4) +
  annotate( geom="text", label = "WASHINGTON", x = -121.8, y = 48.5, color = "grey22", size = 4) +
  annotate( geom="text", label = "Pacific Ocean", x = -127.7, y = 49.3, fontface = "italic", color = "darkblue", size = 4) +
  annotate( geom="text", label = "Salish Sea", x = -123.8, y = 49.3, fontface = "italic", color = "darkblue", size = 3, angle=320) +
  scale_fill_manual( values=pal, guide=NULL) +
  scale_colour_manual( values=pal, guide=NULL) +
  coord_sf(xlim = bounds[c(1,3)] + c(0,1) , ylim = bounds[c(2,4)]) +
  scale_x_continuous(breaks=seq(-128,-122,1)) +
  xlab(NULL) +
  ylab(NULL) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering, width=unit(0.3, "in")) +
  theme_classic() +
  theme(panel.background = element_rect(fill="aliceblue") )
dev.off()

ggplot(fma) +
  geom_sf(data=fma) +
  geom_sf_label(aes(label=MNGMNTR))

# plot with leaflet
# cent <- as.data.frame(st_coordinates(st_centroid(scc)))
# cent$X[1] <- -122 # adjust Howe Sound Burrard Inlet label position
# m <-   leaflet(scc) %>%
#   #addTiles() %>%
#   addWMSTiles(baseUrl = "https://services.arcgisonline.com/arcgis/rest/services/Ocean/World_Ocean_Base/MapServer/tile/{z}/{y}/{x}.png",
#               layers = "1", options = WMSTileOptions(format = "image/png", transparent = TRUE)) %>%
#   addPolygons(color= pal, fillColor=pal) %>%
#   addPolygons(color="black", fillColor=NULL) %>%
#   addLabelOnlyMarkers(lat=cent$Y, lng=cent$X, label=scc$CU_name,
#                       labelOptions = labelOptions(textOnly = TRUE, permanent=TRUE, direction="centre" ,
#                                                    style= list(
#                                                      "color" = pal[1:7],
#                                                   #   "font-family" = "serif",
#                                                      "font-style" = "bold",
#                                                   #   "box-shadow" = "3px 3px rgba(0,0,0,0.25)",
#                                                      "font-size" = "14px"
#                                                   #   "border-color" = "rgba(0,0,0,0.5)"
#                                                   ) ))
# m 
# mapview::mapshot(m, file="Figures/fig_chum_CU_map.png")
# 

