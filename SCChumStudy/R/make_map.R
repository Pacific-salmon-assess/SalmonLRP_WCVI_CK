# Make a map of the 7 Inside South Coast Chum Conservation Units

library(bcmaps)
library(ggplot2)
library(sf)
library(stars)
library(ggspatial)
#devtools::install_github('Chrisjb/basemapR')
#library(basemapR)

#install.packages('bcmapsdata', repos='https://bcgov.github.io/drat/', lib="C:/R/R-4.0.3/library" ) # only have to run once

# Download chum CU data shapefile from:
# https://open.canada.ca/data/en/dataset/f86c0867-d38d-4072-bd08-57cbbcbafa46
chum_cu <- st_read("DataIn/chum_CU_boundary_shapefile/Chum Salmon CU Boundary_En.shp")
# get south coast chum CUs 
sc_cus <- c("Southern Coastal Streams", "Northeast Vancouver Island", "Upper Knight", "Loughborough", "Bute Inlet", "Georgia Strait", "Howe Sound-Burrard Inlet")
# Get just the 5 CUs without CU-level infilling
sc_cus_no_CU_infilled <- c("Southern Coastal Streams", "Northeast Vancouver Island", "Loughborough","Georgia Strait", "Howe Sound-Burrard Inlet")
#scc <- chum_cu[chum_cu$CU_name %in% sc_cus_no_CU_infilled, ] # get sf object of just south coast chum CUs 
scc <- chum_cu[chum_cu$CU_name %in% sc_cus, ] # get sf object of just south coast chum CUs 

# get bc_maps layers
bc <- bc_bound_hres()
riv5 <- watercourses_5M()
nb <- bc_neighbours()
#elv <- cded_stars(scc) # optional digital elevation data - takes a while

# get bounds 
bounds <- as.numeric(st_bbox(scc))

# Get base map
#rosm::osm.types() # check available open street map layers

# plot
png("Figures/fig_chum_CU_map.png", width=8, height=7, units="in", res=300)
ggplot(data = scc) + 
  geom_sf(data=bc, fill="antiquewhite") +
  geom_sf(data=nb[nb$name=="Washington",], fill="antiquewhite") +
  geom_sf(data=scc, aes(fill=CU_name)) + 
  #annotation_map_tile(zoom=5, type="hillshade") + # this takes a while
  #annotation_map_tile(zoom=10, type="osm", labels=) + # this takes a while
  #geom_sf(data=riv5, colour="dodgerblue") + # poor resolution river layer
  #base_map(bbox=st_bbox(scc), increase_zoom = 0, basemap="hydda") + # unsuccessful base map. Proxy issue?
  geom_sf_label(data=scc, aes(label = CU_name, colour=CU_name), size=3) +
  annotate( geom="text", label = "British Columbia", x = -121.8, y = 49.5, fontface = "italic", color = "grey22", size = 4) +
  annotate( geom="text", label = "Washington", x = -121.8, y = 48.5, fontface = "italic", color = "grey22", size = 4) +
  annotate( geom="text", label = "Pacific Ocean", x = -127.7, y = 49.3, fontface = "italic", color = "darkblue", size = 4) +
  annotate( geom="text", label = "Salish Sea", x = -123.5, y = 49.2, fontface = "italic", color = "darkblue", size = 3, angle=320) +
  scale_fill_discrete(guide=NULL) +
  scale_colour_discrete(guide=NULL) +
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
