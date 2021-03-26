# Make a map of the 7 Inside South Coast Chum Conservation Units

library(bcmaps)
library(ggplot2)
library(sf)
library(stars)
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

# plot
png("Figures/fig_chum_CU_map.png", width=7, height=7, units="in", res=300)
ggplot(data = scc) + 
  geom_sf(data=bc, fill="olivedrab") +
  geom_sf(data=nb[nb$name=="Washington",], fill="olivedrab") +
  geom_sf(data=scc, aes(fill=CU_name)) +
  geom_sf_label(data=scc, aes(label = CU_name, colour=CU_name), size=5) +
  #geom_sf(data=riv5, colour="dodgerblue") +
  scale_fill_discrete(guide=NULL) +
  scale_colour_discrete(guide=NULL) +
  coord_sf(xlim = bounds[c(1,3)], ylim = bounds[c(2,4)]) +
  theme_void() +
  theme(plot.background = element_rect(fill="paleturquoise1"))
dev.off()
