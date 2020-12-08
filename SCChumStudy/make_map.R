# Make a map of the 7 Inside South Coast Chum Conservation Units

library(bcmaps)
library(bcmapsdata)
library(ggplot2)
library(sf)
#install.packages('bcmapsdata', repos='https://bcgov.github.io/drat/', lib="C:/R/R-4.0.3/library" ) # only have to run once

# Download chum CU data shapefile from:
# https://open.canada.ca/data/en/dataset/f86c0867-d38d-4072-bd08-57cbbcbafa46

bc <- bc_bound() # get BC province boundary
rivers <- watercourses_5M # get rivers layer
chum_cu <- st_read("spatial/Chum Salmon CU Boundary_En.shp")
# get south coast chum CUs 
sc_cus <- c("Southern Coastal Streams", "Northeast Vancouver Island", "Upper Knight", "Loughborough", "Bute Inlet", "Georgia Strait", "Howe Sound-Burrard Inlet")
scc <- chum_cu[chum_cu$CU_name %in% sc_cus, ] # get sf object of just south coast chum CUs 

# get bounds 
bounds <- as.numeric(st_bbox(scc))
# plot
png("Figures/fig_CU_map.png", width=6, height=4, units="in", res=300)
ggplot(data = scc) + 
  geom_sf(data=bc) +
  geom_sf(data=scc, aes(fill=CU_name)) +
  #geom_sf(data=rivers) +
  coord_sf(xlim = bounds[c(1,3)], ylim = bounds[c(2,4)]) +
  scale_fill_discrete(name=NULL) +
  theme_void()
dev.off()
