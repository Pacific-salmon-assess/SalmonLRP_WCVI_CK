# Show Time-series of Correlations among spawners

#read.csv for Inlet_Sum
setwd(wcviCKDir)
# Inlet_Sum <- read.csv("DataIn/Inlet_Sum.csv")

# Inlet_Sum.df <- as.data.frame(Inlet_Sum) %>%
#   add_column( Years = 1953:2020)
#
# # Plot time-series, 1953-2020
# Inlet_Sum.df.long <- Inlet_Sum.df  %>% pivot_longer(names(Inlet_Sum.df)[1:5], names_to="Inlet_Name", values_to="Spawners")

Inlet_Sum.df.long <- read.csv("DataIn/Inlet_Sum.csv")
ggplot(Inlet_Sum.df.long, aes(BroodYear, Spawners, group=Inlet_Name, colour=Inlet_Name))+geom_line()

# from 1996 onwards, but includes years when some CUs are missing (which are not included in corr)
#Inlet_Sum.df.long <- Inlet_Sum.df  %>% filter(BroodYear>1996) %>% pivot_longer(names(Inlet_Sum.df)[1:5], names_to="Inlet_Name", values_to="Spawners")

SpawnerInlets <- ggplot(data =  Inlet_Sum.df.long %>% filter(BroodYear>1996),
       aes(BroodYear, Spawners, group=Inlet_Name, colour=Inlet_Name))+geom_line()
ggsave("Figures/SpawnerInlets.png", plot = SpawnerInlets, width = 6,
       height = 4, units = "in")
#with na.omit when some CUs are missingn (But adds lines over ommited years which is deceiving)
#Inlet_Sum.df.long <- Inlet_Sum.df  %>% na.omit() %>% pivot_longer(names(Inlet_Sum.df)[1:5], names_to="Inlet_Name", values_to="Spawners")

ggplot(data = Inlet_Sum.df.long %>% select(c(BroodYear, Inlet_Name, Spawners)) %>% filter(BroodYear>1996) %>% na.omit(), aes(BroodYear, Spawners, group=Inlet_Name, colour=Inlet_Name))+geom_line()


meanS_df <- Inlet_Sum.df.long %>% group_by(Inlet_Name) %>% summarize(meanS=mean(Spawners, na.rm=T), sdS=sd(Spawners, na.rm=T))
Inlet_Sum.df.long2 <- Inlet_Sum.df.long %>% left_join(meanS_df)
Inlet_Sum.df.long2 <- Inlet_Sum.df.long2 %>% mutate(Standardized_Spawners=(Spawners-meanS)/sdS)

ggplot(Inlet_Sum.df.long2 %>% filter(BroodYear>1996), aes(BroodYear, Standardized_Spawners, group=Inlet_Name, colour=Inlet_Name))+geom_line()


Years <- 1995:2000
tri <- matrix(NA, nrow=10,ncol=length(Years))
med <- NA



for (i in 1:length(Years)){
  Year.Start <- Years[i]
  Inlet_Sum.df <-   Inlet_Sum.df.long %>%
    filter(BroodYear %in% c(Year.Start : (Year.Start + 20) ) ) %>%
    select(-c(CU_Name, Inlet_ID, Recruits)) %>%
    # pivot_wider(id_cols=c(Inlet_Name, BroodYear), names_from=Inlet_Name,
    #             values_from=Spawners) %>%
    pivot_wider(id_cols=c(BroodYear), names_from=Inlet_Name,
                values_from=Spawners) %>%
    select(-BroodYear) %>% na.omit()

  co <- cor(Inlet_Sum.df)
  tri[,i] <- co[lower.tri(co)==TRUE]
  med[i] <- median(co[lower.tri(co)==TRUE])
}


dum <- as.data.frame(tri)
names(dum) <- as.character(Years)
tri_longer <- dum %>% pivot_longer(col=names(dum), names_to="StartYear", values_to="Correlation")
tri_longer$StartYear <- factor(tri_longer$StartYear, levels = names(dum),ordered = TRUE)

# Add overall correlation (over all years)
Inlet_Sum.df.all <-   Inlet_Sum.df.long %>%
  select(-c(CU_Name, Inlet_ID, Recruits)) %>%
  # pivot_wider(id_cols=c(Inlet_Name, BroodYear), names_from=Inlet_Name,
  #             values_from=Spawners) %>% select(-BroodYear) %>% na.omit()
  pivot_wider(id_cols=c(BroodYear), names_from=Inlet_Name,
              values_from=Spawners) %>% select(-BroodYear) %>% na.omit()
co <- cor(Inlet_Sum.df.all)
tri.all <- data.frame(StartYear="All", Correlation= co[lower.tri(co)==TRUE])
tri_longer <- tri_longer %>% add_row(tri.all)

# Plot of running 20-year correlations among inlets, over the time series 1953-2020

RunningCor <- ggplot(tri_longer, aes(x=StartYear, y=Correlation)) +
  geom_boxplot() +
  # ylab("Corrélation") +
  # xlab("Année Début") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1))
# ggsave("Figures/RunningCorrelations.png", plot = RunningCor,
#        width = 6, height = 4, units = "in")
# ggsave("Figures/RunningCorrelationsFR.png", plot = RunningCor,
#        width = 6, height = 4, units = "in")
