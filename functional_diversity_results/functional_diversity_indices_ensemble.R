# -------------------------------------------------------------------------------------------------
#Script to calculate model ensemble for functional diversity indices, Jonas Wydler, 17.01.2021
#
# -------------------------------------------------------------------------------------------------
Sys.setenv(LANGUAGE= 'en')
Sys.setenv(LANGUAGE= 'en')

library("tidyverse")
library("viridis")

#global var
months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
world2 <- map_data("world2") # world coastline for maps 
world1 <- map_data("world")


wd_fd_maps <- ("/data/") #directory where you want your maps to go
wd_ensmb <- ("/data/") #directory where you want to save the ensemble data

("/data/functional_diversity_results/functional_diversitiy_indices/") 

wd_fd_gam <- ("/data/functional_diversity_results/functional_diversitiy_indices/GAM") 
wd_fd_glm <- ("/data/functional_diversity_results/functional_diversitiy_indices/GLM") 
wd_fd_ann <- ("/data/functional_diversity_results/functional_diversitiy_indices/ANN") 

wd_clim <- ("/data/global_monthly_clims_1d/")
# # -------------------------------------------------------------------------------------------------

setwd(wd_clim)
#env variables, I use this just for the grid

apr <- read.table("glob_stack_month_apr_18_09_20.txt", h = T, sep = ";")
jul <- read.table("glob_stack_month_jul_18_09_20.txt", h = T, sep = ";")
oct <- read.table("glob_stack_month_oct_18_09_20.txt", h = T, sep = ";")
jan <- read.table("glob_stack_month_jan_18_09_20.txt", h = T, sep = ";")
feb <- read.table("glob_stack_month_feb_18_09_20.txt", h = T, sep = ";")
mar <- read.table("glob_stack_month_mar_18_09_20.txt", h = T, sep = ";")
may <- read.table("glob_stack_month_may_18_09_20.txt", h = T, sep = ";")
jun <- read.table("glob_stack_month_jun_18_09_20.txt", h = T, sep = ";")
aug <- read.table("glob_stack_month_aug_18_09_20.txt", h = T, sep = ";")
sep <- read.table("glob_stack_month_sep_18_09_20.txt", h = T, sep = ";")
nov <- read.table("glob_stack_month_nov_18_09_20.txt", h = T, sep = ";")
dec <- read.table("glob_stack_month_dec_18_09_20.txt", h = T, sep = ";")


apr <- apr[-which(apr$SSS < 20),] 
apr <- apr[-which(apr$Bathy > -175),] 
jul <- jul[-which(jul$SSS < 20),] 
jul <- jul[-which(jul$Bathy > -175),] 
oct <- oct[-which(oct$SSS < 20),] 
oct <- oct[-which(oct$Bathy > -175),] 
jan <- jan[-which(jan$SSS < 20),] 
jan <- jan[-which(jan$Bathy > -175),] 
feb <- feb[-which(feb$SSS < 20),] 
feb <- feb[-which(feb$Bathy > -175),] 
mar <- mar[-which(mar$SSS < 20),] 
mar <- mar[-which(mar$Bathy > -175),] 
may <- may[-which(may$SSS < 20),] 
may <- may[-which(may$Bathy > -175),] 
jun <- jun[-which(jun$SSS < 20),] 
jun <- jun[-which(jun$Bathy > -175),] 
aug <- aug[-which(aug$SSS < 20),] 
aug <- aug[-which(aug$Bathy > -175),] 
sep <- sep[-which(sep$SSS < 20),] 
sep <- sep[-which(sep$Bathy > -175),] 
nov <- nov[-which(nov$SSS < 20),] 
nov <- nov[-which(nov$Bathy > -175),] 
dec <- dec[-which(dec$SSS < 20),] 
dec <- dec[-which(dec$Bathy > -175),]





data_clim_annual <- rbind(jan,feb,mar,apr,may,jun,jul,aug,sep,oct,nov,dec)

data_clim_annual$id <- as.factor(data_clim_annual$id)

data_clim_annual_new <- data.frame(data_clim_annual %>% group_by(id) %>% 
                                     summarize(x = unique(x), y = unique(y), SST = mean(SST, na.rm = T),
                                               logNO3 = mean(logNO3, na.rm = T), logChl = mean(logChl, na.rm = T),
                                               logEKE = mean(logEKE, na.rm = T), Sistar = mean(Sistar, na.rm = T),
                                               Nstar = mean(Nstar, na.rm = T), dO2 = mean(dO2, na.rm = T),
                                               MLD = mean(MLD, na.rm = T), PAR = mean(PAR, na.rm = T),
                                               Wind = mean(Wind, na.rm = T), pCO2 = mean(pCO2, na.rm = T),
                                               SSS = mean(SSS, na.rm = T), logPO4 = mean(logPO4, na.rm = T)
                                     ))


data_clim_annual_new$x.2 <- data_clim_annual_new$x
data_clim_annual_new[data_clim_annual_new$x < 0 ,"x.2"] <- (data_clim_annual_new[data_clim_annual_new$x < 0 ,"x"]) + 360

data_clim_annual_new$cell_id <- factor(paste(data_clim_annual_new$x.2, data_clim_annual_new$y, sep = "_"))
data_clim_annual_new <- data_clim_annual_new[c("cell_id", "x.2", "y", "SST", "logNO3", "logChl", "logEKE", "Sistar", "Nstar",
                                               "dO2",  "MLD", "PAR", "Wind", "pCO2", "SSS", "logPO4")]

colnames(data_clim_annual_new)[2] <- "x"
# # -------------------------------------------------------------------------------------------------
setwd(wd_fd_gam)
load("df_annual_gam.Rda")
df_gam <- df_annual_gam
rm(df_annual_gam)


df_gam <- data.frame(df_gam %>% group_by(cell_id) %>%
                                   summarize(x = unique(x), y = unique(y), FEve_gam = mean(FEve, na.rm = T),
                                             FDis_gam = mean(FDis, na.rm = T), RaoQ_gam = mean(RaoQ, na.rm = T)))


df_gam$x.2 <- df_gam$x
df_gam[df_gam$x < 0 ,"x.2"] <- (df_gam[df_gam$x < 0 ,"x"]) + 360

df_gam$cell_id <- as.character(paste(df_gam$x.2, df_gam$y, sep = "_"))

df_gam <- df_gam[c("cell_id", "x.2", "y", "FEve_gam", "FDis_gam", "RaoQ_gam")]
colnames(df_gam)[2] <- "x"
                                  
# 
# # -------------------------------------------------------------------------------------------------
# 
# 
setwd(wd_fd_glm)
load("dat_convex_annual_GLM.Rda")
df_glm <- table_df
rm(table_df)

df_glm <- data.frame(df_glm %>% group_by(cell_id) %>%
                       summarize(x = unique(x), y = unique(y), FEve_glm = mean(FEve, na.rm = T),
                                 FDis_glm = mean(FDis, na.rm = T), RaoQ_glm = mean(RaoQ, na.rm = T)  ))



df_glm$x.2 <- df_glm$x
df_glm[df_glm$x < 0 ,"x.2"] <- (df_glm[df_glm$x < 0 ,"x"]) + 360

df_glm$cell_id <- as.character(paste(df_glm$x.2, df_glm$y, sep = "_"))

df_glm <- df_glm[c("cell_id", "x.2", "y", "FEve_glm", "FDis_glm", "RaoQ_glm")]
colnames(df_glm)[2] <- "x"
# # -------------------------------------------------------------------------------------------------
# 
setwd(wd_fd_ann)
load("dat_convex_annual_ANN.Rda")
df_ann <- table_df
rm(table_df)
df_ann <- data.frame(df_ann %>% group_by(cell_id) %>%
                       summarize(x = unique(x), y = unique(y), FEve_ann = mean(FEve, na.rm = T),
                                 FDis_ann = mean(FDis, na.rm = T), RaoQ_ann = mean(RaoQ, na.rm = T) ))

  
                    
df_ann$x.2 <- df_ann$x
df_ann[df_ann$x < 0 ,"x.2"] <- (df_ann[df_ann$x < 0 ,"x"]) + 360

df_ann$cell_id <- as.character(paste(df_ann$x.2, df_ann$y, sep = "_"))

df_ann <- df_ann[c("cell_id", "x.2", "y", "FEve_ann", "FDis_ann", "RaoQ_ann")]
colnames(df_ann)[2] <- "x"
# # -------------------------------------------------------------------------------------------------

setwd(wd_ensmb)

dat_ensm_part0  <- merge(data_clim_annual_new[c("cell_id","x","y")], df_gam, all.x = T, by = c("cell_id","x","y"))
dat_ensm_part1  <- merge(dat_ensm_part0, df_ann, all.x = T, by = c("cell_id","x","y"))
dat_ensm_part2  <- merge(dat_ensm_part1, df_glm, all.x = T, by = c("cell_id","x","y"))


data_ensm <- data.frame(dat_ensm_part2 %>% group_by(cell_id) %>% 
                          summarize(x = unique(x), y = unique(y), FEve = mean(c(FEve_ann, FEve_gam, FEve_glm), na.rm = F), 
                                    FDis = mean(c(FDis_gam, FDis_glm, FDis_ann), na.rm = F), RaoQ = mean(c(RaoQ_glm, RaoQ_ann, RaoQ_gam), na.rm =F)
                                    
                          ))
df2 <- data_ensm



g_FEve <- ggplot() + geom_raster(aes(x = x, y = y, fill = FEve), data = df2) +
  # play around with palettes. You can also use scale_fill_distiller() and
  geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
  coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                           panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
  theme(legend.text=element_text(size=6), legend.title = element_text(size=8)) +
  scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
  scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
  scale_fill_viridis(na.value = "white") + labs(fill = paste0("FEve")) +
  #labs(title = paste0("Copepod functional group ", fg, ", (", fg_text,"), #spp: ",number_spp,", annual mean")) +
  theme(plot.title = element_text(size=8)) + theme(legend.text=element_text(size=5)) +
  theme(legend.title = element_text(size=8))

setwd(wd_fd_maps)
ggsave(plot = g_FEve, filename = paste0("FEve_famd_ward_narmF.jpg"),dpi = 300, width = 7, height = 4)



g_FDis <- ggplot() + geom_raster(aes(x = x, y = y, fill = FDis), data = df2) +
  # play around with palettes. You can also use scale_fill_distiller() and
  geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
  coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                           panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
  theme(legend.text=element_text(size=6), legend.title = element_text(size=8)) +
  scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
  scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
  scale_fill_viridis(na.value = "white") + labs(fill = paste0("FDis")) +
  #labs(title = paste0("Copepod functional group ", fg, ", (", fg_text,"), #spp: ",number_spp,", annual mean")) +
  theme(plot.title = element_text(size=8)) + theme(legend.text=element_text(size=5)) +
  theme(legend.title = element_text(size=8))

setwd(wd_fd_maps)
ggsave(plot = g_FDis, filename = paste0("FDis_famd_ward_narmF.jpg"),dpi = 300, width = 7, height = 4)



g_RaoQ <- ggplot() + geom_raster(aes(x = x, y = y, fill = RaoQ), data = df2) +
  # play around with palettes. You can also use scale_fill_distiller() and
  geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
  coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                           panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
  theme(legend.text=element_text(size=6), legend.title = element_text(size=8)) +
  scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
  scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
  scale_fill_viridis(na.value = "white") + labs(fill = paste0("RaoQ")) +
  #labs(title = paste0("Copepod functional group ", fg, ", (", fg_text,"), #spp: ",number_spp,", annual mean")) +
  theme(plot.title = element_text(size=8)) + theme(legend.text=element_text(size=6)) +
  theme(legend.title = element_text(size=8))

setwd(wd_fd_maps)
ggsave(plot = g_RaoQ, filename = paste0("RaoQ_famd_ward_narmF.jpg"),dpi = 300, width = 7, height = 4)


panel_1 <- ggarrange(g_FEve, g_FDis, g_RaoQ,
                     labels = c( "FEve", "FDis", "RaoQ"),
                     font.label =  list(size = 12, face = "bold", color ="black"),
                     ncol = 2, nrow = 2)

panel_1 <- annotate_figure(
  panel_1,
  top = paste0("Functional diversity indices, ensemble, annual"),
  fig.lab.size  = 10
)

setwd(wd_fd_maps)
ggsave(plot = panel_1, filename = paste0("Functional_diversity_indices_ensemble_annual_narmF.jpg"),dpi = 300, width = 16, height = 9)

setwd(wd_ensmb)
data_fd_indices_ensemble <- df2
save(data_fd_indices_ensemble, file = "data_fd_indices_ensemble_narmF.RData")






