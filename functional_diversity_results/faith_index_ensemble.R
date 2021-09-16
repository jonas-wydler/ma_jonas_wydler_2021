# -------------------------------------------------------------------------------------------------
#map fd, 17.01.2021
#
#Update 01.02.2021: The current plan is to use ses.pd (faith index) for PDichness and the the indices from the FD packages for everything else.
#Here we compute the ensemble to use it for correlations with env. properties. PDic and FDic need the convex hulls to work which sadly doesn't seem
#to.
#updated 17.02.2021
#updated 20.3 2nd go
# -------------------------------------------------------------------------------------------------
Sys.setenv(LANGUAGE= 'en')

library("tidyverse")
library("reshape2")
library("biomod2")
library("viridis")
#library("vegan")
library("FactoMineR")
library("dplyr")
library("parallel")
library("FD")
library("missMDA")
library("ggpubr")

#global var
months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
world2 <- map_data("world2") # world coastline for maps 
world1 <- map_data("world")

wd_fd_maps <- ("/data/") #directory where you want your maps to go
wd_ensmb <- ("/data/") #directory where you want to save the ensemble data

wd_fd_dat_nonses_faith_gam <- ("/data/functional_diversity_results/functional_richness_faith_non_ses/faith_nonses_bin_famd_ward_GAM/")
wd_fd_dat_nonses_faith_glm <- ("/data/functional_diversity_results/functional_richness_faith_non_ses/faith_nonses_bin_famd_ward_GLM/")
wd_fd_dat_nonses_faith_ann <- ("/data/functional_diversity_results/functional_richness_faith_non_ses/faith_nonses_bin_famd_ward_ANN/")

wd_fd_dat_ses_faith_gam <- ("/data/functional_diversity_results/functional_richness_faith_ses/faith_ses_bin_famd_ward_GAM/")
wd_fd_dat_ses_faith_glm <- ("/data/functional_diversity_results/functional_richness_faith_ses/faith_ses_bin_famd_ward_GLM/")
wd_fd_dat_ses_faith_ann <- ("/data/functional_diversity_results/functional_richness_faith_ses/faith_ses_bin_famd_ward_ANN)
wd_clim <- ("/data/global_monthly_clims_1d/")


# # -------------------------------------------------------------------------------------------------
setwd(wd_clim)
#env variables

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


# setwd(wd_fd_gam)
# load("df_annual_gam.Rda")
# df_gam <- df_annual_gam
# rm(df_annual_gam)
# 
# 
# df_gam <- data.frame(df_gam %>% group_by(cell_id) %>% 
#                                    summarize(x = unique(x), y = unique(y), FEve = mean(FEve, na.rm = T), 
#                                              FDis = mean(FDis, na.rm = T), RaoQ = mean(RaoQ, na.rm = T)
#                                              
#                                    ))
# 
# # -------------------------------------------------------------------------------------------------
# 
# 
# setwd(wd_fd_glm)
# load("dat_convex_annual_GLM.Rda")
# df_glm <- table_df
# rm(table_df)
# 
# df_glm <- data.frame(df_glm %>% group_by(cell_id) %>% 
#                        summarize(x = unique(x), y = unique(y), FEve = mean(FEve, na.rm = T), 
#                                  FDis = mean(FDis, na.rm = T), RaoQ = mean(RaoQ, na.rm = T)
#                                  
#                        ))
# 
# # -------------------------------------------------------------------------------------------------
# 
# setwd(wd_fd_ann)
# load("dat_convex_annual_ANN.Rda")
# df_ann <- table_df
# rm(table_df)
# df_ann <- data.frame(df_ann %>% group_by(cell_id) %>% 
#                        summarize(x = unique(x), y = unique(y), FEve = mean(FEve, na.rm = T), 
#                                  FDis = mean(FDis, na.rm = T), RaoQ = mean(RaoQ, na.rm = T)
#                                  
#                        ))
# 
# # -------------------------------------------------------------------------------------------------
# gc()
# 
# data_ensm <- rbind(df_ann,df_gam,df_glm)
# data_ensm <- data.frame(data_ensm %>% group_by(cell_id) %>% 
#                        summarize(x = unique(x), y = unique(y), FEve = mean(FEve, na.rm = T), 
#                                  FDis = mean(FDis, na.rm = T), RaoQ = mean(RaoQ, na.rm = T)
#                                  
#                        ))
# 
# 
# 
# df2 <- data_ensm
# df2$x.2 <- df2$x
# df2[df2$x < 0 ,"x.2"] <- (df2[df2$x < 0 ,"x"]) + 360
# 
# 
# 
# g_FEve <- ggplot() + geom_raster(aes(x = x.2, y = y, fill = FEve), data = df2) +
#   # play around with palettes. You can also use scale_fill_distiller() and
#   geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
#   coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#                            panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
#   theme(legend.text=element_text(size=6), legend.title = element_text(size=8)) +
#   scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
#   scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
#   scale_fill_viridis(na.value = "white") + labs(fill = paste0("FEve")) +
#   #labs(title = paste0("Copepod functional group ", fg, ", (", fg_text,"), #spp: ",number_spp,", annual mean")) +
#   theme(plot.title = element_text(size=8)) + theme(legend.text=element_text(size=5)) +
#   theme(legend.title = element_text(size=8))
# 
# 
# g_FDis <- ggplot() + geom_raster(aes(x = x.2, y = y, fill = FDis), data = df2) +
#   # play around with palettes. You can also use scale_fill_distiller() and
#   geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
#   coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#                            panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
#   theme(legend.text=element_text(size=6), legend.title = element_text(size=8)) +
#   scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
#   scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
#   scale_fill_viridis(na.value = "white") + labs(fill = paste0("FDis")) +
#   #labs(title = paste0("Copepod functional group ", fg, ", (", fg_text,"), #spp: ",number_spp,", annual mean")) +
#   theme(plot.title = element_text(size=8)) + theme(legend.text=element_text(size=5)) +
#   theme(legend.title = element_text(size=8))
# 
# g_RaoQ <- ggplot() + geom_raster(aes(x = x.2, y = y, fill = RaoQ), data = df2) +
#   # play around with palettes. You can also use scale_fill_distiller() and
#   geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
#   coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#                            panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
#   theme(legend.text=element_text(size=6), legend.title = element_text(size=8)) +
#   scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
#   scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
#   scale_fill_viridis(na.value = "white") + labs(fill = paste0("RaoQ")) +
#   #labs(title = paste0("Copepod functional group ", fg, ", (", fg_text,"), #spp: ",number_spp,", annual mean")) +
#   theme(plot.title = element_text(size=8)) + theme(legend.text=element_text(size=5)) +
#   theme(legend.title = element_text(size=8))
# 
# panel_1 <- ggarrange(g_FEve, g_FDis, g_RaoQ,
#                      labels = c( "FEve", "FDis", "RaoQ"),
#                      font.label =  list(size = 12, face = "bold", color ="black"),
#                      ncol = 2, nrow = 2)
# 
# panel_1 <- annotate_figure(
#   panel_1,
#   top = paste0("Functional diversity indices, ensemble, annual"),
#   fig.lab.size  = 10
# )
# 
# setwd(wd_fd_maps)
# ggsave(plot = panel_1, filename = paste0("Functional_diversity_indices_ensemble_annual.jpg"),dpi = 600, width = 16, height = 9)
# 
# setwd(wd_ensmb)
# data_fd_indices_ensemble <- df2
# save(data_fd_indices_ensemble, file = "data_fd_indices_ensemble.RData")
# 

# -------------------------------------------------------------------------------------------------
#famd ward without ses
setwd(wd_fd_dat_nonses_faith_gam)

# load("dat_binary_GAM_apr_famd_ward_v2.Rda")
# dat_gam_apr <- df; rm(df);
# 
# load("dat_binary_GAM_aug_famd_ward_v2.Rda")
# dat_gam_aug <- df; rm(df);
# 
# load("dat_binary_GAM_dec_famd_ward_v2.Rda")
# dat_gam_dec <- df; rm(df);
# 
# load("dat_binary_GAM_feb_famd_ward_v2.Rda")
# dat_gam_feb <- df; rm(df);
# 
# load("dat_binary_GAM_jan_famd_ward_v2.Rda")
# dat_gam_jan <- df; rm(df);
# 
# load("dat_binary_GAM_jul_famd_ward_v2.Rda")
# dat_gam_jul <- df; rm(df);
# 
# load("dat_binary_GAM_jun_famd_ward_v2.Rda")
# dat_gam_jun <- df; rm(df);
# 
# load("dat_binary_GAM_mar_famd_ward_v2.Rda")
# dat_gam_mar <- df; rm(df);
# 
# load("dat_binary_GAM_may_famd_ward_v2.Rda")
# dat_gam_may <- df; rm(df);
# 
# load("dat_binary_GAM_nov_famd_ward_v2.Rda")
# dat_gam_nov <- df; rm(df);
# 
# load("dat_binary_GAM_oct_famd_ward_v2.Rda")
# dat_gam_oct <- df; rm(df);
# 
# load("dat_binary_GAM_sep_famd_ward_v2.Rda")
# dat_gam_sept <- df; rm(df);

setwd(wd_fd_dat_nonses_faith_gam)


load("dat_faith_annual_GAM_binary_famd_ward_v2.Rda")
dat_gam <- table_df;rm(table_df)
dat_gam$cell_id <- as.character(paste(dat_gam$x.2, dat_gam$y, sep = "_"))
 

dat_gam <- data.frame(dat_gam %>% group_by(cell_id) %>% 
                        summarize(x = unique(x.2), y = unique(y), SR_gam = mean(SR, na.rm = T), FR_gam = mean(PD, na.rm = T)))




setwd(wd_fd_dat_nonses_faith_glm)


# load("dat_binary_GLM_apr_famd_ward_v2.Rda")
# dat_glm_apr <- df; rm(df);
# 
# load("dat_binary_GLM_aug_famd_ward_v2.Rda")
# dat_glm_aug <- df; rm(df);
# 
# load("dat_binary_GLM_dec_famd_ward_v2.Rda")
# dat_glm_dec <- df; rm(df);
# 
# load("dat_binary_GLM_feb_famd_ward_v2.Rda")
# dat_glm_feb <- df; rm(df);
# 
# load("dat_binary_GLM_jan_famd_ward_v2.Rda")
# dat_glm_jan <- df; rm(df);
# 
# load("dat_binary_GLM_jul_famd_ward_v2.Rda")
# dat_glm_jul <- df; rm(df);
# 
# load("dat_binary_GLM_jun_famd_ward_v2.Rda")
# dat_glm_jun <- df; rm(df);
# 
# load("dat_binary_GLM_mar_famd_ward_v2.Rda")
# dat_glm_mar <- df; rm(df);
# 
# load("dat_binary_GLM_may_famd_ward_v2.Rda")
# dat_glm_may <- df; rm(df);
# 
# load("dat_binary_GLM_nov_famd_ward_v2.Rda")
# dat_glm_nov <- df; rm(df);
# 
# load("dat_binary_GLM_oct_famd_ward_v2.Rda")
# dat_glm_oct <- df; rm(df);
# 
# load("dat_binary_GLM_sep_famd_ward_v2.Rda")
# dat_glm_sept <- df; rm(df);


load("dat_faith_annual_GLM_binary_famd_ward_v2.Rda")
dat_glm <- table_df;rm(table_df)
dat_glm$cell_id <- as.character(paste(dat_glm$x.2, dat_glm$y, sep = "_"))

dat_glm <- data.frame(dat_glm %>% group_by(cell_id) %>% 
                        summarize(x = unique(x.2), y = unique(y),SR_glm = mean(SR, na.rm = T), FR_glm = mean(PD, na.rm = T)))



setwd(wd_fd_dat_nonses_faith_ann)


# load("dat_binary_ANN_apr_famd_ward_v2.Rda")
# dat_ann_apr <- df; rm(df);
# 
# load("dat_binary_ANN_aug_famd_ward_v2.Rda")
# dat_ann_aug <- df; rm(df);
# 
# load("dat_binary_ANN_dec_famd_ward_v2.Rda")
# dat_ann_dec <- df; rm(df);
# 
# load("dat_binary_ANN_feb_famd_ward_v2.Rda")
# dat_ann_feb <- df; rm(df);
# 
# load("dat_binary_ANN_jan_famd_ward_v2.Rda")
# dat_ann_jan <- df; rm(df);
# 
# load("dat_binary_ANN_jul_famd_ward_v2.Rda")
# dat_ann_jul <- df; rm(df);
# 
# load("dat_binary_ANN_jun_famd_ward_v2.Rda")
# dat_ann_jun <- df; rm(df);
# 
# load("dat_binary_ANN_mar_famd_ward_v2.Rda")
# dat_ann_mar <- df; rm(df);
# 
# load("dat_binary_ANN_may_famd_ward_v2.Rda")
# dat_ann_may <- df; rm(df);
# 
# load("dat_binary_ANN_nov_famd_ward_v2.Rda")
# dat_ann_nov <- df; rm(df);
# 
# load("dat_binary_ANN_oct_famd_ward_v2.Rda")
# dat_ann_oct <- df; rm(df);
# 
# load("dat_binary_ANN_sep_famd_ward_v2.Rda")
# dat_ann_sept <- df; rm(df);
# 



load("dat_faith_annual_ANN_binary_famd_ward_v2.Rda")
dat_ann <- table_df;rm(table_df)
dat_ann$cell_id <- as.character(paste(dat_ann$x.2, dat_ann$y, sep = "_"))


dat_ann <- data.frame(dat_ann %>% group_by(cell_id) %>% 
                        summarize(x = unique(x.2), y = unique(y),SR_ann = mean(SR, na.rm = T), FR_ann = mean(PD, na.rm = T)))



#dat_esm <- bind_rows(dat_gam, dat_glm, dat_ann)
#dat_esm_monthly <- rbind(dat_ann_jan, dat_ann_feb, dat_ann_mar, dat_ann_apr, dat_ann_may, dat_ann_jun, dat_ann_jul, 
                             #dat_ann_aug, dat_ann_sept, dat_ann_oct, dat_ann_nov, dat_ann_dec)
                             
                             # dat_gam_jan, dat_gam_feb, dat_gam_mar, dat_gam_apr, dat_gam_may, dat_gam_jun, dat_gam_jul,
                             # dat_gam_aug, dat_gam_sept, dat_gam_oct, dat_gam_nov, dat_gam_dec,
                             # dat_glm_jan, dat_glm_feb, dat_glm_mar, dat_glm_apr, dat_glm_may, dat_glm_jun, dat_glm_jul,
                             # dat_glm_aug, dat_glm_sept, dat_glm_oct, dat_glm_nov, dat_glm_dec,
                             # )

#dat_gam$PD_GAM <- dat_gam$PD
#dat_ann$PD_ANN <- dat_ann$PD
#dat_glm$PD_GLM <- dat_glm$PD

setwd(wd_ensmb)

dat_ensm_part0  <- merge(data_clim_annual_new[c("cell_id","x","y")], dat_gam, all.x = T, by = c("cell_id","x","y"))
dat_ensm_part1  <- merge(dat_ensm_part0, dat_ann, all.x = T, by = c("cell_id","x","y"))
dat_ensm_part2  <- merge(dat_ensm_part1, dat_glm, all.x = T, by = c("cell_id","x","y"))



dat_esm2 <- data.frame(dat_ensm_part2 %>% group_by(cell_id) %>% 
                          summarize(x = unique(x), y = unique(y), SR =  mean(c(SR_gam,SR_glm,SR_ann), na.rm = F), 
                                    FR =  mean(c(FR_gam, FR_glm, FR_ann), na.rm = F)))



df2 <- dat_esm2
df2$x.2 <- df2$x
df2[df2$x < 0 ,"x.2"] <- (df2[df2$x < 0 ,"x"]) + 360


g1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FR), data = dat_esm2) +
  # play around with palettes. You can also use scale_fill_distiller() and
  geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
  coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                           panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
  theme(legend.text=element_text(size=8)) +
  theme(legend.title = element_text(size=12)) + 
  scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
  scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
  scale_fill_viridis(na.value = "white") + labs(fill = "FR") 


setwd(wd_fd_maps)
data_fd_fait_nonses_ward_ensemble <- dat_esm2
save(data_fd_fait_nonses_ward_ensemble, file = "data_fd_fait_nonses_ward_ensembleNaRmF.RData")
ggsave(plot = g1, filename = paste0("FR_binary_annual_ann_famd_ward_narmF.jpg"),dpi = 300, width = 7, height = 4)
# # -------------------------------------------------------------------------------------------------

#famd ward with ses
setwd(wd_fd_dat_ses_faith_gam)
load("ses_dat_faith_annual_GAM_binary_famd_ward_v2.Rda")
dat_gam <- table_df;rm(table_df)

df2 <- dat_gam
df2$x.2 <- df2$x
df2[df2$x < 0 ,"x.2"] <- (df2[df2$x < 0 ,"x"]) + 360

dat_gam <- df2[c("cell_id", "x.2", "y", "pd.obs.z")]
colnames(dat_gam)[2] <- "x"
dat_gam$cell_id <- as.character(paste(dat_gam$x, dat_gam$y, sep = "_"))
colnames(dat_gam)[2] <- "x"
colnames(dat_gam)[4] <- "PD_GAM"
# dat_gam <- data.frame(dat_gam %>% group_by(cell_id) %>% 
#                         summarize(x = unique(x), y = unique(y), PD_GAM = mean(pd.obs.z, na.rm = T)))
# 



setwd(wd_fd_dat_ses_faith_glm)
load("ses_dat_faith_annual_GLM_binary_famd_ward_v2.Rda")
dat_glm <- table_df;rm(table_df)

df2 <- dat_glm
df2$x.2 <- df2$x
df2[df2$x < 0 ,"x.2"] <- (df2[df2$x < 0 ,"x"]) + 360

dat_glm <- df2[c("cell_id", "x.2", "y", "pd.obs.z")]
dat_glm$cell_id <- as.character(paste(dat_glm$x, dat_glm$y, sep = "_"))
colnames(dat_glm)[2] <- "x"
colnames(dat_glm)[4] <- "PD_GLM"
# dat_glm <- data.frame(dat_glm %>% group_by(cell_id) %>% 
#                         summarize(x = unique(x), y = unique(y), PD_GLM = mean(pd.obs.z, na.rm = T)))







setwd(wd_fd_dat_ses_faith_ann)
load("ses_dat_faith_annual_ANN_binary_famd_ward_v2.Rda")
dat_ann <- table_df;rm(table_df)

df2 <- dat_ann
df2$x.2 <- df2$x
df2[df2$x < 0 ,"x.2"] <- (df2[df2$x < 0 ,"x"]) + 360

dat_ann <- df2[c("cell_id", "x.2", "y", "pd.obs.z")]
colnames(dat_ann)[2] <- "x"
dat_ann$cell_id <- as.character(paste(dat_ann$x, dat_ann$y, sep = "_"))
colnames(dat_ann)[4] <- "PD_ANN"
# dat_ann <- data.frame(dat_ann %>% group_by(cell_id) %>% 
#                         summarize(x = unique(x), y = unique(y), PD_ANN = mean(pd.obs.z, na.rm = T)))
# 





setwd(wd_ensmb)

dat_ensm_part0  <- merge(data_clim_annual_new[c("cell_id","x","y")], dat_gam[c("cell_id","x","y", "PD_GAM")], all.x = T, by = c("cell_id","x","y"))
dat_ensm_part1  <- merge(dat_ensm_part0, dat_ann[c("cell_id","x","y", "PD_ANN")], all.x = T, by = c("cell_id","x","y"))
dat_ensm_part2  <- merge(dat_ensm_part1, dat_glm[c("cell_id","x","y", "PD_GLM")], all.x = T, by = c("cell_id","x","y"))




dat_esm2 <- data.frame(dat_ensm_part2 %>% group_by(cell_id) %>% 
                         summarize(x = unique(x), y = unique(y), sesFR =  mean(c(PD_GAM,PD_GLM,PD_ANN), na.rm = F)))

                      
                      # dat_esm <- data.frame(dat_esm %>% group_by(cell_id) %>% 
                      #   summarize(x = unique(x), y = unique(y), sesfr = mean(pd.obs.z, na.rm = T), 
                      #             p = mean(pd.obs.p, na.rm = T)
                      #   ))



df2 <- dat_esm2
df2$x.2 <- df2$x
df2[df2$x < 0 ,"x.2"] <- (df2[df2$x < 0 ,"x"]) + 360


g3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = sesFR), data = dat_esm2) +
  # play around with palettes. You can also use scale_fill_distiller() and
  geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
  coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                           panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
  theme(legend.text=element_text(size=6), legend.title = element_text(size=8)) +
  scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
  scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
  scale_fill_gradient2(name = "ses FR",space = "Lab", mid = "white", low = "#3288BD", high = "#D53E4F", na.value = "white", aesthetics = "fill",guide = "colourbar")+
  #labs(title = paste0("ses FR, based on hsi, annual mean")) +
  theme(plot.title = element_text(size=12)) + theme(legend.text=element_text(size=8)) +
  theme(legend.title = element_text(size=12)) 





setwd(wd_fd_maps)
data_fd_fait_ses_ward_ensemble <- df2
save(data_fd_fait_ses_ward_ensemble, file = "data_fd_fait_ses_ward_ensemble_naf.RData")

ggsave(plot = g3, filename = paste0("FR_binary_annual_ensemble_famd_ward_ses_v2_naf.jpg"),dpi = 300, width = 7, height = 4)
# # -------------------------------------------------------------------------------------------------



