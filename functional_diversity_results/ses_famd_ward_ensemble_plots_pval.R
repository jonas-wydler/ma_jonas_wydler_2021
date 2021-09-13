## -------------------------------------------------------------------------------------------------
#map fd, 17.02.2021
#
#calculating the freq of signif values for ses faith index
#updated to ward on 20.3
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


wd_dat_gam <- ("/data/functional_diversity_results/ses/faith_ses_bin_famd_ward_GAM//")
wd_dat_glm <- ("/data/functional_diversity_results//ses/faith_ses_bin_famd_ward_GLM//")
wd_dat_ann <- ("/data/functional_diversity_results/ses/faith_ses_bin_famd_ward_ANN//")
wd_map <- ("/data/")
wd_clim <- wd_clim <- ("/data/global_monthly_clims_1d/")

# -------------------------------------------------------------------------------------------------
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

#GAM
setwd(wd_dat_gam)

load("ses_dat_binary_GAM_jan_famd_ward_v2.Rda")
dat_gam_jan <- df;rm(df);
dat_gam_jan$pd.obs.p <- ifelse((dat_gam_jan$pd.obs.p < 0.05 | dat_gam_jan$pd.obs.p > 0.95), 1, 0)

load("ses_dat_binary_GAM_feb_famd_ward_v2.Rda")
dat_gam_feb <- df;rm(df);
dat_gam_feb$pd.obs.p <- ifelse((dat_gam_feb$pd.obs.p < 0.05 | dat_gam_feb$pd.obs.p > 0.95), 1, 0)


load("ses_dat_binary_GAM_mar_famd_ward_v2.Rda")
dat_gam_mar <- df;rm(df);
dat_gam_mar$pd.obs.p <- ifelse((dat_gam_mar$pd.obs.p < 0.05 | dat_gam_mar$pd.obs.p > 0.95), 1, 0)


load("ses_dat_binary_GAM_apr_famd_ward_v2.Rda")
dat_gam_apr <- df;rm(df);
dat_gam_apr$pd.obs.p <- ifelse((dat_gam_apr$pd.obs.p < 0.05 | dat_gam_apr$pd.obs.p > 0.95), 1, 0)


load("ses_dat_binary_GAM_may_famd_ward_v2.Rda")
dat_gam_may <- df;rm(df);
dat_gam_may$pd.obs.p <- ifelse((dat_gam_may$pd.obs.p < 0.05 | dat_gam_may$pd.obs.p > 0.95), 1, 0)


load("ses_dat_binary_GAM_jun_famd_ward_v2.Rda")
dat_gam_jun <- df;rm(df);
dat_gam_jun$pd.obs.p <- ifelse((dat_gam_jun$pd.obs.p < 0.05 | dat_gam_jun$pd.obs.p > 0.95), 1, 0)


load("ses_dat_binary_GAM_jul_famd_ward_v2.Rda")
dat_gam_jul <- df;rm(df);
dat_gam_jul$pd.obs.p <- ifelse((dat_gam_jul$pd.obs.p < 0.05 | dat_gam_jul$pd.obs.p > 0.95), 1, 0)


load("ses_dat_binary_GAM_aug_famd_ward_v2.Rda")
dat_gam_aug <- df;rm(df);
dat_gam_aug$pd.obs.p <- ifelse((dat_gam_aug$pd.obs.p < 0.05 | dat_gam_aug$pd.obs.p > 0.95), 1, 0)

load("ses_dat_binary_GAM_sep_famd_ward_v2.Rda")
dat_gam_sep <- df;rm(df);
dat_gam_sep$pd.obs.p <- ifelse((dat_gam_sep$pd.obs.p < 0.05 | dat_gam_sep$pd.obs.p > 0.95), 1, 0)

load("ses_dat_binary_GAM_oct_famd_ward_v2.Rda")
dat_gam_oct <- df;rm(df);
dat_gam_oct$pd.obs.p <- ifelse((dat_gam_oct$pd.obs.p < 0.05 | dat_gam_oct$pd.obs.p > 0.95), 1, 0)

load("ses_dat_binary_GAM_nov_famd_ward_v2.Rda")
dat_gam_nov <- df;rm(df);
dat_gam_nov$pd.obs.p <- ifelse((dat_gam_nov$pd.obs.p < 0.05 | dat_gam_nov$pd.obs.p > 0.95), 1, 0)

load("ses_dat_binary_GAM_dec_famd_ward_v2.Rda")
dat_gam_dec <- df;rm(df);gc();
dat_gam_dec$pd.obs.p <- ifelse((dat_gam_dec$pd.obs.p < 0.05), 1, 0)


dat_gam_annual <- rbind(dat_gam_jan, dat_gam_feb, dat_gam_mar, dat_gam_apr, dat_gam_may, dat_gam_jun, dat_gam_jul, dat_gam_aug,
                       dat_gam_sep, dat_gam_oct, dat_gam_nov, dat_gam_dec)

dat_gam_annual <- data.frame(dat_gam_annual %>% group_by(cell_id) %>% 
                        summarize(x = unique(x), y = unique(y), ses = mean(pd.obs.z, na.rm = T), 
                                  p = mean(pd.obs.p, na.rm = T)
                        ))

df1 <- dat_gam_annual
df1$x.2 <- df1$x
df1[df1$x < 0 ,"x.2"] <- (df1[df1$x < 0 ,"x"]) + 360


g3 <- ggplot() + geom_raster(aes(x = x.2, y = y, fill = p), data = df1) +
  # play around with palettes. You can also use scale_fill_distiller() and
  geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
  coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                           panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
  theme(legend.text=element_text(size=8)) +
  scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
  scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
  scale_fill_viridis(na.value = "white") + labs(fill = "Freq. of signif. values:") 

setwd(wd_map)
ggsave(plot = g3, filename = "Freq_of_signif_values_ses_famd_ward_GAM.jpg", width  = 7, height = 4, dpi = 300)

# -------------------------------------------------------------------------------------------------
#ANN
setwd(wd_dat_ann)

load("ses_dat_binary_ANN_jan_famd_ward_v2.Rda")
dat_ann_jan <- df;rm(df);
dat_ann_jan$pd.obs.p <- ifelse((dat_ann_jan$pd.obs.p < 0.05 | dat_ann_jan$pd.obs.p > 0.95), 1, 0)

load("ses_dat_binary_ANN_feb_famd_ward_v2.Rda")
dat_ann_feb <- df;rm(df);
dat_ann_feb$pd.obs.p <- ifelse((dat_ann_feb$pd.obs.p < 0.05 | dat_ann_feb$pd.obs.p > 0.95), 1, 0)


load("ses_dat_binary_ANN_mar_famd_ward_v2.Rda")
dat_ann_mar <- df;rm(df);
dat_ann_mar$pd.obs.p <- ifelse((dat_ann_mar$pd.obs.p < 0.05 | dat_ann_mar$pd.obs.p > 0.95), 1, 0)


load("ses_dat_binary_ANN_apr_famd_ward_v2.Rda")
dat_ann_apr <- df;rm(df);
dat_ann_apr$pd.obs.p <- ifelse((dat_ann_apr$pd.obs.p < 0.05 | dat_ann_apr$pd.obs.p > 0.95), 1, 0)


load("ses_dat_binary_ANN_may_famd_ward_v2.Rda")
dat_ann_may <- df;rm(df);
dat_ann_may$pd.obs.p <- ifelse((dat_ann_may$pd.obs.p < 0.05 | dat_ann_may$pd.obs.p > 0.95), 1, 0)


load("ses_dat_binary_ANN_jun_famd_ward_v2.Rda")
dat_ann_jun <- df;rm(df);
dat_ann_jun$pd.obs.p <- ifelse((dat_ann_jun$pd.obs.p < 0.05 | dat_ann_jun$pd.obs.p > 0.95), 1, 0)


load("ses_dat_binary_ANN_jul_famd_ward_v2.Rda")
dat_ann_jul <- df;rm(df);
dat_ann_jul$pd.obs.p <- ifelse((dat_ann_jul$pd.obs.p < 0.05 | dat_ann_jul$pd.obs.p > 0.95), 1, 0)


load("ses_dat_binary_ANN_aug_famd_ward_v2.Rda")
dat_ann_aug <- df;rm(df);
dat_ann_aug$pd.obs.p <- ifelse((dat_ann_aug$pd.obs.p < 0.05 | dat_ann_aug$pd.obs.p > 0.95), 1, 0)

load("ses_dat_binary_ANN_sep_famd_ward_v2.Rda")
dat_ann_sep <- df;rm(df);
dat_ann_sep$pd.obs.p <- ifelse((dat_ann_sep$pd.obs.p < 0.05 | dat_ann_sep$pd.obs.p > 0.95), 1, 0)

load("ses_dat_binary_ANN_oct_famd_ward_v2.Rda")
dat_ann_oct <- df;rm(df);
dat_ann_oct$pd.obs.p <- ifelse((dat_ann_oct$pd.obs.p < 0.05 | dat_ann_oct$pd.obs.p > 0.95), 1, 0)

load("ses_dat_binary_ANN_nov_famd_ward_v2.Rda")
dat_ann_nov <- df;rm(df);
dat_ann_nov$pd.obs.p <- ifelse((dat_ann_nov$pd.obs.p < 0.05 | dat_ann_nov$pd.obs.p > 0.95), 1, 0)

load("ses_dat_binary_ANN_dec_famd_ward_v2.Rda")
dat_ann_dec <- df;rm(df);gc();
dat_ann_dec$pd.obs.p <- ifelse((dat_ann_dec$pd.obs.p < 0.05), 1, 0)


dat_ann_annual <- rbind(dat_ann_jan, dat_ann_feb, dat_ann_mar, dat_ann_apr, dat_ann_may, dat_ann_jun, dat_ann_jul, dat_ann_aug,
                       dat_ann_sep, dat_ann_oct, dat_ann_nov, dat_ann_dec)

dat_ann_annual <- data.frame(dat_ann_annual %>% group_by(cell_id) %>% 
                              summarize(x = unique(x), y = unique(y), ses = mean(pd.obs.z, na.rm = T), 
                                        p = mean(pd.obs.p, na.rm = T)
                              ))

df1 <- dat_ann_annual
df1$x.2 <- df1$x
df1[df1$x < 0 ,"x.2"] <- (df1[df1$x < 0 ,"x"]) + 360


g3 <- ggplot() + geom_raster(aes(x = x.2, y = y, fill = p), data = df1) +
  # play around with palettes. You can also use scale_fill_distiller() and
  geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
  coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                           panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
  theme(legend.text=element_text(size=8)) +
  scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
  scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
  scale_fill_viridis(na.value = "white") + labs(fill = "Freq. of signif. values:") 

setwd(wd_map)
ggsave(plot = g3, filename = "Freq_of_signif_values_ses_famd_ward_ANN.jpg", width  = 7, height = 4, dpi = 300)

# -------------------------------------------------------------------------------------------------
#GLM
setwd(wd_dat_glm)

load("ses_dat_binary_GLM_jan_famd_ward_v2.Rda")
dat_glm_jan <- df;rm(df);
dat_glm_jan$pd.obs.p <- ifelse((dat_glm_jan$pd.obs.p < 0.05 | dat_glm_jan$pd.obs.p > 0.95), 1, 0)

load("ses_dat_binary_GLM_feb_famd_ward_v2.Rda")
dat_glm_feb <- df;rm(df);
dat_glm_feb$pd.obs.p <- ifelse((dat_glm_feb$pd.obs.p < 0.05 | dat_glm_feb$pd.obs.p > 0.95), 1, 0)


load("ses_dat_binary_GLM_mar_famd_ward_v2.Rda")
dat_glm_mar <- df;rm(df);
dat_glm_mar$pd.obs.p <- ifelse((dat_glm_mar$pd.obs.p < 0.05 | dat_glm_mar$pd.obs.p > 0.95), 1, 0)


load("ses_dat_binary_GLM_apr_famd_ward_v2.Rda")
dat_glm_apr <- df;rm(df);
dat_glm_apr$pd.obs.p <- ifelse((dat_glm_apr$pd.obs.p < 0.05 | dat_glm_apr$pd.obs.p > 0.95), 1, 0)


load("ses_dat_binary_GLM_may_famd_ward_v2.Rda")
dat_glm_may <- df;rm(df);
dat_glm_may$pd.obs.p <- ifelse((dat_glm_may$pd.obs.p < 0.05 | dat_glm_may$pd.obs.p > 0.95), 1, 0)


load("ses_dat_binary_GLM_jun_famd_ward_v2.Rda")
dat_glm_jun <- df;rm(df);
dat_glm_jun$pd.obs.p <- ifelse((dat_glm_jun$pd.obs.p < 0.05 | dat_glm_jun$pd.obs.p > 0.95), 1, 0)


load("ses_dat_binary_GLM_jul_famd_ward_v2.Rda")
dat_glm_jul <- df;rm(df);
dat_glm_jul$pd.obs.p <- ifelse((dat_glm_jul$pd.obs.p < 0.05 | dat_glm_jul$pd.obs.p > 0.95), 1, 0)


load("ses_dat_binary_GLM_aug_famd_ward_v2.Rda")
dat_glm_aug <- df;rm(df);
dat_glm_aug$pd.obs.p <- ifelse((dat_glm_aug$pd.obs.p < 0.05 | dat_glm_aug$pd.obs.p > 0.95), 1, 0)

load("ses_dat_binary_GLM_sep_famd_ward_v2.Rda")
dat_glm_sep <- df;rm(df);
dat_glm_sep$pd.obs.p <- ifelse((dat_glm_sep$pd.obs.p < 0.05 | dat_glm_sep$pd.obs.p > 0.95), 1, 0)

load("ses_dat_binary_GLM_oct_famd_ward_v2.Rda")
dat_glm_oct <- df;rm(df);
dat_glm_oct$pd.obs.p <- ifelse((dat_glm_oct$pd.obs.p < 0.05 | dat_glm_oct$pd.obs.p > 0.95), 1, 0)

load("ses_dat_binary_GLM_nov_famd_ward_v2.Rda")
dat_glm_nov <- df;rm(df);
dat_glm_nov$pd.obs.p <- ifelse((dat_glm_nov$pd.obs.p < 0.05 | dat_glm_nov$pd.obs.p > 0.95), 1, 0)

load("ses_dat_binary_GLM_dec_famd_ward_v2.Rda")
dat_glm_dec <- df;rm(df);gc();
dat_glm_dec$pd.obs.p <- ifelse((dat_glm_dec$pd.obs.p < 0.05), 1, 0)


dat_glm_annual <- rbind(dat_glm_jan, dat_glm_feb, dat_glm_mar, dat_glm_apr, dat_glm_may, dat_glm_jun, dat_glm_jul, dat_glm_aug,
                       dat_glm_sep, dat_glm_oct, dat_glm_nov, dat_glm_dec)

dat_glm_annual <- data.frame(dat_glm_annual %>% group_by(cell_id) %>% 
                              summarize(x = unique(x), y = unique(y), ses = mean(pd.obs.z, na.rm = T), 
                                        p = mean(pd.obs.p, na.rm = T)
                              ))

df1 <- dat_glm_annual
df1$x.2 <- df1$x
df1[df1$x < 0 ,"x.2"] <- (df1[df1$x < 0 ,"x"]) + 360


g3 <- ggplot() + geom_raster(aes(x = x.2, y = y, fill = p), data = df1) +
  # play around with palettes. You can also use scale_fill_distiller() and
  geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
  coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                           panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
  theme(legend.text=element_text(size=8)) +
  scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
  scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
  scale_fill_viridis(na.value = "white") + labs(fill = "Freq. of signif. values:") 

setwd(wd_map)
ggsave(plot = g3, filename = "Freq_of_signif_values_ses_famd_ward_GLM.jpg", width  = 7, height = 4, dpi = 300)



# -------------------------------------------------------------------------------------------------
#ensemble

df2 <- dat_gam_annual
df2$x.2 <- df2$x
df2[df2$x < 0 ,"x.2"] <- (df2[df2$x < 0 ,"x"]) + 360
dat_gam <- df2[c("cell_id", "x.2", "y", "p")]
dat_gam$cell_id <- as.character(paste(dat_gam$x, dat_gam$y, sep = "_"))
colnames(dat_gam)[2] <- "x"
colnames(dat_gam)[4] <- "PD_GAM"


df2 <- dat_ann_annual
df2$x.2 <- df2$x
df2[df2$x < 0 ,"x.2"] <- (df2[df2$x < 0 ,"x"]) + 360
dat_ann <- df2[c("cell_id", "x.2", "y", "p")]
dat_ann$cell_id <- as.character(paste(dat_ann$x, dat_ann$y, sep = "_"))
colnames(dat_ann)[2] <- "x"
colnames(dat_ann)[4] <- "PD_ANN"


df2 <- dat_glm_annual
df2$x.2 <- df2$x
df2[df2$x < 0 ,"x.2"] <- (df2[df2$x < 0 ,"x"]) + 360
dat_glm <- df2[c("cell_id", "x.2", "y", "p")]
dat_glm$cell_id <- as.character(paste(dat_glm$x, dat_glm$y, sep = "_"))
colnames(dat_glm)[2] <- "x"
colnames(dat_glm)[4] <- "PD_GLM"


dat_ensm_part0  <- merge(data_clim_annual_new[c("cell_id","x","y")], dat_gam[c("cell_id","x","y", "PD_GAM")], all.x = T, by = c("cell_id","x","y"))
dat_ensm_part1  <- merge(dat_ensm_part0, dat_ann[c("cell_id","x","y", "PD_ANN")], all.x = T, by = c("cell_id","x","y"))
dat_ensm_part2  <- merge(dat_ensm_part1, dat_glm[c("cell_id","x","y", "PD_GLM")], all.x = T, by = c("cell_id","x","y"))




dat_esm2 <- data.frame(dat_ensm_part2 %>% group_by(cell_id) %>% 
                         summarize(x = unique(x), y = unique(y), p =  mean(c(PD_GAM,PD_GLM,PD_ANN), na.rm = F)))



df1 <- dat_esm2
df1$x.2 <- df1$x
df1[df1$x < 0 ,"x.2"] <- (df1[df1$x < 0 ,"x"]) + 360

g4 <- ggplot() + geom_raster(aes(x = x.2, y = y, fill = p), data = df1) +
  # play around with palettes. You can also use scale_fill_distiller() and
  geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
  coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                           panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
  theme(legend.text=element_text(size=8)) +
  scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
  scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
  scale_fill_viridis(na.value = "white") + labs(fill = "Freq. of signif. values:") 

setwd(wd_map)
sesFRsignif_data_06_04 <- df1
save(sesFRsignif_data_06_04, file = "sesFRsignif_data_06_04.RData")

ggsave(plot = g4, filename = "Freq_of_signif_values_ses_famd_ward_ensemble_narm.jpg", width  = 7, height = 4, dpi = 300)

# -------------------------------------------------------------------------------------------------
#test plots
# setwd(wd_dat_gam)
# 
# #ensemble
# load("ses_dat_faith_annual_GAM_binary_famd_ward_v2.Rda")
# dat_gam <- table_df;rm(table_df)
# dat_gam <- data.frame(dat_gam %>% group_by(cell_id) %>%
#                        summarize(x = unique(x), y = unique(y), ses = mean(pd.obs.z, na.rm = T),
#                                   p = mean(pd.obs.p, na.rm = T)
#                                                         ))
# 
# setwd(wd_dat_glm)
# load("ses_dat_faith_annual_GLM_binary_famd_ward_v2.Rda")
# dat_glm <- table_df;rm(table_df)
# dat_glm <- data.frame(dat_glm %>% group_by(cell_id) %>%
#                         summarize(x = unique(x), y = unique(y), ses = mean(pd.obs.z, na.rm = T),
#                                   p = mean(pd.obs.p, na.rm = T)
#                         ))
# 
# setwd(wd_dat_ann)
# load("ses_dat_faith_annual_ANN_binary_famd_ward_v2.Rda")
# dat_ann <- table_df;rm(table_df)
# dat_ann <- data.frame(dat_ann %>% group_by(cell_id) %>%
#                         summarize(x = unique(x), y = unique(y), ses = mean(pd.obs.z, na.rm = T),
#                                   p = mean(pd.obs.p, na.rm = T)
#                         ))
# 
# 
# 
# dat_ensm <- rbind(dat_gam, dat_glm, dat_ann)
# 
# #set p values below 0.05 to 1 and else to 0 to later calc a freq of signif
# dat_ensm$p <- ifelse(dat_ensm$p < 0.05, 1, 0)
# 
# 
# dat_ensm <- data.frame(dat_ensm %>% group_by(cell_id) %>%
#                         summarize(x = unique(x), y = unique(y), ses = mean(ses, na.rm = T),
#                                   p = mean(p, na.rm = T)
#                         ))
# 
# df2 <- dat_ensm
# df2$x.2 <- df2$x
# df2[df2$x < 0 ,"x.2"] <- (df2[df2$x < 0 ,"x"]) + 360
# 
# 
# g1 <- ggplot() + geom_raster(aes(x = x.2, y = y, fill = p), data = df2) +
#   # play around with palettes. You can also use scale_fill_distiller() and
#   geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
#   coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#                            panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
#   theme(legend.text=element_text(size=8)) +
#   scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
#   scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
#   scale_fill_viridis(na.value = "white") + labs(fill = "")


# # -------------------------------------------------------------------------------------------------
