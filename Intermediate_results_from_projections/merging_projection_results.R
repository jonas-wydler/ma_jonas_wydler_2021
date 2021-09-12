### 27/01/2021: R script to combine projection results, needed for definition of regions
# the second part of the script is used to plot the pca (figure 27 in thesis) that was used for the ocean regions in the thesis
#upated:12.9.2021
### ---------------------------------------------------------------------

Sys.setenv(LANGUAGE= 'en')

library("tidyverse") #needed
library("maps") #needed
library("FactoMineR") #needed
library("missMDA") #needed
library("ggrepel") #needed

### ---------------------------------------------------------------------
wd_script2 <- ("/data/manuscript_figure_4/") #directory where you want to save intermediate results
wd_script <- ("/data/") #directory where you stored the datasets for the ecoystem properties (e.g. fish catch data etc.)
wd_dat <- ("C:/Users/Jonas/polybox/part3/data_for_ensemble/richness/sumhsidplyrnarmf/")
wd_fatih <- ("C:/Users/Jonas/polybox/part3/2ndgo/")
wd_size <- ("C:/Users/Jonas/ma/sdm/part3/data_for_ensemble/size/")
wd_sac <- ("C:/Users/Jonas/ma/sdm/part3/data_for_ensemble/sac_groups/")
wd_myel <- ("C:/Users/Jonas/ma/sdm/part3/data_for_ensemble/myel_groups/")
wd_thropic <- ("C:/Users/Jonas/ma/sdm/part3/data_for_ensemble/trophic_groups/")
wd_feeding <- ("C:/Users/Jonas/ma/sdm/part3/data_for_ensemble/feeding_groups/")
wd_biomass <- ("C:/Users/Jonas/ma/sdm/part4/Zooplankton_biomass_COPEPOD/")
wd_clim <- ("/data/global_monthly_clims_1d/")

wd_newmaps2 <- ("C:/Users/Jonas/polybox/part4/redoforglmnarmf/")

### ---------------------------------------------------------------------
world2 <- map_data("world2") # world coastline for maps 
world1 <- map_data("world")


### ---------------------------------------------------------------------
#load species richness
setwd(wd_dat)
load("rich_GAM_dplyr_na_f.RData") 
ann_rich_GAM <- an.rich; rm(an.rich); gc();

load("rich_ANN_dplyr_na_f.RData")
ann_rich_ANN <- an.rich; rm(an.rich); gc();

load("rich_GLM_dplyr_na_f.RData")
ann_rich_GLM <- an.rich; rm(an.rich); gc();

an.rich <- rbind(ann_rich_GAM, ann_rich_ANN, ann_rich_GLM)

an.rich <- data.frame(an.rich %>% group_by(cell_id) %>% summarize(x = unique(x), y = unique(y), an.rich = mean(an.rich, na.rm = T)))
rm(ann_rich_GLM, ann_rich_ANN, ann_rich_GAM);gc();

### ---------------------------------------------------------------------
#fix the grid for the richness data
df <- an.rich
df$x.2 <- df$x
df[df$x < 0 ,"x.2"] <- (df[df$x < 0 ,"x"]) + 360

table_df_an <- df[c("cell_id","x.2","y","an.rich")]
colnames(table_df_an)[2] <-  "x"
table_df_an$cell_id <- factor(paste(table_df_an$x, table_df_an$y, sep = "_"))
data_rich <- table_df_an
#data_rich$log_rich <- log(data_rich$an.rich)
data_rich <- data_rich[c("cell_id", "x", "y", "an.rich")]

data_rich[sapply(data_rich, is.nan)] <- NA

rm(table_df_an)
rm(an.rich)
rm(df)
gc()
### ---------------------------------------------------------------------
#load functional richness estimate (according to Faith's index; normal and standard effect size (ses))
setwd(paste0(wd_fatih))
load("data_fd_fait_nonses_ward_ensembleNaRmF.RData")
#colnames(data_fd_fait_nonses_ward_ensemble)[5] <- "FR"
data_fr <- data_fd_fait_nonses_ward_ensemble
#data_fr$cell_id <- factor(paste(round(data_fr$x.2,1), round(data_fr$y,1), sep = "_"))
data_fr <- data_fr[c("cell_id", "x","y", "SR", "FR")]
#colnames(data_fr)[2] <- "x"
rm(data_fd_fait_nonses_ward_ensemble)

setwd(paste0(wd_fatih))
load("data_fd_fait_ses_ward_ensemble_naf.RData")
#colnames(data_fd_fait_ses_ward_ensemble)[4] <- "sesFR"
data_fr2 <- data_fd_fait_ses_ward_ensemble
#data_fr2$cell_id <- factor(paste(round(data_fr2$x.2,1), round(data_fr2$y,1), sep = "_"))
data_fr2 <- data_fr2[c("cell_id", "x","y", "sesFR")]
#colnames(data_fr2)[2] <- "x"
#rm(data_fd_fait_ses_ward_ensemble)



data_fr <- merge(data_fr,data_fr2, all.x = T, by = c("cell_id","x","y"))



### ---------------------------------------------------------------------
#load functional diversity indices FEve, FDis, RaoQ
setwd(wd_fatih)
load("data_fd_indices_ensemble_narmF.RData")
data_fd <- data_fd_indices_ensemble
data_fd <- data_fd[c("cell_id", "x", "y", "FEve", "FDis", "RaoQ")]
#data_fd$cell_id <- factor(paste(round(data_fd$x.2,1), round(data_fd$y,1), sep = "_"))
#colnames(data_fd)[2] <- "x"
rm(data_fd_indices_ensemble)
### ---------------------------------------------------------------------
#load CW median body size projections
setwd(wd_size)
load("annual_size_ANN_naF")
load("annual_size_GAM_naF")
load("annual_size_GLM_naF")

data_size <- rbind(annual_size_GAM_naF, annual_size_ANN_naF, annual_size_GLM_naF)

data_size <- data.frame(data_size %>% group_by(cell_id) %>% summarize(x = unique(x), y = unique(y), size = mean(med.size, na.rm = T)))

data_size$x.2 <- data_size$x
data_size[data_size$x < 0 ,"x.2"] <- (data_size[data_size$x < 0 ,"x"]) + 360

data_size$cell_id <- factor(paste(data_size$x.2, data_size$y, sep = "_"))
data_size <- data_size[c("cell_id", "x.2", "y", "size")]
colnames(data_size)[2] <- "x"

rm(annual_size_ANN_naF, annual_size_GAM_naF, annual_size_GLM_naF);gc();

### ---------------------------------------------------------------------
#load ecosystem proberties
setwd(wd_script)
load("data_particles.RData")

colnames(data_particles)
data_particles$log_Annual_NPP_v2 <- log(data_particles$Annual_NPP_v2)
data_particles$log_NPP <- log(data_particles$NPP)
data_particles$log_FPOC100m <- log(data_particles$FPOC100m)
data_particles$log_FPOCex <- log(data_particles$FPOCex)
colnames(data_particles)
data_particles <- data_particles[c("cell_id", "x", "y", "log_Annual_NPP_v2", "log_NPP", "log_FPOC100m", "log_FPOCex", "e")]

load("data_plankton.RData")

load("data_fish.RData")
colnames(data_fish)[4] <- "logged_fish_estimate"

load("data_biodiv.RData")
colnames(data_biodiv)[2] <- "x"

### ---------------------------------------------------------------------
#add cwm proportion of sac spawning 
setwd(wd_sac)
load("df_sac_ANN.Rda")
load("df_sac_GAM.Rda")
load("df_sac_GLM.Rda")

data_sac <- rbind(df_sac_GAM, df_sac_ANN, df_sac_GLM)

data_sac <- data.frame(data_sac %>% group_by(cell_id) %>% summarize(x = unique(x), y = unique(y), sac_prop = mean(relfrq, na.rm = T)))

data_sac$x.2 <- data_sac$x
data_sac[data_sac$x < 0 ,"x.2"] <- (data_sac[data_sac$x < 0 ,"x"]) + 360

data_sac$cell_id <- factor(paste(data_sac$x.2, data_sac$y, sep = "_"))
data_sac <- data_sac[c("cell_id", "x.2", "y", "sac_prop")]
colnames(data_sac)[2] <- "x"

rm(df_sac_ANN, df_sac_GAM, df_sac_GLM);gc();

### ---------------------------------------------------------------------
#add cwm proportion of myelination
setwd(wd_myel)
load("df_myel_ANN.Rda")
load("df_myel_GAM.Rda")
load("df_myel_GLM.Rda")

data_myel <- rbind(df_myel_GAM, df_myel_ANN, df_myel_GLM)

data_myel <- data.frame(data_myel %>% group_by(cell_id) %>% summarize(x = unique(x), y = unique(y), myel_prop = mean(relfrq, na.rm = T)))

data_myel$x.2 <- data_myel$x
data_myel[data_myel$x < 0 ,"x.2"] <- (data_myel[data_myel$x < 0 ,"x"]) + 360

data_myel$cell_id <- factor(paste(data_myel$x.2, data_myel$y, sep = "_"))
data_myel <- data_myel[c("cell_id", "x.2", "y", "myel_prop")]
colnames(data_myel)[2] <- "x"


rm(df_myel_ANN, df_myel_GAM, df_myel_GLM);gc();
### ---------------------------------------------------------------------
#add cwm proportion of feeding modes
setwd(wd_feeding)
#ambush
load("df_ambush_GAM.Rda")
load("df_ambush_GLM.Rda")
load("df_ambush_ANN.Rda")

data_ambush <- rbind(df_ambush_GAM, df_ambush_ANN, df_ambush_GLM)

data_ambush <- data.frame(data_ambush %>% group_by(cell_id) %>% summarize(x = unique(x), y = unique(y), ambush_prop = mean(relfrq, na.rm = T)))

data_ambush$x.2 <- data_ambush$x
data_ambush[data_ambush$x < 0 ,"x.2"] <- (data_ambush[data_ambush$x < 0 ,"x"]) + 360

data_ambush$cell_id <- factor(paste(data_ambush$x.2, data_ambush$y, sep = "_"))
data_ambush <- data_ambush[c("cell_id", "x.2", "y", "ambush_prop")]
colnames(data_ambush)[2] <- "x"


rm(df_ambush_ANN, df_ambush_GAM, df_ambush_GLM);gc();

#add cwm proportion of cruise feeding
load("df_cruise_ANN.Rda")
load("df_cruise_GAM.Rda")
load("df_cruise_GLM.Rda")

data_cruise <- rbind(df_cruise_GAM, df_cruise_ANN, df_cruise_GLM)

data_cruise <- data.frame(data_cruise %>% group_by(cell_id) %>% summarize(x = unique(x), y = unique(y), cruise_prop = mean(relfrq, na.rm = T)))

data_cruise$x.2 <- data_cruise$x
data_cruise[data_cruise$x < 0 ,"x.2"] <- (data_cruise[data_cruise$x < 0 ,"x"]) + 360

data_cruise$cell_id <- factor(paste(data_cruise$x.2, data_cruise$y, sep = "_"))
data_cruise <- data_cruise[c("cell_id", "x.2", "y", "cruise_prop")]
colnames(data_cruise)[2] <- "x"


rm(df_cruise_ANN, df_cruise_GAM, df_cruise_GLM);gc();

#add cwm proportion of current-ambush feeding
load("df_current_ambush_GAM.Rda")
load("df_current_ambush_GLM.Rda")
load("df_current_ambush_ANN.Rda")

data_current_ambush <- rbind(df_current_ambush_GAM, df_current_ambush_ANN, df_current_ambush_GLM)

data_current_ambush <- data.frame(data_current_ambush %>% group_by(cell_id) %>% summarize(x = unique(x), y = unique(y), current_ambush_prop = mean(relfrq, na.rm = T)))

data_current_ambush$x.2 <- data_current_ambush$x
data_current_ambush[data_current_ambush$x < 0 ,"x.2"] <- (data_current_ambush[data_current_ambush$x < 0 ,"x"]) + 360

data_current_ambush$cell_id <- factor(paste(data_current_ambush$x.2, data_current_ambush$y, sep = "_"))
data_current_ambush <- data_current_ambush[c("cell_id", "x.2", "y", "current_ambush_prop")]
colnames(data_current_ambush)[2] <- "x"


rm(df_current_ambush_ANN, df_current_ambush_GAM, df_current_ambush_GLM);gc();

#add cwm proportion of current feeding
load("df_current_ANN.Rda")
load("df_current_GAM.Rda")
load("df_current_GLM.Rda")

data_current <- rbind(df_current_GAM, df_current_ANN, df_current_GLM)

data_current <- data.frame(data_current %>% group_by(cell_id) %>% summarize(x = unique(x), y = unique(y), current_prop = mean(relfrq, na.rm = T)))

data_current$x.2 <- data_current$x
data_current[data_current$x < 0 ,"x.2"] <- (data_current[data_current$x < 0 ,"x"]) + 360

data_current$cell_id <- factor(paste(data_current$x.2, data_current$y, sep = "_"))
data_current <- data_current[c("cell_id", "x.2", "y", "current_prop")]
colnames(data_current)[2] <- "x"


rm(df_current_ANN, df_current_GAM, df_current_GLM);gc();

#add cwm proportion of Current-cruise feeding
load("df_current_cruise_ANN.Rda")
load("df_current_cruise_GAM.Rda")
load("df_current_cruise_GLM.Rda")

data_current_cruise <- rbind(df_current_cruise_GAM, df_current_cruise_ANN, df_current_cruise_GLM)

data_current_cruise <- data.frame(data_current_cruise %>% group_by(cell_id) %>% summarize(x = unique(x), y = unique(y), current_cruise_prop = mean(relfrq, na.rm = T)))

data_current_cruise$x.2 <- data_current_cruise$x
data_current_cruise[data_current_cruise$x < 0 ,"x.2"] <- (data_current_cruise[data_current_cruise$x < 0 ,"x"]) + 360

data_current_cruise$cell_id <- factor(paste(data_current_cruise$x.2, data_current_cruise$y, sep = "_"))
data_current_cruise <- data_current_cruise[c("cell_id", "x.2", "y", "current_cruise_prop")]
colnames(data_current_cruise)[2] <- "x"


rm(df_current_cruise_ANN, df_current_cruise_GAM, df_current_cruise_GLM);gc();

#add cwm proportion of Current-cruise feeding
load("df_particle_GLM.Rda")
load("df_particles_ANN.Rda")
load("df_particles_GAM.Rda")

data_particle_feeding <- rbind(df_particle_GLM, df_particles_ANN, df_particles_GAM)

data_particle_feeding <- data.frame(data_particle_feeding %>% group_by(cell_id) %>% summarize(x = unique(x), y = unique(y), particle_feeding = mean(relfrq, na.rm = T)))

data_particle_feeding$x.2 <- data_particle_feeding$x
data_particle_feeding[data_particle_feeding$x < 0 ,"x.2"] <- (data_particle_feeding[data_particle_feeding$x < 0 ,"x"]) + 360

data_particle_feeding$cell_id <- factor(paste(data_particle_feeding$x.2, data_particle_feeding$y, sep = "_"))
data_particle_feeding <- data_particle_feeding[c("cell_id", "x.2", "y", "particle_feeding")]
colnames(data_particle_feeding)[2] <- "x"


rm(df_particle_GLM, df_particles_ANN, df_particles_GAM);gc();


### ---------------------------------------------------------------------
#add cwm proportion of carnivores
setwd(wd_thropic)
load("df_carni_GAM.Rda")
df_carni_GAM <-  df_carni
load("df_carni_GLM.Rda")
df_carni_GLM <-  df_carni

load("df_carni_ANN.Rda")
df_carni_ANN <-  df_carni


data_carni <- rbind(df_carni_GAM, df_carni_GLM, df_carni_ANN)
data_carni <- data.frame(data_carni %>% group_by(cell_id) %>% summarize(x = unique(x), y = unique(y), carni_prop = mean(relfrq, na.rm = T)))

data_carni$x.2 <- data_carni$x
data_carni[data_carni$x < 0 ,"x.2"] <- (data_carni[data_carni$x < 0 ,"x"]) + 360

data_carni$cell_id <- factor(paste(data_carni$x.2, data_carni$y, sep = "_"))
data_carni <- data_carni[c("cell_id", "x.2", "y", "carni_prop")]
colnames(data_carni)[2] <- "x"
rm(df_carni_GAM, df_carni_GLM, df_carni_ANN);gc();

#add cwm proportion of herbivores
setwd(wd_thropic)
load("df_herbi_GAM.Rda")
df_herbi_GAM <-  df_herbi

load("df_herbi_GLM.Rda")
df_herbi_GLM <-  df_herbi

load("df_herbi_ANN.Rda")
df_herbi_ANN <-  df_herbi


data_herbi <- rbind(df_herbi_GAM, df_herbi_GLM, df_herbi_ANN)
data_herbi <- data.frame(data_herbi %>% group_by(cell_id) %>% summarize(x = unique(x), y = unique(y), herbi_prop = mean(relfrq, na.rm = T)))

data_herbi$x.2 <- data_herbi$x
data_herbi[data_herbi$x < 0 ,"x.2"] <- (data_herbi[data_herbi$x < 0 ,"x"]) + 360

data_herbi$cell_id <- factor(paste(data_herbi$x.2, data_herbi$y, sep = "_"))
data_herbi <- data_herbi[c("cell_id", "x.2", "y", "herbi_prop")]
colnames(data_herbi)[2] <- "x"


rm(df_herbi_GAM, df_herbi_GLM, df_herbi_ANN);gc();

#add cwm proportion of detritivores
setwd(wd_thropic)
load("df_detri_GAM.Rda")
df_detri_GAM <-  df_detri

load("df_detri_GLM.Rda")
df_detri_GLM <-  df_detri

load("df_detri_ANN.Rda")
df_detri_ANN <-  df_detri


data_detri <- rbind(df_detri_GAM, df_detri_GLM, df_detri_ANN)
data_detri <- data.frame(data_detri %>% group_by(cell_id) %>% summarize(x = unique(x), y = unique(y), detri_prop = mean(relfrq, na.rm = T)))

data_detri$x.2 <- data_detri$x
data_detri[data_detri$x < 0 ,"x.2"] <- (data_detri[data_detri$x < 0 ,"x"]) + 360

data_detri$cell_id <- factor(paste(data_detri$x.2, data_detri$y, sep = "_"))
data_detri <- data_detri[c("cell_id", "x.2", "y", "detri_prop")]
colnames(data_detri)[2] <- "x"


rm(df_detri_GAM, df_detri_GLM, df_detri_ANN, df_detri, df_herbi, df_carni);gc();


#add cwm proportion of omnivores
setwd(wd_thropic)
load("df_omni_GAM.Rda")
df_omni_GAM <-  df_omni

load("df_omni_GLM.Rda")
df_omni_GLM <-  df_omni

load("df_omni_ANN.Rda")
df_omni_ANN <-  df_omni


data_omni <- rbind(df_omni_GAM, df_omni_GLM, df_omni_ANN)
data_omni <- data.frame(data_omni %>% group_by(cell_id) %>% summarize(x = unique(x), y = unique(y), omni_prop = mean(relfrq, na.rm = T)))

data_omni$x.2 <- data_omni$x
data_omni[data_omni$x < 0 ,"x.2"] <- (data_omni[data_omni$x < 0 ,"x"]) + 360

data_omni$cell_id <- factor(paste(data_omni$x.2, data_omni$y, sep = "_"))
data_omni <- data_omni[c("cell_id", "x.2", "y", "omni_prop")]
colnames(data_omni)[2] <- "x"


rm(df_omni_GAM, df_omni_GLM, df_omni_ANN, df_omni);gc();
### ---------------------------------------------------------------------
#add zooplankton biomass estimate
setwd(wd_biomass)
data_biomass <- read.csv2("table_all_zooplankton_biomasses_fields_1d_COPEPOD-MAREDAT_20_01_21.txt", h = T, sep = "\t")
data_biomass$x <- as.numeric(data_biomass$x)
data_biomass$y <- as.numeric(data_biomass$y)
data_biomass$Biomass <- as.numeric(data_biomass$Biomass)

data_biomass <- data.frame(data_biomass %>% group_by(id) %>% summarize(x = unique(x), y = unique(y), biomass = mean(Biomass, na.rm = T)))

data_biomass$x.2 <- data_biomass$x
data_biomass[data_biomass$x < 0 ,"x.2"] <- (data_biomass[data_biomass$x < 0 ,"x"]) + 360

data_biomass$cell_id <- factor(paste(data_biomass$x.2, data_biomass$y, sep = "_"))
data_biomass <- data_biomass[c("cell_id", "x.2", "y", "biomass")]
colnames(data_biomass)[2] <- "x"
### ---------------------------------------------------------------------
#add env variables
setwd(wd_clim)

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


### ---------------------------------------------------------------------
#MERGE
data_total <- data_rich
summary(data_total)

data_total  <- merge(data_rich,data_fd, all.x = T, by = c("cell_id","x","y"))
data_total  <- merge(data_total,data_fr, all.x = T, by = c("cell_id","x","y"))
data_total2 <- merge(data_total,data_particles, all.x = T, by = c("cell_id","x","y"))
data_total3 <- merge(data_total2,data_plank, all.x = T, by = c("cell_id","x","y"))
data_total4 <- merge(data_total3,data_size, all.x = T, by = c("cell_id","x","y"))

data_total5 <- merge(data_total4,data_sac, all.x = T, by = c("cell_id","x","y"))
data_total5 <- merge(data_total5,data_myel, all.x = T, by = c("cell_id","x","y"))

data_total6 <- merge(data_total5,data_ambush, all.x = T, by = c("cell_id","x","y"))
data_total6 <- merge(data_total6,data_cruise, all.x = T, by = c("cell_id","x","y"))
data_total6 <- merge(data_total6,data_current, all.x = T, by = c("cell_id","x","y"))
data_total6 <- merge(data_total6,data_current_ambush, all.x = T, by = c("cell_id","x","y"))
data_total6 <- merge(data_total6,data_current_cruise, all.x = T, by = c("cell_id","x","y"))

data_total7 <- merge(data_total6,data_carni, all.x = T, by = c("cell_id","x","y"))
data_total7 <- merge(data_total7,data_herbi, all.x = T, by = c("cell_id","x","y"))
data_total7 <- merge(data_total7,data_detri, all.x = T, by = c("cell_id","x","y"))
data_total7 <- merge(data_total7,data_omni, all.x = T, by = c("cell_id","x","y"))

#also merge ecosystem properties
data_total8 <- merge(data_total7,data_fish, all.x = T, by = c("cell_id","x","y"))
data_total9 <- merge(data_total8,data_biomass, all.x = T, by = c("cell_id","x","y"))
data_total10 <- merge(data_total9,data_biodiv, all.x = T, by = c("cell_id","x","y"))


data_total_clim <- merge(data_total10,data_clim_annual_new, all.x = T, by = c("cell_id","x","y"))
data_total_everything <- data_total_clim
setwd(wd_script2)
colnames(data_total_everything)[17] <- "body_size"
colnames(data_total_everything)[4] <- "SR"
colnames(data_total_everything)[8] <- "SR_binary"

save(data_total_everything, file = "data_total_for_boxplots.RData")

# ### ---------------------------------------------------------------------
# #pca for ocean regions in thesis (figure 15)
setwd(wd_script2)
load("data_total_for_boxplots.RData")

data_pca_total <- na.omit(data_total_everything[c(1:28,35:47)])
colnames(data_pca_total) #check these
#data_pca_traits <- data_pca_total[c(4:7,9,17:28)]

ncp <- estim_ncp(na.omit(data_pca_traits))
plot(ncp$criterion ~ seq(1:17), ylab = "MSEP criterion", xlab = "the number of components retained for the PCA")
colnames(data_pca_total[2:41])
colnames(data_pca_total)
data_pca_total <- data_pca_total[-c(7,11,12,13,14,15,16)]
#data_pca_total <- data_pca_total[-c(7,11,13)]

colnames(data_pca_total[2:length(data_pca_total)])

data_pca_total2 <- data_pca_total[c(4:33)]
pca1 <- PCA(data_pca_total2, scale.unit = T, ncp = 5, graph = T,  quanti.sup = c(4,6,19:30))#check these
#pca1 <- PCA(data_pca_traits, scale.unit = T, ncp = 5, graph = T)#check these

pca_coords <- as.data.frame(pca1$ind$coord)
pca_coords <- cbind(data_pca_total[1:3],pca_coords[c(1:5)])

# ### ---------------------------------------------------------------------
### Better PCA plot from fabio
# Extract data from a FactoMineR::PCA object
augment.PCA <- function(x, dims = c(1:5), which="col") {
  .get <- function(x, element, dims) {
    y <- as.data.frame(x[[element]]$coord[,dims])
    if (nrow(y) == 0) {
      y <- NULL
    } else {
      y$type <- element
    }
    return(y)
  }
  if (which == "col") {
    y <- rbind(.get(x, "var", dims), .get(x, "quanti.sup", dims))
  } else {
    y <- rbind(.get(x, "ind", dims), .get(x, "quali.sup", dims))
  }
  y$var <- row.names(y)
  row.names(y) <- NULL
  return(y)
}

pcad <- augment.PCA(pca1)
pcad$var
pcad$type #
pcad$type2 <- NA
pcad[pcad$type == 'var',"type2"] <- "vars"
pcad[pcad$var %in% c("SST","logNO3","logChl","Sistar","Nstar","PAR","logEKE", "Wind", "MLD", "pCO2"),"type2"] <- "Predictor"
pcad[pcad$var %in% c("e","slope_psi","dO2","log_NPP","log_FPOCex", "latitude","x", "y","SSS","logPO4","sesFR","SR_binary"),"type2"] <- "Covariate"

p1 <- ggplot(pcad) + coord_fixed() + scale_x_continuous(breaks=0) + scale_y_continuous(breaks=0) +
  annotate(geom="path", colour="black", x=cos(seq(0, 2*pi, length.out=100)), y=sin(seq(0, 2*pi, length.out=100)), colour = "grey30") +
  geom_segment(aes(x=0, xend = Dim.1, y=0, yend = Dim.2, color = factor(type2)), arrow=arrow(angle=20, length=unit(0.01, "npc"))) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") + geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") + 
  scale_colour_manual(name = "", values = c("#b2182b","#2166ac","black"), guide = F) + 
  geom_text_repel(aes(x=Dim.1, y=Dim.2, color=type2, label=var), 
                  data=pcad, segment.alpha= 0.5) +
  xlab("pca1 (64.20%)") + ylab("pca2 (19.51%)") + theme_bw()

setwd(wd_newmaps2)
ggsave(plot = p1, filename = paste("pca_1_2_v3.jpg"),dpi = 600, width = 10, height = 10)


p2 <- ggplot(pcad) + coord_fixed() + scale_x_continuous(breaks=0) + scale_y_continuous(breaks=0) +
  annotate(geom="path", colour="black", x=cos(seq(0, 2*pi, length.out=100)), y=sin(seq(0, 2*pi, length.out=100)), colour = "grey30") +
  geom_segment(aes(x=0, xend = Dim.2, y=0, yend = Dim.3, colour = factor(type2)), arrow=arrow(angle=20, length=unit(0.01,"npc"))) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") + geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") + 
  scale_colour_manual(name = "", values = c("#b2182b","#2166ac","black"), guide = F) + 
  geom_text_repel(aes(x=Dim.2, y=Dim.3, colour = type2, label = var), 
                  data=pcad, segment.alpha=0.5) +
  xlab("pca2 (19.4%)") + ylab("pca3 (6.17%)") + theme_bw()

setwd(wd_newmaps2)
ggsave(plot = p2, filename = paste("pca_2_3_v2.jpg"),dpi = 600, width = 10, height = 10)

p3 <- ggplot(pcad) + coord_fixed() + scale_x_continuous(breaks=0) + scale_y_continuous(breaks=0) +
  annotate(geom="path", colour="black", x=cos(seq(0, 2*pi, length.out=100)), y=sin(seq(0, 2*pi, length.out=100)), colour = "grey30") +
  geom_segment(aes(x=0, xend = Dim.3, y=0, yend = Dim.4, colour = factor(type2)), arrow=arrow(angle=20, length=unit(0.01,"npc"))) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") + geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") + 
  scale_colour_manual(name = "", values = c("#b2182b","#2166ac","black"), guide = F) + 
  geom_text_repel(aes(x=Dim.3, y=Dim.4, colour = type2, label = var), 
                  data=pcad, segment.alpha=0.5) +
  xlab("pca3 (6.17%)") + ylab("pca4 (5.17)%") + theme_bw()

setwd(wd_newmaps2)
ggsave(plot = p3, filename = paste("pca_3_4_v2.jpg"),dpi = 600, width = 10, height = 10)

p4 <- ggplot(pcad) + coord_fixed() + scale_x_continuous(breaks=0) + scale_y_continuous(breaks=0) +
  annotate(geom="path", colour="black", x=cos(seq(0, 2*pi, length.out=100)), y=sin(seq(0, 2*pi, length.out=100)), colour = "grey30") +
  geom_segment(aes(x=0, xend = Dim.4, y=0, yend = Dim.5, colour = factor(type2)), arrow=arrow(angle=20, length=unit(0.01,"npc"))) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") + geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") + 
  scale_colour_manual(name = "", values = c("#b2182b","#2166ac","black"), guide = F) + 
  geom_text_repel(aes(x=Dim.4, y=Dim.5, colour = type2, label = var), 
                  data=pcad, segment.alpha=0.5) +
  xlab("pca4 (5.17%)") + ylab("pca5 (2.19%)") + theme_bw()

setwd(wd_newmaps2)
ggsave(plot = p4, filename = paste("pca_4_5_v2.jpg"),dpi = 600, width = 10, height = 10)
