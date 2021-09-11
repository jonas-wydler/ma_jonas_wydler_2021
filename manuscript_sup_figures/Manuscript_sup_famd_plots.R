# -------------------------------------------------------------------------------------------------
#FAMD plots, Jonas Wydler, 03/03/2021
# - Redo famd and create vectorimages for manuscript
# -------------------------------------------------------------------------------------------------
Sys.setenv(LANGUAGE= 'en')

library("tidyverse")
library("vegan")
library("FactoMineR")
library("factoextra")
#library("clustMixType")
library("missMDA")
library("ggrepel")
library("broom")
library("factoextra")
library("ggpubr")
library("svglite")
### -----------------------------------------------------------------------
wd_dat <- ("/data/") #where the trait table is
wd_plots <- ("/data/") #where you want the plots to go
### -----------------------------------------------------------------------
#load data
setwd(wd_dat)

fct <- read.csv("table_funct_traits_copepods_v2.csv", h = T, sep = ";", dec = ",")
fct <- fct[,c(3:20)]

### Count NA
names <- colnames(fct)[c(7,8,9,10,15)] ; names
fct$na_count <- apply(fct[,names], 1, function(x) sum(is.na(x)))
# Drop species with missing body length or 2 or more missing traits
fct <- fct[!is.na(fct$max_body_length),]
fct <- fct[fct$na_count < 2,]
### -----------------------------------------------------------------------

fct <- fct[]

#saving as factors for FAMD
fct$Spawning <- as.factor(fct$Spawning)
fct$Myelination <- as.factor(fct$Myelination)
fct$Omnivore <- as.factor(fct$Omnivore)
fct$Carnivore <- as.factor(fct$Carnivore)
fct$Herbivore <- as.factor(fct$Herbivore)
fct$Detritivore <- as.factor(fct$Detritivore)
fct$Current <- as.factor(fct$Current)
fct$Cruise <- as.factor(fct$Cruise)
fct$Ambush <- as.factor(fct$Ambush)
fct$Trophism <- as.factor(fct$Trophism)
fct$Feeding_mode <- as.factor(fct$Feeding_mode)
fct$Trophism <- as.factor(fct$Trophism)
fct$Feeding_mode <- as.factor(fct$Feeding_mode)
### -----------------------------------------------------------------------
#(298/385,3)*100

traits <- fct[,c(7:9,11:14,16:18)]
#traits <- fct[,c("max_body_length", "Myelination", "Spawning","Trophism", "Feeding_mode")]
colnames(traits)[1] <- "Body_size"

colnames(traits) <- c("Body_size", "Myelination", "Spawning_mode", "Omnivory", "Carnivory", "Herbivory", "Detritivory",
"Current_feeding", "Cruise_feeding", "Ambush_feeding")


### -----------------------------------------------------------------------

#FAMD <- FAMD(na.omit(traits), graph = T)
#FAMD <- FAMD(traits, graph = T)

#determine optimal number of retained pcs
#npfamd <- estim_ncpFAMD(as.data.frame(traits))

#imputation of missing variables
compfamd <- imputeFAMD(traits, npc = 4) #npc here can be npfamd from above
FAMD <- FAMD(traits, ncp = 4, tab.disj = compfamd$tab.disj, graph = T)

famd_sp <- data.frame(FAMD$ind$coord[,1:4])
colnames(famd_sp) <- c("FAMD1","FAMD2","FAMD3","FAMD4")

eigFAMD<- data.frame(perc = FAMD$eig[,"percentage of variance"], nb = c(1:nrow(FAMD$eig)) )

famd1 <- paste0("FAMD 1 (",floor(eigFAMD$perc[1]*100)/100,"%)")
famd2 <- paste0("FAMD 2 (",floor(eigFAMD$perc[2]*100)/100,"%)")
famd3 <- paste0("FAMD 3 (",floor(eigFAMD$perc[3]*100)/100,"%)")
famd4 <- paste0("FAMD 4 (",floor(eigFAMD$perc[4]*100)/100,"%)")

fviz_famd_var(FAMD, repel = TRUE, axes = c(1,2))
fviz_famd_var(FAMD, repel = TRUE, axes = c(3,4))

round(FAMD$var$contrib,1)

FAMD$var
FAMD$eig

#save plots
setwd(wd_plots)
ggsave(plot = fviz_famd_var(FAMD, repel = TRUE, axes = c(1,2)), filename = paste0("FG_famd_1_2_var_plot.svg"),width = 7, height = 7, dpi = 300)
ggsave(plot = fviz_famd_var(FAMD, repel = TRUE, axes = c(3,4)), filename = paste0("FG_famd_3_4_var_plot.svg"),width = 7, height = 7, dpi = 300)





