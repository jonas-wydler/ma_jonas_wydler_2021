### -----------------------------------------------------------------------
#03.06.2021 Check impact of dimension reduction according to methods by Mouillot et al. 2021
#jonas wydler
### -----------------------------------------------------------------------
library("dimRed")
library("coRanking")
library("tidyverse")
library("missMDA")
library("FactoMineR")
library("ape")
library("FD")
Sys.setenv(LANGUAGE= 'en')
### -----------------------------------------------------------------------
wd_trait_dat <- ("wd_data")
### -----------------------------------------------------------------------
#load data
setwd(wd_dat)

fct <- read.csv("table_funct_traits_copepods_v2.csv", h = T, sep = ";", dec = ",")
fct <- fct[,c(3:20)]

### Count NA
names <- colnames(fct)[c(7,8,9,10,15)] ; names
fct$na_count <- apply(fct[,names], 1, function(x) sum(is.na(x)))
# Drop some NA, discuss with Fabio
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
#select the traits for the analysis
#all
traits <- fct[,c(7:9,11:14,16:18)]

#categories only to investigate trait omission
#traits <- fct[,c("max_body_length", "Myelination", "Spawning","Trophism", "Feeding_mode")]
colnames(traits)[1] <- "Body_size"

colnames(traits)
#
#traits2 <- traits[,c("Myelination", "Spawning", "Trophism","Feeding_mode")]
#traits2 <- traits[,c("Spawning", "Myelination", "Omnivore", "Carnivore", 
                     #"Herbivore","Detritivore", "Current", "Cruise", "Ambush")]
#colnames(traits2)

#colnames(traits) <- c("Body_size", "Myelination", "Spawning_mode", "Omnivory", "Carnivory", "Herbivory", "Detritivory",
#"Current_feeding", "Cruise_feeding", "Ambush_feeding")
### -----------------------------------------------------------------------
#high dimensional data
gow_dis <- gowdis(traits)


### -----------------------------------------------------------------------
#low dimensional data


#FAMD <- FAMD(na.omit(traits), graph = T)
#FAMD <- FAMD(traits, graph = T)

npfamd <- estim_ncpFAMD(as.data.frame(traits))
np <- 4
compfamd <- imputeFAMD(traits, npc = np)
FAMD <- FAMD(traits, ncp = np, graph = F)
FAMD <- FAMD(traits, ncp = np, tab.disj = compfamd$tab.disj, graph = T)

famd_sp <- data.frame(FAMD$ind$coord[,1:np])
colnames(famd_sp) <- c("FAMD1","FAMD2","FAMD3","FAMD4")[1:np]

eigFAMD<- data.frame(perc = FAMD$eig[,"percentage of variance"], nb = c(1:nrow(FAMD$eig)))

famd1 <- paste0("FAMD 1 (",floor(eigFAMD$perc[1]*100)/100,"%)")
famd2 <- paste0("FAMD 2 (",floor(eigFAMD$perc[2]*100)/100,"%)")
famd3 <- paste0("FAMD 3 (",floor(eigFAMD$perc[3]*100)/100,"%)")
famd4 <- paste0("FAMD 4 (",floor(eigFAMD$perc[4]*100)/100,"%)")

famd_dist <- dist(famd_sp, method = "euclidean")

#MCA (for case without body size)
# npfamd <- estim_ncpMCA(as.data.frame(traits2))
# compmca <- imputeMCA(don = traits2, ncp = 4)
# MCA <- MCA(X = traits2, ncp = 4, graph = F, tab.disj = compmca$tab.disj)
# mca_sp <- data.frame(MCA$ind$coord[,1:4])
# colnames(mca_sp) <- c("MCA1","MCA2","MCA3","MCA4")[1:4]
# 
# eigMCA<- data.frame(perc = MCA$eig[,"percentage of variance"], nb = c(1:nrow(MCA$eig)) )
# 
# mca1 <- paste0("MCA 1 (",floor(eigMCA$perc[1]*100)/100,"%)")
# mca2 <- paste0("MCA 2 (",floor(eigMCA$perc[2]*100)/100,"%)")
# mca3 <- paste0("MCA 3 (",floor(eigMCA$perc[3]*100)/100,"%)")
# mca4 <- paste0("MCA 4 (",floor(eigMCA$perc[4]*100)/100,"%)")
# 
# 
# mca_dist <- dist(mca_sp, method = "euclidean")

### -----------------------------------------------------------------------
#Method by Mouillot et al. 2021
#calculate coranking matrix
Q = coranking(Xi = gow_dis, X = famd_dist, input_Xi = c("dist"), input_X = c("dist"), use = "R")
coRanking::imageplot(Q)
NX        <- coRanking::R_NX(Q)
#calculate AUC
AUC       <- coRanking::AUC_ln_K(NX)
AUC


