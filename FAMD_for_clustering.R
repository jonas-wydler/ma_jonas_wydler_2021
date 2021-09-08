#Clustering functional traits using FAMD, 12.10.2020, Jonas Wydler

setwd("wd_data")
Sys.setenv(LANG = "en")

# --------------------------------------------------------------------------------------------------------------------
library("tidyverse")
library("FactoMineR")
library("factoextra")
library("clustMixType")
library("missMDA")

# --------------------------------------------------------------------------------------------------------------------
### Load the functional traits table and explore
fct <- read.csv("table_funct_traits_copepods_v2.csv", h = T, sep = ";", dec = ",")
dim(fct) ; head(fct);  colnames(fct)
str(fct)
fct <- fct[,c(3:20)]

### Count NA
names <- colnames(fct)[c(7,8,9,10,15)] ; names
fct$na_count <- apply(fct[,names], 1, function(x) sum(is.na(x)))
# Drop species with missing body length and more than two missing traits
fct <- fct[!is.na(fct$max_body_length),]
fct <- fct[fct$na_count < 2,]
# ---------------------------------------------------------------------------------------------------------------------
###FAMD

#not sure why we do this
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



#Determining nr of dimensons
npfamd <- estim_ncpFAMD(fct[,c(7:9,11:14,16:18)])

#FAMD with missing values
compfamd <- imputeFAMD(fct[,c(7:9,11:14,16:18)], npc = 4)
FAMD <- FAMD(fct[,c(7:9,11:14,16:18)], tab.disj = compfamd$tab.disj)
#FAMD <- FAMD(fct[,c(7:9,11:14,16:18)])
eigFAMD<- data.frame(perc = FAMD$eig[,"percentage of variance"], nb = c(1:nrow(FAMD$eig)) )


# --------------------------------------------------------------------------------------------------------------------

famd1 <- paste0("FAMD 1 (",floor(eigFAMD$perc[1]*100)/100,"%)")
famd2 <- paste0("FAMD 2 (",floor(eigFAMD$perc[2]*100)/100,"%)")
famd3 <- paste0("FAMD 3 (",floor(eigFAMD$perc[3]*100)/100,"%)")
famd4 <- paste0("FAMD 4 (",floor(eigFAMD$perc[4]*100)/100,"%)")


famd <- data.frame(FAMD$ind$coord[,1:4])


famd_sp <- data.frame(FAMD$ind$coord[,1:4])
colnames(famd_sp) <- c("FAMD1","FAMD2","FAMD3","FAMD4")
famd_all_temp <- rbind(famd_sp)

famd_all_sp <- famd_all_temp[ order(row.names(famd_all_temp)), ]

famd_dist <- dist(famd_sp, method = "euclidean")
save(famd_sp, file = "famd_withNA.RData")

plot(hclust(famd_dist,"average"), hang = -1)
