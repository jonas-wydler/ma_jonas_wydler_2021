#Jonas Wydler;  dendrogram for paper; 23.05.2021
#Careful!! For the plot in the manuscript we changed the species names manually 
#according to scheme described in FG_name_change.txt
#
#-----------------------------------------------------------------
wd_trait_dat <- ("wd_data")
Sys.setenv(LANG = "en")
#-----------------------------------------------------------------
library('dendextend')
library("tidyverse")
library("FactoMineR")
library("missMDA")
library("DendSer ")
library("svglite")
library("RColorBrewer")
#-----------------------------------------------------------------
#Read in trait data
setwd(wd_trait_dat)
fct <- read.csv("table_funct_traits_copepods_v2.csv", h = T, sep = ";", dec = ",")
fct <- fct[,c(3:20)]
names <- colnames(fct)[c(7,8,9,10,15)] ; names
fct$na_count <- apply(fct[,names], 1, function(x) sum(is.na(x)))
#-----------------------------------------------------------------
# Drop species with missing body size info
fct <- fct[!is.na(fct$max_body_length),]
#-----------------------------------------------------------------
# Drop species with more than two missing traits
fct <- fct[fct$na_count < 2,]
#-----------------------------------------------------------------
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
#-----------------------------------------------------------------
#FAMD
compfamd <- imputeFAMD(fct[,c(7:9,11:14,16:18)], npc = 4)
FAMD <- FAMD(fct[,c(7:9,11:14,16:18)], tab.disj = compfamd$tab.disj, graph = F)
famd <- data.frame(FAMD$ind$coord[,1:4])
famd_sp <- data.frame(FAMD$ind$coord[,1:4])
colnames(famd_sp) <- c("FAMD1","FAMD2","FAMD3","FAMD4")
famd_all_temp <- rbind(famd_sp)
famd_all_sp <- famd_all_temp[ order(row.names(famd_all_temp)), ]
famd_dist <- dist(famd_all_temp, method = "euclidean")
#-----------------------------------------------------------------
#Clustering
fit_famd_ward <- hclust(famd_dist, method = "ward.D2")
kk <- 11
groups <- cutree(fit_famd_ward, k = kk) 
fct$FG <- groups
colnames(fct)
trait_dat <- fct[c("Species", "n", "max_body_length", "Myelination", "Spawning", "Trophism",
                   "Omnivore", "Carnivore", "Herbivore", "Detritivore", "Feeding_mode",
                   "Current", "Cruise", "Ambush", "FG")]
colnames(trait_dat)[1] <- "species"
colnames(trait_dat)[3] <- "body_size"
#-----------------------------------------------------------------
#plot dendrogram
dend <- fit_famd_ward %>% as.dendrogram 
fit_famd_ward %>% color_branches(k = 11) %>% set("branches_lwd", 2.5)  %>% plot()
setwd(wd_plots)
colors <- c('#8b4513', '#008000', '#4682b4', '#4b0082', 
            '#ff0000', '#ffd700', '#00ff00', '#00ffff', 
            '#0000ff', '#ff1493', '#ffe4b5')
ggsave(plot = plot(dend %>% color_branches(k = 11, col =  colors, groupLabels = T) %>% set("branches_lwd", 2.5), horiz = TRUE), filename = paste0("test.svg"),width = 6, height = 12, dpi = 300)

