#Jonas Wydler;  dendrogram for paper; 23.05.2021
#Careful about the FG names, in order to have a sequential order on the famd+ward plot
#(figure 1 in the manuscript) we used the order assinged here and changed it everywhere else in the text
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
library("FD")
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
###plot dendrogram famd ward
dend <- fit_famd_ward %>% as.dendrogram 

fit_famd_ward %>% color_branches(k = 11) %>% set("branches_lwd", 2.5)  %>% plot()
setwd(wd_plots)
colors <- c('#8b4513', '#008000', '#4682b4', '#4b0082', 
            '#ff0000', '#ffd700', '#00ff00', '#00ffff', 
            '#0000ff', '#ff1493', '#ffe4b5')
ggsave(plot = plot(dend %>% color_branches(k = 11, col =  colors, groupLabels = T) %>% set("branches_lwd", 2.5), horiz = TRUE), filename = paste0("dend_famd_ward.svg"),width = 6, height = 12, dpi = 300)

plot = plot(dend %>% color_branches(k = 11, col =  colors, groupLabels = T) %>% set("branches_lwd", 2.5), horiz = TRUE)

table_traits <- as.data.frame(trait_dat %>% group_by(cell_id) %>% summarize())

table_traits_subset <- subset(trait_dat, FG == 6)
table(table_traits_subset$Feeding_mode)  
length((table_traits_subset$Feeding_mode))

#-----------------------------------------------------------------
###plot dendrogram famd average
fit_famd_avg <- hclust(famd_dist, method = "average")
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
###plot dendrogram famd ward
#careful not directly comparable as the groups do not necessarily align with dendrogramgs derived from other methods
dend2 <- fit_famd_avg %>% as.dendrogram 

fit_famd_avg %>% color_branches(k = 11) %>% set("branches_lwd", 2.5)  %>% plot()
setwd(wd_plots)
colors <- c('#8b4513', '#008000', '#4682b4', '#4b0082', 
            '#ff0000', '#ffd700', '#00ff00', '#00ffff', 
            '#0000ff', '#ff1493', '#ffe4b5')

ggsave(plot = plot(dend2  %>% set("branches_lwd", 2.5), horiz = TRUE), filename = paste0("dend_famd_avg.png"),width = 6, height = 12, dpi = 300)
plot = plot(dend2  %>% set("branches_lwd", 2.5), horiz = TRUE)

#-----------------------------------------------------------------
###plot dendrograms for gower distance 
#careful not directly comparable as the groups do not necessarily align with dendrogramgs derived from other methods

# Compute Gower's distance matrix, with all species having 0 or just 1 NA and then just 0 NA
gow <- gowdis(fct[,c(7:9,11:14,16:18)])# maybe we dont need to check for another na
fit_gow_ward <- hclust(gow, method = "ward.D2")
dend3 <- fit_gow_ward %>% as.dendrogram
plot = plot(dend3  %>% set("branches_lwd", 2.5), horiz = TRUE)
ggsave(plot = plot(dend3  %>% set("branches_lwd", 2.5), horiz = TRUE), filename = paste0("dend_gow_ward.png"),width = 6, height = 12, dpi = 300)


# Compute Gower's distance matrix, with all species having 0 or just 1 NA and then just 0 NA
gow <- gowdis(fct[,c(7:9,11:14,16:18)])# maybe we dont need to check for another na
fit_gow_avg <- hclust(gow, method = "average")
dend4 <- fit_gow_avg %>% as.dendrogram
plot = plot(dend4  %>% set("branches_lwd", 2.5), horiz = TRUE)
ggsave(plot = plot(dend4  %>% set("branches_lwd", 2.5), horiz = TRUE), filename = paste0("dend_gow_avg.png"),width = 6, height = 12, dpi = 300)


