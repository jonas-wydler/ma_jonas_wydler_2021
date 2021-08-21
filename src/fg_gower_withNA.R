setwd("C:/Users/Jonas/ma/r/current//")
Sys.setenv(LANG = "en")



#library("mdendro")#this is not good
library('dendextend')
library(stats)
library(ape)
library(ggdendro)
library("tidyverse")
library("vegan")
library("FactoMineR")
library("factoextra")
library("clustMixType")
library("missMDA")
library(ggdendro)
#library(ggdark)

library("tidyverse")
library("raster")
library("reshape2")
library("scales")
library("maps")
library("RColorBrewer")
library("viridis")

require("cluster")
require("FD") 
require("fpc") 

#load datasets

fct <- read.csv("table_funct_traits_copepods_v2.csv", h = T, sep = ";", dec = ",")
fct <- fct[,c(3:20)]



### Count NA
names <- colnames(fct)[c(7,8,9,10,15)] ; names
fct$na_count <- apply(fct[,names], 1, function(x) sum(is.na(x)))
# Drop some NA, discuss with Fabio
fct <- fct[!is.na(fct$max_body_length),]
fct <- fct[fct$na_count < 2,]
# ---------------------------------------------------------------------------------------------------------------------

#not sure why we do this
fct <- fct[]
fct2 <- fct 
traits <- fct
### Convert bolleans to 1/0, keep max_body_length and convert SC to 'ordered' (so gowdis understands that soize class 2 is smaller thansize class 3 for instance)


# test with 1 column
traits$Spawning <- traits$Spawning == "Sac"

traits$Myelination <- as.integer(as.logical(traits$Myelination))
traits$Spawning <- as.integer(as.logical(traits$Spawning)) # 1 means sac
traits$Omnivore <- as.integer(as.logical(traits$Omnivore))
traits$Carnivore <- as.integer(as.logical(traits$Carnivore))
traits$Herbivore <- as.integer(as.logical(traits$Herbivore))
traits$Detritivore <- as.integer(as.logical(traits$Detritivore))
traits$Current <- as.integer(as.logical(traits$Current))
traits$Cruise <- as.integer(as.logical(traits$Cruise))
traits$Ambush <- as.integer(as.logical(traits$Ambush))


# Change size colname
colnames(traits)[7] <- "size"
# Examine structure
str(traits) # ok

# # Add a NA-count column to select species ith enugh information for analyses
# names <- colnames(traits)[c(7,8,9,10,15)] ; names
# traits$na_count <- apply(traits[,names], 1, function(x) sum(is.na(x)))
# (summary(factor(traits$na_count)) / length(unique(traits$Species)))*100

# Compute Gower's distance matrix, with all species having 0 or just 1 NA and then just 0 NA
gow <- gowdis(traits[,c(7,8,9,11,12,13,14,16,17,18)])# maybe we dont need to check for another na
class(gow)

# Draw hierachical dendrogram based on 2 aggregation link: Ward's (like in Benedetti et al. 2016) and UGPMA (Mérigot et al. 2008)
fit_gow_ward <- hclust(gow, method = "ward.D2") 
plot(fit_gow_ward, hang = -1) # display dendogram with names

fit_gow_avg <- hclust(gow, method = "average") 
plot(fit_gow_avg, hang = -1) # display dendogram with names

kk <- 11
groups <- cutree(fit_gow_ward, k = kk)
plot(fit_gow_ward, hang = -1)
rect.hclust(tree = fit_gow_ward, k = kk, which = 1:kk, border = 1:kk, cluster = groups)

### Cut tree to define groups...look at their composition in temrs of species/genus/family and choose those that make most sense
kk <- 14
groups <- cutree(fit_gow_ward, k = kk) 
traits$FG <- NA
traits[traits$na_count < 2,"FG"] <- groups



########
dend <-as.dendrogram(fit_gow_ward)
clust <- cutree(fit_gow_ward, k=kk)
clust.cutree <- dendextend:::cutree(dend, k=kk, order_clusters_as_data = FALSE)
idx <- order(as.numeric(names(clust.cutree)))
clust.cutree <- clust.cutree[idx]
tbl <- table(clust, clust.cutree)
lbls <- apply(tbl,2,which.max)
dend1 <- color_branches(dend, k = kk, groupLabels = lbls)
plot(dend1)

par(bg = 'grey')


d1 <- subset(traits, traits$FG==1)
d2 <- subset(traits, traits$FG==2)
d3 <- subset(traits, traits$FG==3)
d4 <- subset(traits, traits$FG==4)
d5 <- subset(traits, traits$FG==5)
d6 <- subset(traits, traits$FG==6)
d7 <- subset(traits, traits$FG==7)
d8 <- subset(traits, traits$FG==8)
d9 <- subset(traits, traits$FG==9)
d10 <- subset(traits, traits$FG==10)
d11 <- subset(traits, traits$FG==11)
d12 <- subset(traits, traits$FG==12)
d13 <- subset(traits, traits$FG==13)
d14 <- subset(traits, traits$FG==14)

##copy traits to analyze
traits2 <- traits
traits2$FG <- traits$FG
##

####To check summary of traits we have to have factors
traits2$Myelination <- as.factor(as.logical(traits$Myelination))
traits2$Omnivore <- as.factor(as.logical(traits$Omnivore))
traits2$Carnivore <- as.factor(as.logical(traits$Carnivore))
traits2$Herbivore <- as.factor(as.logical(traits$Herbivore))
traits2$Detritivore <- as.factor(as.logical(traits$Detritivore))
traits2$Current <- as.factor(as.logical(traits$Current))
traits2$Cruise <- as.factor(as.logical(traits$Cruise))
traits2$Ambush <- as.factor(as.logical(traits$Ambush))
traits2$Feeding_mode <- as.factor((traits$Feeding_mode))
traits2$Trophism <- as.factor((traits$Trophism))
traits2$Spawning <- as.factor(traits$Spawning)
traits2$FG <- as.factor((traits$FG))


###
d1 <- subset(traits2, traits2$FG==1)
d2 <- subset(traits2, traits2$FG==2)
d3 <- subset(traits2, traits2$FG==3)
d4 <- subset(traits2, traits2$FG==4)
d5 <- subset(traits2, traits2$FG==5)
d6 <- subset(traits2, traits2$FG==6)
d7 <- subset(traits2, traits2$FG==7)
d8 <- subset(traits2, traits2$FG==8)
d9 <- subset(traits2, traits2$FG==9)
d10 <- subset(traits2, traits2$FG==10)
d11 <- subset(traits2, traits2$FG==11)
d12 <- subset(traits2, traits2$FG==12)
d13 <- subset(traits2, traits2$FG==13)
d14 <- subset(traits2, traits2$FG==14)

dtest <- subset(traits,(traits$FG == 7 | traits$FG == 11))
summary(d6)

summary(comparedf(d1,d2))

###


### Define a function/for loop that looks at proportion of traits categories per FG (donut plots) and then retunrs then into a panel plot
# g <- 3
for(g in unique(na.omit(traits$FG)) ) {
  
  subset <- traits[traits$FG == g & !is.na(traits$FG),]
  # Save n species in group g
  nsp <- nrow(subset)
  message(paste("Making plots for FG == ",g," || Nsp = ",nsp, sep = ""))
  # Provide origianl feeding_mode and trophic guild from 'fct2' (makes plotting easier)
  subset$Feeding <- NA
  subset$Trophism <- NA
  subset[,"Feeding"] <- fct[which(fct$Species %in% unique(subset$Species)),"Feeding_mode"]
  subset[,"Trophism"] <- fct[which(fct$Species %in% unique(subset$Species)),"Trophism"]
  # Tally % of traits categories like above
  
  # Plot contrib of Myelination
  tally.myel <- data.frame(subset %>% group_by(Myelination) %>% summarize(Nsp = n(), Psp = (n()/nsp)*100))
  plot.myel <- ggplot(tally.myel, aes(x = 2, y = Psp, fill = factor(Myelination))) +
    geom_bar(stat = "identity", color = "black") + coord_polar(theta = "y", start = 0) +
    scale_fill_brewer(name = "Myelination", palette = "Paired") + theme_void() + xlim(0.5,2.5)
  # Plot contrib of spwaning strategies
  tally.spawn <- data.frame(subset %>% group_by(Spawning) %>% summarize(Nsp = n(), Psp = (n()/nsp)*100))
  plot.spawn <- ggplot(tally.spawn, aes(x = 2, y = Psp, fill = factor(Spawning))) +
    geom_bar(stat = "identity", color = "black") + coord_polar(theta = "y", start = 0) +
    scale_fill_brewer(name = "Spawning mode", palette = "Paired") + theme_void() + xlim(0.5,2.5)
  # Plot contrib of trophism
  tally.tc <- data.frame(subset %>% group_by(Trophism) %>% summarize(Nsp = n(), Psp = (n()/nsp)*100))
  plot.tc <- ggplot(tally.tc, aes(x = 2, y = Psp, fill = factor(Trophism))) +
    geom_bar(stat = "identity", color = "black") + coord_polar(theta = "y", start = 0) +
    scale_fill_brewer(name = "Trophic guild", palette = "Paired") + theme_void() + xlim(0.5,2.5)
  # Plot contrib feeding modes
  tally.feed <- data.frame(subset %>% group_by(Feeding) %>% summarize(Nsp = n(), Psp = (n()/nsp)*100))
  plot.feed <- ggplot(tally.feed, aes(x = 2, y = Psp, fill = factor(Feeding))) +
    geom_bar(stat = "identity", color = "black") + coord_polar(theta = "y", start = 0) +
    scale_fill_brewer(name = "Feeding mode", palette = "Paired") + theme_void() + xlim(0.5,2.5) 
  
  # And a last donut plot for genera or families
  tally.genus <- data.frame(subset %>% group_by(Genus) %>% summarize(Nsp = n(), Psp = (n()/nsp)*100))
  plot.genus <- ggplot(tally.genus, aes(x = 2, y = Psp, fill = factor(Genus))) +
    geom_bar(stat = "identity", color = "black") + coord_polar(theta = "y", start = 0) +
    theme_void() + xlim(0.5,2.5) 
  
  plot.sizes <- ggplot(data = traits[!is.na(traits$FG),], aes(x = factor(FG), y = size, fill = factor(FG))) +
    geom_boxplot(colour = "black") + xlab("Copepod Functional Group") + ylab("Mean body size (mm)") +
    theme_classic() 
  
  require("ggpubr")
  panel <- ggarrange(plot.feed, plot.tc, plot.myel, plot.spawn, plot.genus, plot.sizes, labels = letters[1:6], ncol = 2, nrow = 3)
  # Annotate panel with basic info (group name, Nsp, method)
  panel2 <- annotate_figure(p = panel, top = text_grob(paste("FG",g, "(",nsp," spp) - (Gower+Ward)", sep = ""), color = "black", face = "bold", size = 12))
  ggsave(plot = panel2, filename = paste("plots_contrib_traits_Gower+Ward_FG_k=14#",g,"_withNA.jpg", sep = ""), width = 10, height = 15, dpi = 300)
  
}

### WORKS A LOT BETTER !!! :-)

### Assess clusters validity through C-H index
scores <- function(k) {
  km <- cutree(fit_gow_avg, k)
  ss <- silhouette(km, gow)
  ch <- calinhara(x = gow, clustering = km, cn = k)
  data.frame(kval = k, mean = mean(ss[,3]), sd = sd(ss[,3]), CHI = ch)
} # eo FUN
k <- c(2:15)
avg_silh <- lapply(k,scores) ; ddf_silh <- bind_rows(avg_silh)
quartz()
ggplot() + geom_path(aes(x = factor(kval), y = CHI), linetype = "dashed", data = ddf_silh, group=1) + 
  geom_point(aes(x = factor(kval), y = CHI), pch = 21, colour = "black", fill = "#b8e186", data = ddf_silh) + 
  xlab("Number of clusters (k)") + ylab("Calinski-Harabasz index") +
  theme_classic() 

# 9 or 10 too
kk <-13
groups <- cutree(fit_gow_avg, k = kk)
plot(fit_gow_avg, hang = -1)
rect.hclust(tree = fit_gow_avg, k = kk, which = 1:kk, border = 1:kk, cluster = groups)

### Cut tree to define groups...look at their composition in temrs of species/genus/family and choose those that make most sense

traits$FG <- NA
traits[traits$na_count < 2,"FG"] <- groups



########
dend <-as.dendrogram(fit_gow_avg)
clust <- cutree(fit_gow_avg, k=kk)
clust.cutree <- dendextend:::cutree(dend, k=kk, order_clusters_as_data = FALSE)
idx <- order(as.numeric(names(clust.cutree)))
clust.cutree <- clust.cutree[idx]
tbl <- table(clust, clust.cutree)
lbls <- apply(tbl,2,which.max)
dend1 <- color_branches(dend, k = kk, groupLabels = lbls)
plot(dend1)

par(bg = 'grey')


d1 <- subset(traits, traits$FG==1)
d2 <- subset(traits, traits$FG==2)
d3 <- subset(traits, traits$FG==3)
d4 <- subset(traits, traits$FG==4)
d5 <- subset(traits, traits$FG==5)
d6 <- subset(traits, traits$FG==6)
d7 <- subset(traits, traits$FG==7)
d8 <- subset(traits, traits$FG==8)
d9 <- subset(traits, traits$FG==9)
d10 <- subset(traits, traits$FG==10)
d11 <- subset(traits, traits$FG==11)
d12 <- subset(traits, traits$FG==12)
d13 <- subset(traits, traits$FG==13)
d14 <- subset(traits, traits$FG==14)

##copy traits to analyze
traits2$FG <- traits$FG
##

####To check summary of traits we have to have factors
traits2$Myelination <- as.factor(as.logical(traits$Myelination))
traits2$Omnivore <- as.factor(as.logical(traits$Omnivore))
traits2$Carnivore <- as.factor(as.logical(traits$Carnivore))
traits2$Herbivore <- as.factor(as.logical(traits$Herbivore))
traits2$Detritivore <- as.factor(as.logical(traits$Detritivore))
traits2$Current <- as.factor(as.logical(traits$Current))
traits2$Cruise <- as.factor(as.logical(traits$Cruise))
traits2$Ambush <- as.factor(as.logical(traits$Ambush))
traits2$Feeding_mode <- as.factor((traits$Feeding_mode))
traits2$Trophism <- as.factor(traits$Trophism)
traits2$Spawning <- as.factor(traits$Spawning)
traits2$FG <- as.factor(traits$FG)

###
d1 <- subset(traits2, traits2$FG==1)
d2 <- subset(traits2, traits2$FG==2)
d3 <- subset(traits2, traits2$FG==3)
d4 <- subset(traits2, traits2$FG==4)
d5 <- subset(traits2, traits2$FG==5)
d6 <- subset(traits2, traits2$FG==6)
d7 <- subset(traits2, traits2$FG==7)
d8 <- subset(traits2, traits2$FG==8)
d9 <- subset(traits2, traits2$FG==9)
d10 <- subset(traits2, traits2$FG==10)
d11 <- subset(traits2, traits2$FG==11)
d12 <- subset(traits2, traits2$FG==12)
d13 <- subset(traits2, traits2$FG==13)
d14 <- subset(traits2, traits2$FG==14)

dtest <- subset(traits,(traits$FG == 7 | traits$FG == 11))
summary(d8)


summary(comparedf(d1,d2))

###
# Same as above
for(g in unique(na.omit(traits$FG)) ) {
  
  subset <- traits[traits$FG == g & !is.na(traits$FG),]
  # Save n species in group g
  nsp <- nrow(subset)
  message(paste("Making plots for FG == ",g," || Nsp = ",nsp, sep = ""))
  # Provide origianl feeding_mode and trophic guild from 'fct2' (makes plotting easier)
  subset$Feeding <- NA
  subset$Trophism <- NA
  subset[,"Feeding"] <- fct[which(fct$Species %in% unique(subset$Species)),"Feeding_mode"]
  subset[,"Trophism"] <- fct[which(fct$Species %in% unique(subset$Species)),"Trophism"]
  # Tally % of traits categories like above
  
  # Plot contrib of Myelination
  tally.myel <- data.frame(subset %>% group_by(Myelination) %>% summarize(Nsp = n(), Psp = (n()/nsp)*100))
  plot.myel <- ggplot(tally.myel, aes(x = 2, y = Psp, fill = factor(Myelination))) +
    geom_bar(stat = "identity", color = "black") + coord_polar(theta = "y", start = 0) +
    scale_fill_brewer(name = "Myelination", palette = "Paired") + theme_void() + xlim(0.5,2.5)
  # Plot contrib of spwaning strategies
  tally.spawn <- data.frame(subset %>% group_by(Spawning) %>% summarize(Nsp = n(), Psp = (n()/nsp)*100))
  plot.spawn <- ggplot(tally.spawn, aes(x = 2, y = Psp, fill = factor(Spawning))) +
    geom_bar(stat = "identity", color = "black") + coord_polar(theta = "y", start = 0) +
    scale_fill_brewer(name = "Spawning mode", palette = "Paired") + theme_void() + xlim(0.5,2.5)
  # Plot contrib of trophism
  tally.tc <- data.frame(subset %>% group_by(Trophism) %>% summarize(Nsp = n(), Psp = (n()/nsp)*100))
  plot.tc <- ggplot(tally.tc, aes(x = 2, y = Psp, fill = factor(Trophism))) +
    geom_bar(stat = "identity", color = "black") + coord_polar(theta = "y", start = 0) +
    scale_fill_brewer(name = "Trophic guild", palette = "Paired") + theme_void() + xlim(0.5,2.5)
  # Plot contrib feeding modes
  tally.feed <- data.frame(subset %>% group_by(Feeding) %>% summarize(Nsp = n(), Psp = (n()/nsp)*100))
  plot.feed <- ggplot(tally.feed, aes(x = 2, y = Psp, fill = factor(Feeding))) +
    geom_bar(stat = "identity", color = "black") + coord_polar(theta = "y", start = 0) +
    scale_fill_brewer(name = "Feeding mode", palette = "Paired") + theme_void() + xlim(0.5,2.5) 
  
  # And a last donut plot for genera or families
  tally.genus <- data.frame(subset %>% group_by(Genus) %>% summarize(Nsp = n(), Psp = (n()/nsp)*100))
  plot.genus <- ggplot(tally.genus, aes(x = 2, y = Psp, fill = factor(Genus))) +
    geom_bar(stat = "identity", color = "black") + coord_polar(theta = "y", start = 0) +
    theme_void() + xlim(0.5,2.5) 
  
  plot.sizes <- ggplot(data = traits[!is.na(traits$FG),], aes(x = factor(FG), y = size, fill = factor(FG))) +
    geom_boxplot(colour = "black") + xlab("Copepod Functional Group") + ylab("Mean body size (mm)") +
    theme_classic() 
  
  require("ggpubr")
  panel <- ggarrange(plot.feed, plot.tc, plot.myel, plot.spawn, plot.genus, plot.sizes, labels = letters[1:6], ncol = 2, nrow = 3)
  # Annotate panel with basic info (group name, Nsp, method)
  panel2 <- annotate_figure(p = panel, top = text_grob(paste("FG",g, "(",nsp," spp) - (Gower+Avg)", sep = ""), color = "black", face = "bold", size = 12))
  ggsave(plot = panel2, filename = paste("plots_contrib_traits_Gower+avg_FG_k = 13#",g,"_withNA.jpg", sep = ""), width = 10, height = 15, dpi = 300)
  
}

######dendrogram labels not in the right order, see below
# dend <- as.dendrogram(fit_gow_avg)
# 
# 
# dend <- set(dend,"labels_cex", 0.5)
# dend <- set(dend, "branches_lwd", 2)
# cols_branches <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a')
# dend <- color_branches(dend, k = 10, col = cols_branches, groupLabels = unique(traits$FG))
# plot(dend)
# #gdend <- dend %>% ggdendrogram()
# 
# #gdend + dark_theme_gray() 
# #circlize_dendrogram(dend, dend_track_height = 0.85, group_labels = TRUE)
# 
# 
# a <- as.numeric(as.factor(traits[,10])[order.dendrogram(dend)])
# labels_colors(dend) <- a
# labels(dend) <- as.character(traits[,1])[order.dendrogram(dend)]
# 
# legend(x = 250,y = 7, legend= unique(traits$Trophism[order.dendrogram(dend)]),
#        fill=unique(a), cex=0.8, bty = "n")
# 
# par(bg = 'grey')




###test to use FG as a new variable for another mca
# traits2 <- traits[,-20]
# 
# kk = 10
# groups3 <- cutree(fit_gow_ward, k = kk)
# traits2$FGgowward <- groups3
# groups4 <- cutree(fit_gow_avg, k = kk) 
# traits2$FGgowavg <- groups4
# 
# save(traits2, file = "traits2.RData")


