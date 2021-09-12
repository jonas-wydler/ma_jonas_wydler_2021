#-----------------------------------------------------------
#Script to invesigate sampling effort and plot spatial and temporal patterns of sampling effort 
#Jonas Wydler 27.07.2021. Updated 13.09.2021
#-----------------------------------------------------------

library("tidyverse")
library("RColorBrewer") # nice color palettes (discrete and continuous)
library("viridis") # even better palettes for continuous data like SST etc.
library("cmocean") # very nice palettes as well
library("svglite")
library("ggpubr")
### Load the results 

world2 <- map_data("world2") # world coastline for maps 

setwd("/data/data_species_v3.2v5.1_thinned2/")
#setwd("")
res <- lapply(dir(), function(f) {
  dat <- read.table(f, sep = "\t",  h = T)
  return(dat)
}
) # eo lapply
ddf <- bind_rows(res)
# Tally per species etc.
tally <- data.frame(ddf %>% group_by(species) %>% summarize(n = n() ) )

tally <- tally[order(tally$n, decreasing = T),] 
tally
# Check new threshold perhaps
nrow(tally[tally$n >= 30,]) # 357/386 species --> 92%
nrow(tally[tally$n >= 40,]) # 333/386 --> 86%
nrow(tally[tally$n >= 50,]) # 304/386 --> 79% #why 304 and not 343 like other data

### Make maps of sampling effort and raster of effort per month, like previously
ddf$id <- factor( paste(ddf$x, ddf$y, sep = "_") )
ddf <- ddf[order(ddf$id),]
# ddf[ddf$id == "146_-39",]
effort <- data.frame(ddf %>% group_by(id) %>% summarize(x = unique(x), y = unique(y), n = n(), rich = length(unique(species))) ) # eo ddf
effort$logn <- log(effort$n)
dim(effort)
summary(effort)

### Make bins for nicer color palette
effort$logn_bin <- factor(cut_interval(effort$logn,9))
levels(effort$logn_bin)
levels <- str_replace_all(levels(effort$logn_bin), ",", "-")
levels <- gsub("\\[|\\]", "", levels)
levels <- gsub("\\(|\\)", "", levels)
levels
levels(effort$logn_bin) <- levels
# Rotate x axis
effort$x2 <- effort$x
effort[effort$x < 0 ,"x2"] <- (effort[effort$x < 0 ,"x"]) + 360

map <- ggplot() + geom_raster(aes(x = x2, y = y, fill = factor(logn_bin)), data = effort) + 
  scale_fill_manual(name = "Sampling effort\n(ln)", values = cmocean('thermal')(10) ) +
  geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
  coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                                        labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
  scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
                     labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
  theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
        panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
pcomplete <- plot(map)
#setwd("C:/Users/Jonas/Documents/ETH/master/ma/sdm/data_species_v3.2v5.1_thinned/")
#ggsave(plot = map, filename = "map_effort_logn_thinned.jpg", dpi = 300, width = 7, height = 4)

### And examine sampling effort per month
ddf$ybin_5d <- ceiling(ddf$y/5)*5
effort3 <- data.frame(ddf %>% group_by(month,ybin_5d) %>% summarize(n = n(), rich = length(unique(species))) )
summary(effort3)
effort3$logn <- log(effort3$n)
effort3$logn_bin <- factor(cut_interval(effort3$logn,9))
levels(effort3$logn_bin)
levels <- str_replace_all(levels(effort3$logn_bin), ",", "-")
levels <- gsub("\\[|\\]", "", levels)
levels <- gsub("\\(|\\)", "", levels)
levels
levels(effort3$logn_bin) <- levels

plot <- ggplot() + geom_raster(aes(x = factor(ybin_5d), y = factor(month), fill = factor(logn_bin)), data = effort3) + 
  scale_fill_manual(name = "Sampling effort\n(ln)", values = cmocean('thermal')(10)) +
  ylab("Month") + xlab("Latitude (5°)") + theme_light() + 
  theme(axis.text.x = element_text(angle=90, hjust=1, size = 5)) 
#
#ggsave(plot = plot, filename = "plot_effort_mon_lat_logn_thinned.jpg", dpi = 300, width = 8, height = 4)
### Ok, same patterns as before

#=====================================================================================

#merge fg dataset to occurences
setwd("C:/Users/Jonas/ma/sdm/")
load("species_FG")
colnames(species_FG) <- c("species","FG_f11")
ddf_fg <- base::merge(ddf,species_FG, by  = "species")



### Make maps of sampling effort and raster of effort per month, for fgs
ddf_fg$id <- factor(paste(ddf_fg$x, ddf_fg$y, sep = "_") )
ddf_fg <- ddf_fg[order(ddf_fg$id),]

plots = list()
for (i in 1:length(unique(ddf_fg$FG_f11))){
  #do effort calc. per group
  ddf_fg_x <- subset(ddf_fg, FG_f11 == i)
  
  # ddf[ddf$id == "146_-39",]
  effort_fg_x <- data.frame(ddf_fg_x %>% group_by(id) %>% summarize(x = unique(x), y = unique(y), n = n(), rich = length(unique(species))) ) # eo ddf
  effort_fg_x$logn <- log(effort_fg_x$n)

  
  ### Make bins for nicer color palette
  effort_fg_x$logn_bin <- factor(cut_interval(effort_fg_x$logn,9))
  levels(effort_fg_x$logn_bin)
  levels <- str_replace_all(levels(effort_fg_x$logn_bin), ",", "-")
  levels <- gsub("\\[|\\]", "", levels)
  levels <- gsub("\\(|\\)", "", levels)
  levels
  levels(effort_fg_x$logn_bin) <- levels
  # Rotate x axis
  effort_fg_x$x2 <- effort_fg_x$x
  effort_fg_x[effort_fg_x$x < 0 ,"x2"] <- (effort_fg_x[effort_fg_x$x < 0 ,"x"]) + 360
  
  
  map_fg_x <- ggplot() + geom_raster(aes(x = x2, y = y, fill = factor(logn_bin)), data = effort_fg_x) + 
    scale_fill_manual(name = "Sampling effort\n(log)", values = cmocean('thermal')(10) ) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(breaks = c(0,60,120,180,240,300,360),
                                          labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(breaks = c(-90,-60,-30,0,30,60,90),
                       labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
          panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
  
  #plots[[i]] <- map_fg_x
  #ggsave(map_fg_x,filename = paste0("occurrence_data_FGs",i,".png"),width=10, height=7, dpi=600)
  
}
lol <- append(plots,list(map))

panel <- ggarrange(plotlist = lol, labels = c(1:11,"a"), ncol = 3, nrow = 4, common.legend = T)

# Annotate panel with basic info (group name, Nsp, method)
#ggsave(plot(panel))

plot(panel)
p1 <- plot(panel)
ggsave(plot(panel),filename = "occurrence_data_FGs.svg",width=10, height=7, dpi=600)

#=====================================================================================
### And examine sampling effort per month



# ddf[ddf$id == "146_-39",]
ddf_fg$ybin_5d <- ceiling(ddf_fg$y/5)*5
effort_fg_heat <- data.frame(ddf_fg %>% group_by(FG_f11,ybin_5d) %>% summarize(n = n(), rich = length(unique(species))) )
summary(effort_fg_heat)

effort_fg_heat$logn <- log(effort_fg_heat$n)
effort_fg_heat$logn_bin <- factor(cut_interval(effort_fg_heat$logn,9))
levels(effort_fg_heat$logn_bin)
levels <- str_replace_all(levels(effort_fg_heat$logn_bin), ",", "-")
levels <- gsub("\\[|\\]", "", levels)
levels <- gsub("\\(|\\)", "", levels)
levels
levels(effort_fg_heat$logn_bin) <- levels

plot_fg_heat <- ggplot() + geom_raster(aes(x = factor(ybin_5d), y = factor(FG_f11), fill = factor(logn_bin)), data = effort_fg_heat) + 
  scale_fill_manual(name = "Sampling effort\n(log)", values = cmocean('thermal')(10)) +
  ylab("FG") + xlab("Latitude (5°)") + theme_light() + 
  theme(axis.text.x = element_text(angle=90, hjust=1, size = 5)) 
#
#ggsave(plot = plot_fg_heat, filename = "Heatmap_FGs.jpg", dpi = 300, width = 8, height = 4)
#ggsave(plot = plot, filename = "plot_effort_mon_lat_logn_thinned.jpg", dpi = 300, width = 8, height = 4)
### Ok, same patterns as before

#=====================================================================================

#look at species richness distribution for > 5, 10, 25 and 50 species observed
ddf$id <- factor( paste(ddf$x, ddf$y, sep = "_") )
ddf <- ddf[order(ddf$id),]
# ddf[ddf$id == "146_-39",]
effort <- data.frame(ddf %>% group_by(id) %>% summarize(x = unique(x), y = unique(y), n = n(), rich = length(unique(species))) ) # eo ddf
effort$logn <- log10(effort$n)

effort5 <- subset(effort, rich > 5)
effort10 <- subset(effort, rich > 10)
effort25 <- subset(effort, rich > 25)
effort50 <- subset(effort, rich > 50)

plotsrich <- list()
vec <- c(5,10,25,50)
for (i in 1:length(vec)){
  #do effort calc. per group
  effortRICH <- subset(effort, rich > vec[i])
  effortRICH$logn <- log10(effortRICH$n)
  
  ### Make bins for nicer color palette
  effortRICH$logn_bin <- factor(cut_interval(effortRICH$logn,9))
  levels(effortRICH$logn_bin)
  levels <- str_replace_all(levels(effortRICH$logn_bin), ",", "-")
  levels <- gsub("\\[|\\]", "", levels)
  levels <- gsub("\\(|\\)", "", levels)
  levels
  levels(effortRICH$logn_bin) <- levels
  # Rotate x axis
  effortRICH$x2 <- effortRICH$x
  effortRICH[effortRICH$x < 0 ,"x2"] <- (effortRICH[effortRICH$x < 0 ,"x"]) + 360
  
  
  map_rich <- ggplot() + geom_raster(aes(x = x2, y = y, fill = factor(logn_bin)), data = effortRICH) + 
    scale_fill_manual(name = "Sampling effort\n(log 10)", values = cmocean('thermal')(10) ) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
                                          labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
    scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
                       labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
          panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
  
  plotsrich[[i]] <- map_rich
  
}


panel <- ggarrange(plotlist =plotsrich, labels = c(">5",">10",">25",">50"), ncol = 2, nrow = 2)
plot(panel)
