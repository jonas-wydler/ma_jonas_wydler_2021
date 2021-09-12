#Jonas Wydler, 11.5.2021, The purpose of this script is to define biomes for CWM Traits
#updated 12.09.2021
### ---------------------------------------------------------------------
Sys.setenv(LANGUAGE= 'en')
#library("raster")
#library("sp")
#library("rgdal")
#library("viridis")
library("FactoMineR")
library("tidyverse") 
#library("cluster")
library("missMDA")
library("fastcluster")
library("svglite") 
#library("ggrepel")
### ---------------------------------------------------------------------
wd_script2 <- ("/data/") #directory where you want to save plots
### ---------------------------------------------------------------------
world2 <- map_data("world2") # world coastline for maps 
world1 <- map_data("world")
### ---------------------------------------------------------------------
setwd(wd_script2)
load("data_total_everything.RData") #load dataset containing the cwm traits 
colnames(data_total_everything)
dat_traits <- data_total_everything[c("cell_id", "x", "y", "body_size", "sac_prop", "myel_prop", "ambush_prop",       
                                      "cruise_prop", "current_prop", "current_ambush_prop", "current_cruise_prop",
                                      "carni_prop", "herbi_prop", "detri_prop", "omni_prop")]

dat_traits <- na.omit(dat_traits)
colnames(dat_traits[4:length(dat_traits)])
estim_ncp(X = dat_traits[4:length(dat_traits)])
pca1 <- PCA(X = dat_traits[2:length(dat_traits)], ncp = 11, graph = T, scale = T, quanti.sup = c(1,2))
pca_coords <- as.data.frame(pca1$ind$coord)
pca_coords <- cbind(dat_traits[1:3],pca_coords[c(1:5)])


### -----------------------------------------------------------------------
eig <- data.frame(perc = pca1$eig[,"percentage of variance"], nb = c(1:nrow(pca1$eig))) # eig
pc1 <- paste0("PC1 (",floor(eig$perc[1]*100)/100,"%)")
pc2 <- paste0("PC2 (",floor(eig$perc[2]*100)/100,"%)")
pc3 <- paste0("PC3 (",floor(eig$perc[3]*100)/100,"%)")
pc4 <- paste0("PC4 (",floor(eig$perc[4]*100)/100,"%)")
pc5 <- paste0("PC5 (",floor(eig$perc[5]*100)/100,"%)")

### -----------------------------------------------------------------------

df_for_clustering <- pca_coords[4:length(pca_coords)]
dist <- dist(df_for_clustering,"euclidean") 
#start_time <- Sys.time()
fit_ward <-  fastcluster::hclust(dist, method = "ward.D2")
#end_time <- Sys.time()
#end_time - start_time

reg <- c(3:8)
i = 5
res <- lapply(reg, function(i){
  km <- cutree(fit_ward, i)#regions = km
  message(paste0(i))
  
  df <- cbind(pca_coords, km)
  df$regions <- as.factor(df$km)
  
  
  g1 <- ggplot() + geom_tile(aes(x = x, y = y, fill = regions), data = df) +
    # play around with palettes. You can also use scale_fill_distiller() and
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                             panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
    theme(legend.text=element_text(size=6), legend.title = element_text(size=8)) +
    scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    scale_fill_brewer(palette = "Paired") + coord_map("mollweide", orientation = c(90,-180,0))
  #setwd(wd_script2)
  #ggsave(plot = g1, filename = paste0("Ocean_",i,"_biomes_traits_moll_v2.png"),dpi = 300, width = 7, height = 4)
  
  return(g1)
}
)
#save pca_coords for cluster stability analysis
setwd(wd_script2)
save(pca_coords, file = "coords_traits_biomes_v2.RData")

#---------------------------------------------------------
#part below is to analyse and plot the pca that was used in the defition of the regions

###########
#function for nicer pca plots
augment.PCA <- function(x, dims = c(1:7), which="col") {
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


pcad1 <- ggplot(pcad) +
  coord_fixed() + scale_x_continuous(breaks=0) + scale_y_continuous(breaks=0) +
  annotate(geom="path", colour="black", x=cos(seq(0, 2*pi, length.out=100)), y=sin(seq(0, 2*pi, length.out=100)), colour = "grey30") +
  geom_segment(aes(x=0, xend = Dim.1, y=0, yend = Dim.2), arrow=arrow(angle=20, length=unit(0.01, "npc"))) +
  geom_text_repel(aes(x=Dim.1, y=Dim.2, label=var), 
                  #data=filter(pcad, (Dim.1^2+Dim.2^2) > 0.2^2), segment.alpha=0.5) +
  ) + xlab(pc1) + ylab(pc2) + theme_bw()

setwd(wd_script2)
ggsave(plot = pcad1, filename = paste("PCA_trait_regions_1_2.svg", sep = ""), width = 9, height = 9, dpi = 300)

