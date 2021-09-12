#Stability analysis for clustering of ocean regions, Jonas Wydler, 11.05.2021
#updated 12.09.2021
### -----------------------------------------------------------------------

Sys.setenv(LANGUAGE= 'en')

library("tidyverse")
library("cluster")
library("fastcluster")
library("clValid")
library("svglite")
library("reshape2")
### -----------------------------------------------------------------------
#First part; calculate cluster stability indices; takes quite some time
### -----------------------------------------------------------------------
#go to wd where the RData file "coords_traits_biomes_v2.RData" is.
wd_script <- ("/data/manuscript_figure_4/")

#load PCA coords for each grid cell
setwd(wd_script)
load("coords_traits_biomes_v2.RData")

coords <- pca_coords
df_for_clustering <- pca_coords[4:length(pca_coords)]
dist <- dist(df_for_clustering,"euclidean") 

clust <- hclust(dist,"ward.D2")

valid.all <- clValid(obj = df_for_clustering, nClust = c(3:10), clMethods = c("hierarchical"),
                     method = "ward", validation = "internal", metric = "euclidean", 
                    verbose = T, maxitems = nrow(df_for_clustering) )

save(valid.all, file = paste("clValid_indices_traits_Biomes_v2", sep = ""))

### -----------------------------------------------------------------------
#second part; plot cluster stability indices
### -----------------------------------------------------------------------

load("clValid_indices_traits_Biomes_v2")
#test <- as.data.frame(valid.all@measNames)

slotNames(valid.all)
## view results
valid.all
summary(valid.all)
optimalScores(valid.all)
plot(valid.all)

#change to numeric for ploting
df_indices <- as.data.frame(valid.all@measures)
dat <- data.frame(n = c(3:10))
dat$Connectivity <- as.numeric(df_indices[1,1:8])
dat$Dunn <- as.numeric(df_indices[2,1:8])
dat$Silhouette <- as.numeric(df_indices[3,1:8])

dat_meltd <- melt(dat, id.vars = c("n"))

dat_Con <- subset(dat_meltd, variable == "Connectivity")
dat_Dunn <- subset(dat_meltd, variable == "Dunn")
dat_Sil <- subset(dat_meltd, variable == "Silhouette")

#plots
g_Con <- ggplot(dat_Con, aes(x=n, y=value, fill = variable)) +
  geom_point(size=2, shape=19) + geom_line() + theme_minimal()

g_dunn <- ggplot(dat_Dunn, aes(x=n, y=value, fill = variable)) +
  geom_point(size=2, shape=19) + geom_line() + theme_minimal()

g_sil <- ggplot(dat_Sil, aes(x=n, y=value, fill = variable)) +
  geom_point(size=2, shape=19) + geom_line() + theme_minimal()

setwd("C:/Users/Jonas/polybox/arbeit_upgroup/trait_biome_map/")
ggsave(g_Con, filename = "conn.svg", width = 12, height = 7, dpi = 300)

setwd("C:/Users/Jonas/polybox/arbeit_upgroup/trait_biome_map/")
ggsave(g_dunn, filename = "dunn.svg", width = 12, height = 7, dpi = 300)

setwd("C:/Users/Jonas/polybox/arbeit_upgroup/trait_biome_map/")
ggsave(g_sil, filename = "sil.svg", width = 12, height = 7, dpi = 300)
