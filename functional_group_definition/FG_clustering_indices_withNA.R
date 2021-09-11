# Rscript, Jonas Wydler, 16.10.2020
# - calculating and plotting clustering stability indices
#updated 08.09.2021
# --------------------------------------------------------------------------------------------------------------------
setwd("wd_data/species_background_data/")
Sys.setenv(LANG = "en")

#Packages
library("clValid")
library("tidyverse")
library("directlabels")
library("grid") 
library("fpc")

#load FAMD/Gower space data
load("famd_withNA.RData")
load("gow_withNA.RData")
famd_dist <- dist(famd_sp, method = "euclidean")

# --------------------------------------------------------------------------------------------------------------------
#clustering for k groups and calculating clustering stability indices
k <- c(4:13)
indices <- function(k, distance, link){
  
  dist = distance
  clusterObj = hclust(distance, method = link)
  nc = k
  cluster <- cutree(clusterObj,nc)
  
  con = connectivity(distance, cluster)
  dun = dunn(distance, cluster)
  ss = silhouette(x = cluster, dist = distance)
  ch = calinhara(distance, cluster, k)
  
  data.frame(nCl = k, MEAN = mean(ss[,3]), CON = con, CHI = ch, DUN = dun)
}


index_famd_ward_temp <- lapply(k, distance = famd_dist,link  = "ward.D2", indices)
index_famd_ward <- bind_rows(index_famd_ward_temp)
index_famd_ward_long <- gather(index_famd_ward, key = index_name, value = index_val, MEAN, CON, CHI, DUN)
index_famd_ward_long$type <- "famd"
index_famd_ward_long$link <- "ward"
index_famd_ward_long$desc <- "Famd & Ward"

index_famd_average_temp <- lapply(k, distance = famd_dist,link  = "average", indices)
index_famd_average <- bind_rows(index_famd_average_temp)
index_famd_average_long <- gather(index_famd_average, key = index_name, value = index_val, MEAN, CON, CHI, DUN)
index_famd_average_long$type <- "famd"
index_famd_average_long$link <- "average"
index_famd_average_long$desc <- "Famd & Average"

index_gow_ward_temp <- lapply(k, distance = gow,link  = "ward.D2", indices)
index_gow_ward <- bind_rows(index_gow_ward_temp)
index_gow_ward_long <- gather(index_gow_ward, key = index_name, value = index_val, MEAN, CON, CHI, DUN)
index_gow_ward_long$type <- "gow"
index_gow_ward_long$link <- "ward"
index_gow_ward_long$desc <- "Gower & Ward"

index_gow_average_temp <- lapply(k, distance = gow,link  = "average", indices)
index_gow_average <- bind_rows(index_gow_average_temp)
index_gow_average_long <- gather(index_gow_average, key = index_name, value = index_val, MEAN, CON, CHI, DUN)
index_gow_average_long$type <- "gow"
index_gow_average_long$link <- "average"
index_gow_average_long$desc <- "Gower & Average"

ind <- bind_rows(index_famd_ward_long,index_famd_average_long,index_gow_ward_long,index_gow_average_long)

ind_chi <- ind[ind$index_name == "CHI",]
chi_Mean_temp <-  aggregate(ind_chi[, 3], list(ind_chi$nCl), mean)
chi_Mean_temp$desc <- "Mean"
colnames(chi_Mean_temp) <- c("nCl","index_val","desc")
ind_chi <- dplyr::bind_rows(ind_chi, chi_Mean_temp)

ind_con <- ind[ind$index_name == "CON",]
con_Mean_temp <-  aggregate(ind_con[, 3], list(ind_con$nCl), mean)
con_Mean_temp$desc <- "Mean"
colnames(con_Mean_temp) <- c("nCl","index_val","desc")
ind_con <- dplyr::bind_rows(ind_con, con_Mean_temp)

ind_dun <- ind[ind$index_name == "DUN",]
dun_Mean_temp <-  aggregate(ind_dun[, 3], list(ind_dun$nCl), mean)
dun_Mean_temp$desc <- "Mean"
colnames(dun_Mean_temp) <- c("nCl","index_val","desc")
ind_dun <- dplyr::bind_rows(ind_dun, dun_Mean_temp)

ind_mean <- ind[ind$index_name == "MEAN",]
mean_Mean_temp <-  aggregate(ind_mean[, 3], list(ind_mean$nCl), mean)
mean_Mean_temp$desc <- "Mean"
colnames(mean_Mean_temp) <- c("nCl","index_val","desc")
ind_mean <- dplyr::bind_rows(ind_mean, mean_Mean_temp)
# --------------------------------------------------------------------------------------------------------------------
#indiviudal plots
colors <- c("#CA0020", "#F4A582", "#92C5DE", "#0571B0","#000000") 


g1v2 <- ggplot(data = ind_chi, aes(x = nCl, y = index_val)) +
  geom_path(data = ind_chi[ind_chi$desc != "Mean",], aes(x = nCl, y = index_val, colour = desc), linetype = "dashed", size = 1.2)+
  geom_point(data = ind_chi[ind_chi$desc != "Mean",], aes(x = nCl, y = index_val, fill = desc, colour = desc), size = 4, pch = 21) + 
  geom_line(aes(x = nCl, y = index_val), data = ind_chi[ind_chi$desc == "Mean",], size= 2, color = "black", linetype = "solid") + 
  theme_classic() + theme(legend.position = "none") + xlab("Number of clusters (k)") + 
  ylab("Calinski-Harabasz index") +
  scale_x_continuous(limits = c(4, 15)) +
  #scale_y_continuous(limits = c(90, 430)) +
  geom_dl(aes(label = desc, colour = desc), method=list(dl.trans(x = x + .3), "last.bumpup"), data = ind_chi)  + 
  scale_colour_manual(values = colors) + scale_fill_manual(values = colors) +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))


g2v2 <- ggplot(data = ind_con, aes(x = nCl, y = index_val))  +
  geom_path(data = ind_con[ind_con$desc != "Mean",], aes(x = nCl, y = index_val, colour = desc), linetype = "dashed", size = 1.2)+
  geom_point(data = ind_con[ind_con$desc != "Mean",], aes(x = nCl, y = index_val, fill = desc, colour = desc), size = 4, pch = 21) + 
  geom_line(aes(x = nCl, y = index_val), data = ind_con[ind_con$desc == "Mean",], size= 2, color = "black", linetype = "solid") + 
  theme_classic() + theme(legend.position = "none") + xlab("Number of clusters (k)") + 
  ylab("Connectivity") +
  scale_x_continuous(limits = c(4, 15)) +
  #scale_y_continuous(limits = c(5, 33)) +
  geom_dl(aes(label = desc, colour = desc), method=list(dl.trans(x = x + .3), "last.bumpup"), data = ind_con) + 
  scale_colour_manual(values = colors) + scale_fill_manual(values = colors) +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 

g3v2 <- ggplot(data = ind_dun, aes(x = nCl, y = index_val))  +
  geom_path(data = ind_dun[ind_dun$desc != "Mean",], aes(x = nCl, y = index_val, colour = desc), linetype = "dashed", size = 1.2)+
  geom_point(data = ind_dun[ind_dun$desc != "Mean",], aes(x = nCl, y = index_val, fill = desc, colour = desc), size = 4, pch = 21) + 
  geom_line(aes(x = nCl, y = index_val), data = ind_dun[ind_dun$desc == "Mean",], size= 2, color = "black", linetype = "solid") + 
  theme_classic() + theme(legend.position = "none") + xlab("Number of clusters (k)") + 
  ylab("Dunn index") +
  scale_x_continuous(limits = c(4, 15)) +
  #scale_y_continuous(limits = c(0,0.35)) +
  geom_dl(aes(label = desc, colour = desc), method=list(dl.trans(x = x + .3), "last.bumpup"), data = ind_dun) + 
  scale_colour_manual(values = colors) + scale_fill_manual(values = colors)+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

g4v2 <- ggplot(data = ind_mean, aes(x = nCl, y = index_val))  +
  geom_path(data = ind_mean[ind_mean$desc != "Mean",], aes(x = nCl, y = index_val, colour = desc), linetype = "dashed", size = 1.2)+
  geom_point(data = ind_mean[ind_mean$desc != "Mean",], aes(x = nCl, y = index_val, fill = desc, colour = desc), size = 4, pch = 21) + 
  geom_line(aes(x = nCl, y = index_val), data = ind_mean[ind_mean$desc == "Mean",], size= 2, color = "black", linetype = "solid") + 
  theme_classic() + theme(legend.position = "none") + xlab("Number of clusters (k)") + 
  ylab("Avg. Silhouette width") +
  scale_x_continuous(limits = c(4, 15)) +
  #scale_y_continuous(limits = c(0.4,0.8)) +
  geom_dl(aes(label = desc, colour = desc), method = list(dl.trans(x = x + .3), "last.bumpup"), data = ind_mean) + 
  scale_colour_manual(values = colors) + scale_fill_manual(values = colors) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

#save panel plot
grid.arrange(cbind(rbind(ggplotGrob(g1v2), ggplotGrob(g3v2)), rbind(ggplotGrob(g2v2), ggplotGrob(g4v2))), top=textGrob("Indices (Including Data with one NA, n = 343)",gp=gpar(fontsize=14,font=3), just = "top"))
ggsave(file = "Combined_Indices_withoutNa.png")



