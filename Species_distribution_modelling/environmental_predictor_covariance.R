#Rscript Master thesis Jonas Wydler, 03.09.2021
# - Investigate covariance of environmental predictors and create corresponding plots
#-----------------------------------
library("tidyverse")
library("parallel") 
library("ggpubr")
library("reshape2") 
#-----------------------------------
### Draw the correlation coeff heatmaps for each net
get_lower_tri <- function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
reorder_cormat <- function(cormat){
  # Utiliser la corrélation entre les variables
  # comme mésure de distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <- cormat[hc$order, hc$order]
}

setwd(paste0("wd_data/species_background_data/")) # go to dir with background data
files <- dir()[grep("backgrounddata",dir())] #; files
#f <- files[100] # for testing
res <- mclapply(files, function(f) {
  
  # Message
  data <- read.table(f, sep = ";", h = T)
  message(paste(unique(data[data$obs == 1,'species']), sep = ""))
  names <- colnames(data)[c(10,11,12,14,18,21,22,23,24,26,28,29,30,31)] # indices of variables to keep
  mydata <- na.omit(data[,names]) #  head(mydata)
  cormat <- round(cor(mydata, method = "spearman"),2)
  upper_tri <- get_upper_tri(cormat)
  melted_cormat <- melt(upper_tri, na.rm = TRUE)
  melted_cormat <- melted_cormat[!(melted_cormat$Var1 == melted_cormat$Var2),]
  
  melted_cormat$species <- unique(data[data$obs == 1,"species"])
  melted_cormat$group <- unique(data[data$obs == 1,"FG_f11"]) #group2 would your FG, let’s check correlation structure overall and then per FG
            return(melted_cormat)
    
    }, mc.cores = 1 # for running on kryo
    
)
#eo mclapply

table <- bind_rows(res) ; rm(res) ; gc()
dim(table) ; head(table)
### Compute mean or median correlation
table$id <- paste(table$Var1, table$Var2, sep = "_")

ddf <- data.frame(table %>% group_by(id) %>% summarize(V1 = unique(Var1), V2 = unique(Var2), mean = round(mean(value),2), median = round(median(value),2)) )
#ddfFG <- data.frame(table %>% group_by(id) %>% summarize(V1 = unique(Var1), V2 = unique(Var2), mean = round(mean(value),2), median = round(median(value),2)) )

head(ddf) ; summary(ddf)

# ggplot(ddf, aes(factor(V2), factor(V1), fill = mean)) + geom_tile(color = "white") +
#     scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white", midpoint = 0, limits = c(-1,1), name = "Mean\nRho") +
#     xlab("") + ylab("") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1),
#         axis.text.y = element_text(vjust = 1, size = 10, hjust = 1)) + coord_fixed() +
#     geom_text(aes(factor(V2), factor(V1), label = mean), color = "black", size = 3) 




ggheatmap <- ggplot(ddf, aes(factor(V2), factor(V1), fill = mean))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Mean\nRho (spearman)") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+coord_fixed()

print(ggheatmap)

ggheatmap + 
  geom_text(aes(factor(V2), factor(V1), label = mean), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))



# 
# ggplot(ddf, aes(factor(V2), factor(V1), fill = median)) + geom_tile(color = "white") +
#     scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white", midpoint = 0, limits = c(-1,1), name = "Median\nRho") +
#     xlab("") + ylab("") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1),
#         axis.text.y = element_text(vjust = 1, size = 10, hjust = 1)) + coord_fixed() +
#     geom_text(aes(factor(V2), factor(V1), label = median), color = "black", size = 3)

ggheatmap <- ggplot(ddf, aes(factor(V2), factor(V1), fill = median))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Median\nRho (spearman)") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+coord_fixed()

print(ggheatmap)

ggheatmap + 
  geom_text(aes(factor(V2), factor(V1), label = median), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))


    
### Also facet_wrap theses plots across FG 
plotscov <- list()
for(i in 1:11){
  subtable <- subset(table,group == i)
  subtable$id <- paste(subtable$Var1, subtable$Var2, sep = "")
  ddf <- data.frame(subtable %>% group_by(id) %>% summarize(V1 = unique(Var1), V2 = unique(Var2), mean = round(mean(value),2), median = round(median(value),2)) )
  
  ggheatmap <- ggplot(ddf, aes(factor(V2), factor(V1), fill = median))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Mean\nRho (spearman)") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 8, hjust = 1))+coord_fixed()
  
  
  
  ggheatmap <- ggheatmap + 
    geom_text(aes(factor(V2), factor(V1), label = median), color = "black", size = 1.5) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.6, 0.7),
      legend.direction = "horizontal",
     
      
      )+ theme(legend.position="none")
 
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                 title.position = "top", title.hjust = 0.5))
  plotscov[[i]] <- ggheatmap
  #
}
panel <- ggarrange(plotlist =plotscov, labels = c(1:11), ncol = 4, nrow = 3)
setwd("YOUR DIR for plots")
ggsave(plot = panel, filename = paste("name.png"),dpi = 900, width = 16, height = 9)

