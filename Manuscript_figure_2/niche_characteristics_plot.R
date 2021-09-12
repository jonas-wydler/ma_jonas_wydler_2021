#Environmental niche characteristics plots for manuscript, Jonas Wydler, 01.06.2021
#updated 12.09.2021
### -----------------------------------------------------------------------
#Packages
library("tidyverse")
library("FactoMineR")
library("parallel")
library("reshape2")
library("svglite")

Sys.setenv(LANGUAGE= 'en')
### -----------------------------------------------------------------------
wd_background_data <- ("/data/species_background_data/") #directory load background data from
wd_niche_data <- ("/data/niche_traits/") #directory to load niche characteristics data for furhter use in PCA
wd_plots <- ("/data/") #directory save plots to
### -----------------------------------------------------------------------
setwd(wd_niche_data)
#Load data
files <- dir()[grep("table_niche_traits_GAM_group",dir())] # files
res <- mclapply(files, function(f) {
  t <- get(load(f))
  return(t)
}, mc.cores = 1
) # eo mclapply - f in files
# Rbind
traits_gam <- bind_rows(res)
rm(res,files) ; gc()
head(traits_gam) ; dim(traits_gam)
### -----------------------------------------------------------------------
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
### -----------------------------------------------------------------------
#restructure the trait dataset
m_traits_gam <- melt(traits_gam, id.vars = c("group","species","var"))
unique(m_traits_gam$variable)
# Create every combination of var x variable in a factor column and dcast
m_traits_gam$new.col <- factor(paste(m_traits_gam$var ,m_traits_gam$variable, sep = "_"))
unique(m_traits_gam$new.col)
# Dcast to put those new.col as columns (return to wide format)
d_traits_gam <- dcast(m_traits_gam[,c(1,2,5,6)], group + species ~ new.col, value.var = "value") 
# Select the columns corresponding to: center2, width2 and weight.pca (weights to give to the columns, not the individuals)
cols2keep <- c(colnames(d_traits_gam)[grep("center2", colnames(d_traits_gam))],
               colnames(d_traits_gam)[grep("width2", colnames(d_traits_gam))],
               colnames(d_traits_gam)[grep("weight.pca", colnames(d_traits_gam))]
) # eo cols2keep
cols2keep
# Retai columns used for PCA
data4pca <- d_traits_gam[,c("group","species",cols2keep)]
dim(data4pca) ; head(data4pca)
### -----------------------------------------------------------------------
### Check how to add weights properly in PCA
# ?PCA
vars.weights <- data.frame(traits_gam %>% group_by(var) %>% summarize(w = mean(weight.pca), probs = mean(prob.range)) )
vars.weights
vars.weights$w2 <- vars.weights$w/max(vars.weights$w)
#vars.weights[order(vars.weights$w2, decreasing = T),] # OK 
### colnames in data4pca and variables in vars.weight$w2 follow the same order, so just double the vars.weight$w2 vector
col.weights <- c(vars.weights$w2,vars.weights$w2) ; col.weights
colnames(data4pca)

data4pca[is.na(data4pca$group),]$group <- 12
data4pca <- na.omit(data4pca)

pca.phy <- PCA(na.omit(data4pca[3:18]), scale.unit = T, ncp = 7, graph = T, col.w = col.weights)
str(pca.phy)
data4pca[,paste("PC",c(1:5),sep="")] <- pca.phy$ind$coord[,c(1:5)]
cbind(data4pca[3:18], pca.phy$ind$coord[,c(1:5)])
### -----------------------------------------------------------------------
eig <- data.frame(perc = pca.phy$eig[,"percentage of variance"], nb = c(1:nrow(pca.phy$eig))) # eig
pca1 <- paste0("PC1 (",floor(eig$perc[1]*100)/100,"%)")
pca2 <- paste0("PC2 (",floor(eig$perc[2]*100)/100,"%)")
pca3 <- paste0("PC3 (",floor(eig$perc[3]*100)/100,"%)")
pca4 <- paste0("PC4 (",floor(eig$perc[4]*100)/100,"%)")
pca5 <- paste0("PC5 (",floor(eig$perc[5]*100)/100,"%)")
pca6 <- paste0("PC6 (",floor(eig$perc[6]*100)/100,"%)")
pca7 <- paste0("PC7 (",floor(eig$perc[7]*100)/100,"%)")
### -----------------------------------------------------------------------
groups.coords <- data.frame(data4pca %>% group_by(group) %>% summarize(meanPC1 = mean(PC1), meanPC2 = mean(PC2), meanPC3 = mean(PC3), meanPC4 = mean(PC4),
                                                                       se1  = 2*(sd(PC1)/length(PC1)),  se2  = 2*(sd(PC2)/length(PC2)), se3  = 2*(sd(PC3)/length(PC3)),  se4  = 2*(sd(PC4)/length(PC4))))


str(groups.coords)
groups.coords$group <- as.factor(groups.coords$group)
data4pca$group <- as.factor(groups.coords$group)
pcad <- augment.PCA(pca.phy)

#plot niche space
cent1 <- ggplot() + geom_point(aes(x = PC1, y = PC2, fill = factor(group)), data = data4pca, pch = 21, size = 2) + 
  geom_point(aes(x = meanPC1, y = meanPC2, fill = factor(group)), data = groups.coords, pch = 21, size = 6) + 
  xlab(pca1) + ylab(pca2) + 
  geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") + 
  theme_classic() +
  geom_errorbar(data=groups.coords, mapping=aes(x = meanPC1, ymin=meanPC2 + se2, ymax=meanPC2 - se2, color = factor(group)), width = 0.08, size=1) + 
  geom_errorbarh(data=groups.coords, mapping=aes(y = meanPC2, xmin=meanPC1 + se1, xmax=meanPC1 - se1, color = factor(group)), height = 0.11, size=1) +
  scale_fill_brewer(name = "", palette = "Paired") + scale_color_brewer(name = "", palette = "Paired")

setwd(wd_plots)
ggsave(plot = cent1, filename = paste("GAM_niche_space_paper_raw_1_2.svg", sep = ""), width = 12, height = 7, dpi = 300)


cent2 <- ggplot() + geom_point(aes(x = PC3, y = PC4, fill = factor(group)), data = data4pca, pch = 21, size = 2) + 
  geom_point(aes(x = meanPC3, y = meanPC4, fill = factor(group)), data = groups.coords, pch = 21, size = 6) + 
  xlab(pca3) + ylab(pca4) + 
  geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") + 
  theme_classic() +
  geom_errorbar(data=groups.coords, mapping=aes(x = meanPC3, ymin=meanPC4 + se4, ymax=meanPC4 - se4, color = factor(group)), width = 0.08, size=1) + 
  geom_errorbarh(data=groups.coords, mapping=aes(y = meanPC4, xmin=meanPC3 + se3, xmax=meanPC3 - se3, color = factor(group)), height = 0.11, size=1) +
  scale_fill_brewer(name = "", palette = "Paired") + scale_color_brewer(name = "", palette = "Paired")

setwd(wd_plots)
ggsave(plot = cent2, filename = paste("GAM_niche_space_paper_raw_3_4.svg", sep = ""), width = 12, height = 7, dpi = 300)

#plot pca variable plot
pcad <- augment.PCA(pca.phy)
pcad$var # Re-name properly
pcad$var <- str_replace_all(pcad$var, "_", " ")
pcad$var <- str_replace_all(pcad$var, "center2", "center")
pcad$var <- str_replace_all(pcad$var, "width2", "width")



#display.brewer.pal(8, "Dark2")

pcad1 <- ggplot(pcad) +
  coord_fixed() + scale_x_continuous(breaks=0) + scale_y_continuous(breaks=0) +
  annotate(geom="path", colour="black", x=cos(seq(0, 2*pi, length.out=100)), y=sin(seq(0, 2*pi, length.out=100)), colour = "grey30") +
  geom_segment(aes(x=0, xend = Dim.1, y=0, yend = Dim.2), arrow=arrow(angle=20, length=unit(0.01, "npc"))) +
  geom_text_repel(aes(x=Dim.1, y=Dim.2, label=var), 
                  #data=filter(pcad, (Dim.1^2+Dim.2^2) > 0.2^2), segment.alpha=0.5) +
  ) + xlab(pca1) + ylab(pca2) + theme_bw()

setwd(wd_plots)
ggsave(plot = pcad1, filename = paste("GAM_niche_space_vars_1_2.svg", sep = ""), width = 9, height = 9, dpi = 300)


pcad2 <- ggplot(pcad) +
  coord_fixed() + scale_x_continuous(breaks=0) + scale_y_continuous(breaks=0) +
  annotate(geom="path", colour="black", x=cos(seq(0, 2*pi, length.out=100)), y=sin(seq(0, 2*pi, length.out=100)), colour = "grey30") +
  geom_segment(aes(x=0, xend = Dim.3, y=0, yend = Dim.4), arrow=arrow(angle=20, length=unit(0.01, "npc"))) +
  geom_text_repel(aes(x=Dim.3, y=Dim.4, label=var), 
                  #data=filter(pcad, (Dim.1^2+Dim.2^2) > 0.2^2), segment.alpha=0.5) +
  ) + xlab(pca3) + ylab(pca4) + theme_bw()

ggsave(plot = pcad2, filename = paste("GAM_niche_space_vars_3_4.svg", sep = ""), width = 9, height = 9, dpi = 300)







