#Rscript predictor rankings, Master thesis Jonas Wydler, 16.11.2020
#updated 03.09.2021
# - Plot predictor rankings/scores
#-----------------------------------
library("tidyverse")
library("parallel") 
library("ggpubr")
library("reshape2") 
library("svglite")
#-----------------------------------
setwd(paste0("wd_data/Ranks/"))

ranks <- dir()[grep("table_ranks_", dir())] ; ranks

setwd(paste0("wd_data/Scores/"))
scores <- dir()[grep("table_scores_", dir())] ; scores

setwd(paste0("wd_data/"))
load("species_FG")
colnames(species_FG) <- c("species","FG_f11")


for (species in species_FG){

  speciesname <-species_FG$species
  #speciesname <- str_replace_all(unique(species_FG$species), "_", ".")
  speciesname <- gsub("\\(|\\)", "", speciesname)
  
  if( grepl(pattern = '-', x = speciesname, fixed = T) ) 
    speciesname <- str_replace_all(speciesname,'-','')
  
  species_FG$species <- speciesname

}
setwd(paste0("wd_data/"))

names <- read.table(file = "FgNames.txt", sep = ";", h = TRUE)
species_FG <- base::merge(species_FG,names, by  = "FG_f11", all.x = TRUE)

#===================================================================
res <- lapply(ranks,function(c){
  setwd(paste0("wd_data/Ranks/"))
  message(paste("Reading Ranks ",c, sep = ""))
  data <- load(c)#carful with name, called m.ranks for me
  
  return(m.ranks)
})
table_ranks <- bind_rows(res); rm(res); rm(m.ranks); gc()

#restlul <- subset(table_ranks, species == "Scaphocalanus_echinatus")

#same for scores
res <- lapply(scores,function(c){
  setwd(paste0("wd_data/Scores/"))
  message(paste("Reading Scores ",c, sep = ""))
  data <- load(c)#carful with name, called scores.tb1 for me
  
  return(scores.tbl)
})
table_scores <- bind_rows(res); rm(res); rm(scores.tbl); gc()
table_scores <- tibble::rownames_to_column(table_scores, "model")


#this part is to get SDM and RUN out of rowname, it's not ideal and these variables should be defined in the model directly
table_scores$temp <- substr(x = as.character(table_scores$model), start = 1, stop =8)
table_scores$SDM <- data.frame(do.call(rbind,strsplit(as.character(table_scores$temp), split = "_")))[,1]
table_scores$RUN <- data.frame(do.call(rbind,strsplit(as.character(table_scores$temp), split = "_")))[,2]
table_scores$temp <-NULL
#===================================================================
#assign FG
ddf_fg <- base::merge(table_ranks,species_FG, by  = "species", all.x = TRUE)
ddf2_fg <- base::merge(table_scores,species_FG, by  = "species", all.x = TRUE)

summary(ddf_fg)

#===================================================================
#Plots
result11 <- as.data.frame(ddf_fg %>% group_by(FG_f11) %>% summarize(fg = unique(FG_f11), count = length(unique(species))))



# a single plot showing the distribution of ranks across all species and all SDMs 
p1 <- ggplot(ddf_fg, aes(x=rank, y=fct_reorder(predictor,rank, .fun = mean, .desc = FALSE, na.rm = TRUE), fill = predictor)) + 
  geom_boxplot() + xlab("") +#  theme(axis.text=element_text(size=12),axis.title=element_text(size=16)) + 
  ylab("") + theme_minimal() + theme(legend.position="none") 

ggsave(plot = p1, filename = paste("plot1_ranks_predictorsv2.svg"),dpi = 300, width = 6, height = 4)

#a panel of plots showing the distribution of ranks across sdms (facet per sdm)
p2 <- ggplot(ddf_fg, aes(x=rank, y=fct_reorder(predictor,rank, .fun = mean, .desc = FALSE, na.rm = TRUE), fill=SDM)) + 
  geom_boxplot() +
  facet_wrap(~SDM, scale="free_y") + xlab("") + 
  ylab("")+  theme_minimal() + theme(legend.position="none")
ggsave(plot = p2, filename = paste("plot2_ranks_predictors_forSDM.svg"),dpi = 600, width = 16, height = 9)

#a larger panel of plots showing the distribution of ranks across FGs 
p3 <- ggplot(ddf_fg, aes(x=rank, y=fct_reorder(predictor,rank, .fun = median, .desc = FALSE, na.rm = TRUE), fill=FG_desc)) + 
  geom_boxplot() +
  facet_wrap(~FG_desc, scale="free_y") + xlab("") + 
  ylab("") + theme_minimal() + theme(legend.position="none")

ggsave(plot = p3, filename = paste("plot3_ranks_predictors_forFG.svg"),dpi = 600, width = 16, height = 9)

#a panel plot showing the distribution of TSS across all species, facet per SDM 
p4 <- ggplot(ddf2_fg, aes(x=TSS, y=fct_reorder(SDM,TSS, .fun = mean, .desc = FALSE, na.rm = TRUE), fill = SDM)) + 
  geom_boxplot()  + xlab("") + 
  ylab("") + theme_minimal() + theme(legend.position="none")

ggsave(plot = p4, filename = paste("plot4_TSS_forSDM.png"),dpi = 600, width = 16, height = 9)

#a larger panel plot showing the distribution of TSS across FG and SDM, facet per SDM and FG 
p5 <- ggplot(ddf2_fg, aes(x=TSS, y=fct_reorder(SDM,TSS, .fun = mean, .desc = FALSE, na.rm = TRUE), fill=FG_desc)) + 
  geom_boxplot() +
  facet_wrap(~FG_desc, scale="free_y") + xlab("") + 
  ylab("") + theme_minimal() + theme(legend.position="none")

ggsave(plot = p5, filename = paste("plot5_TSS_forFG&SDM.svg"),dpi = 600, width = 16, height = 9)
