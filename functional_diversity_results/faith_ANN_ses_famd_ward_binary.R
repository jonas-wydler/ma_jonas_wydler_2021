# -------------------------------------------------------------------------------------------------
#map standardized effect size of functional richness via faith index, 19.01.2021
# Jonas Wydler
#updated 13.09.2021
# -------------------------------------------------------------------------------------------------
Sys.setenv(LANGUAGE= 'en')


library(stats)
library(ape)
library("tidyverse")
library("reshape2")
library("viridis")
library("FactoMineR")
library("dplyr")
library("parallel")
library("FD")
library("missMDA")
library("picante")#this is to calc faith index
#global var
months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
world2 <- map_data("world2") # world coastline for maps 
world1 <- map_data("world")

#kryo
wd_back <- ("/data/species_background_data/")
wd_ct_ANN <- ("/data/community_tables/ct_ANN")
wd_scores <- ("/data/Scores")
wd_omit <- ("/data/")
wd_clim <- ("/data/global_monthly_clims_1d/")
wd_script <- ("/data/")
wd_size <- ("/data/")


#load species to omit based on TSS score
setwd(wd_omit)
TSS_dat <- read.table("species_omit_TSS_v2.txt", h = T, sep = "\t")
TSS_means_omit <- subset(TSS_dat, TSS < 0.3)
TSS_means_omit_ANN <- subset(TSS_means_omit, SDM  == "ANN")

#projections
setwd(wd_ct_ANN)
files <- dir()[grep("table_mon_composition_baseline_ANN_",dir())]


##############################################
#Read traits
setwd(wd_size)
dat_traits <- read.csv2("table_funct_traits_copepods_v2.csv")
dat_traits_fd <- dat_traits[,c(3:20)]
names <- colnames(dat_traits_fd)[c(7,8,9,10,15)] ; names
dat_traits_fd$na_count <- apply(dat_traits_fd[,names], 1, function(x) sum(is.na(x)))


dat_traits_fd <- dat_traits_fd[!is.na(dat_traits_fd$max_body_length),]
dat_traits_fd <- dat_traits_fd[dat_traits_fd$na_count < 2,]

colnames(dat_traits_fd)[7] <- "size"
colnames(dat_traits_fd)[1] <- "species"

#fix names
for (species in dat_traits_fd){
  #message(paste0(species))
  speciesname <-dat_traits_fd$species
  #speciesname <- str_replace_all(unique(species_FG$species), "_", ".")
  speciesname <- gsub("\\(|\\)", "", speciesname)
  
  if( grepl(pattern = '-', x = speciesname, fixed = T) ) 
    speciesname <- str_replace_all(speciesname,'-','')
  
  dat_traits_fd$species <- speciesname
  
}


#load TSS data
setwd(wd_script)
load("df_ANN_v2.Rda")


dat_traits_fd_v2 <- dat_traits_fd[c("species", "n","size", "Myelination", "Spawning", "Omnivore","Carnivore","Herbivore", "Detritivore", "Current", "Cruise","Ambush")]
#str(dat_traits_fd_v2)
dat_traits_fd_v2$Spawning <- as.factor(dat_traits_fd_v2$Spawning)
dat_traits_fd_v2$Myelination <- as.factor(dat_traits_fd_v2$Myelination)
dat_traits_fd_v2$Omnivore <- as.factor(dat_traits_fd_v2$Omnivore)
dat_traits_fd_v2$Carnivore <- as.factor(dat_traits_fd_v2$Carnivore)
dat_traits_fd_v2$Herbivore <- as.factor(dat_traits_fd_v2$Herbivore)
dat_traits_fd_v2$Detritivore <- as.factor(dat_traits_fd_v2$Detritivore)
dat_traits_fd_v2$Current <- as.factor(dat_traits_fd_v2$Current)
dat_traits_fd_v2$Cruise <- as.factor(dat_traits_fd_v2$Cruise)
dat_traits_fd_v2$Ambush <- as.factor(dat_traits_fd_v2$Ambush)

#sorting
dat_traits_fd_v2 <- dat_traits_fd_v2[order(dat_traits_fd_v2$species),]


# -------------------------------------------------------------------------------------------------
#f <- files[1]
res <- mclapply(files, function(f) {
  
  #ptm <- proc.time()
  
  message(paste(f, sep = ""))
  name_month <- substr(f,nchar(f)-6,nchar(f)-4)
  setwd(wd_ct_ANN)
  comm <- read.table(file =f , sep = "\t", h = T)
  #comm <- get(load(f)) # dim(comm) ; colnames(comm)
  # Change colnames (remove brackets or points in the species names if necessary)
  setwd(wd_script)
  colnames(comm) <-  gsub("(\\.|\\.)", "", colnames(comm))
  
  comm <- comm %>% dplyr::select(-c(TSS_means_omit_ANN$species))
  
  # Colnames of copepods should match the species names in the table where you store the max body length data
  
  specieswithtrait<- unique(dat_traits_fd_v2$species) # vector f species names 
  common.spp <- intersect(specieswithtrait, colnames(comm[4:(length(comm)-1)])) # find the common species between those modelled (colnames in the monthly or annual community table and vector of names in the traits table
  
  #transform to 1/0
  species_tss_names <- unique(df_ANN_v2$species)
  common.spp2 <- intersect(species_tss_names, common.spp) 
  common.spp2 <- sort(common.spp2)
  # Subset the community table to retain only the species with trait information 
  subset <- comm[,c("cell_id","x","y",common.spp2)] # dim(subset)
  df_ANN_v2_subset <- subset(df_ANN_v2, is.element(df_ANN_v2$species,common.spp2))
  
  #sorting
  df_ANN_v2_subset$species <- as.character(df_ANN_v2_subset$species)
  df_ANN_v2_subset <- df_ANN_v2_subset[order(df_ANN_v2_subset$species),]
  
  df_sort <- subset[c(common.spp2)]
  df_sort <- df_sort[,order(colnames(df_sort))]
  subset[common.spp2] <- df_sort
  
  #carful, the order of the species names has to be the same, not a perfect solution
  #df_ANN_v2_subset$species == colnames(subset[4:length(subset)]) #this is to test
  
  
  df_ANN_v2_subset$tss_cutoff <- df_ANN_v2_subset$tss_cutoff/1000
  #subset2 <- subset #subset2 is used for binary ct
  subset[4:length(subset)] <- BinaryTransformation(subset[4:length(subset)],df_ANN_v2_subset$tss_cutoff )

  # Melt to long format, to put species names as vector
  
  #m.sub <- melt(subset, id.vars = c("cell_id","x","y"))
  m.sub <- melt(subset, id.vars = c("cell_id","x","y"))
  colnames(m.sub)[c(4,5)] <- c("species","presence")
  
  #df_fd <- data.frame(m.sub %>% group_by(cell_id) %>% 
  #summarize(x = unique(x), y = unique(y),sum_sp = sum(presence)))
  
  traits <- dat_traits_fd_v2[dat_traits_fd_v2$species %in% common.spp2, ]
  traits_v2 <- traits[,c("size", "Myelination","Spawning","Omnivore","Carnivore","Herbivore","Detritivore","Current","Cruise","Ambush")]
  rownames(traits_v2) <- traits$species
  
  #
  subset$na_count <- apply(subset[4:length(subset)], 1, function(x) sum(is.na(x)))
  
  #this part should later be done on the subset with cell ids so it can later be reconstructed
  subset <- subset[subset$na_count < 1 ,]
  subset <- subset[1:length(subset)-1]
  #abun <- subset2[4:length(subset2)]
  abun <- subset[4:length(subset)]
  

  # -------------------------------------------------------------------------------------------------
  #FAMD
  #npfamd <- estim_ncpFAMD(traits_v2)
  pdf(file = NULL)
  compfamd <- imputeFAMD(traits_v2, npc = 4)
  FAMD <- FAMD(traits_v2, tab.disj = compfamd$tab.disj, graph = F)

  eigFAMD<- data.frame(perc = FAMD$eig[,"percentage of variance"], nb = c(1:nrow(FAMD$eig)) )

  famd1 <- paste0("FAMD 1 (",floor(eigFAMD$perc[1]*100)/100,"%)")
  famd2 <- paste0("FAMD 2 (",floor(eigFAMD$perc[2]*100)/100,"%)")
  famd3 <- paste0("FAMD 3 (",floor(eigFAMD$perc[3]*100)/100,"%)")
  famd4 <- paste0("FAMD 4 (",floor(eigFAMD$perc[4]*100)/100,"%)")


  famd <- data.frame(FAMD$ind$coord[,1:4])


  famd_sp <- data.frame(FAMD$ind$coord[,1:4])
  colnames(famd_sp) <- c("FAMD1","FAMD2","FAMD3","FAMD4")
  famd_all_temp <- rbind(famd_sp)

  famd_all_sp <- famd_all_temp[ order(row.names(famd_all_temp)), ]

  famd_dist <- dist(famd_all_sp, method = "euclidean")

  fit_famd_ward <- hclust(famd_dist, method = "ward.D2")
  dev.off()
  # -------------------------------------------------------------------------------------------------
  
  
  #ex2 <- pd(abun, as.phylo(fit_famd_avg), include.root=TRUE)
  ex2 <- ses.pd(samp = abun, tree = as.phylo(fit_famd_ward), null.model = "taxa.labels", runs = 100, include.root = T)
  
    #add geo information
  ex3 <- cbind(subset[1:3], ex2)
  #proc.time() - ptm
  #test <- data.frame(m.sub %>% group_by(cell_id) %>% summarize(x = unique(x), y = unique(y), val = sum(presence, na.rm = F)))
  # 
  df <- ex3
  df$x.2 <- df$x
  df[df$x < 0 ,"x.2"] <- (df[df$x < 0 ,"x"]) + 360
  #Here I rename the month to filler so that it works with changing names

  g1 <- ggplot() + geom_raster(aes(x = x.2, y = y, fill = pd.obs.z), data = df) +
    # play around with palettes. You can also use scale_fill_distiller() and
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                             panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
    theme(legend.text=element_text(size=8)) +
    scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    #scale_fill_viridis(na.value = "white") + labs(fill = "S.e.s PD (Faith index)") 
    scale_fill_gradient2(name = "sesFR",space = "Lab", mid = "white", low = "#3288BD", high = "#D53E4F", na.value = "white", aesthetics = "fill",guide = "colourbar")
    #scale_fill_brewer(palette = "PuOr")
  setwd(wd_script)
  ggsave(plot = g1, filename = paste0("FR_ses_binary_ANN_",name_month,"_famd_ward_v2.jpg"),dpi = 300, width = 7, height = 4)


  #Here I rename the month to filler so that it works with changing names

  g2 <- ggplot() + geom_raster(aes(x = x.2, y = y, fill = pd.obs.p), data = df) +
    # play around with palettes. You can also use scale_fill_distiller() and
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                             panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
    theme(legend.text=element_text(size=8)) +
    scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    #scale_fill_viridis(na.value = "white") + labs(fill = "Quantiles")
    scale_fill_gradient2(name = "quantiles",space = "Lab", mid = "white", low = "#3288BD", high = "#D53E4F", na.value = "white", aesthetics = "fill",guide = "colourbar")
  
  ggsave(plot = g2, filename = paste0("Quantiles_binary_famd_avg_",name_month,"_ANN_v2.jpg"),dpi = 300, width = 7, height = 4)
  save(df, file = paste0("ses_dat_binary_ANN_", name_month,"_famd_ward_v2.Rda"))
  
  return(ex3)
  
}, mc.cores = 12
)
setwd(wd_script)

table_df <- bind_rows(res)
#test <- res[[1]]
save(table_df, file = paste0("ses_dat_faith_annual_ANN_binary_famd_ward_v2.Rda"))

table_df_an <- data.frame(table_df %>% group_by(cell_id) %>% 
                               summarize(x = unique(x), y = unique(y), val = mean(pd.obs.z, na.rm = T)))
df2 <- table_df_an
df2$x.2 <- df2$x
df2[df2$x < 0 ,"x.2"] <- (df2[df2$x < 0 ,"x"]) + 360
#Here I rename the month to filler so that it works with changing names

g3 <- ggplot() + geom_raster(aes(x = x.2, y = y, fill = val), data = df2) +
  # play around with palettes. You can also use scale_fill_distiller() and
  geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
  coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                           panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
  theme(legend.text=element_text(size=6), legend.title = element_text(size=8)) +
  scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
  scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
  scale_fill_gradient2(name = "sesFR",space = "Lab", mid = "white", low = "#3288BD", high = "#D53E4F", na.value = "white", aesthetics = "fill",guide = "colourbar")+
  #labs(title = paste0("Copepod functional group ", fg, ", (", fg_text,"), #spp: ",number_spp,", annual mean")) +
  theme(plot.title = element_text(size=8)) + theme(legend.text=element_text(size=5)) +
  theme(legend.title = element_text(size=8))

ggsave(plot = g3, filename = paste0("FR_ses_binary_annual_ANN_famd_ward_v2.jpg"),dpi = 300, width = 7, height = 4)

# table_df_an2 <- data.frame(table_df %>% group_by(cell_id) %>%
#                             summarize(x = unique(x), y = unique(y), val = mean(pd.obs.p, na.rm = T)))
# df3 <- table_df_an2
# df3$x.3 <- df3$x
# df3[df3$x < 0 ,"x.2"] <- (df3[df3$x < 0 ,"x"]) + 360
# 
# g4 <- ggplot() + geom_raster(aes(x = x.2, y = y, fill = val), data = df3) +
#   # play around with palettes. You can also use scale_fill_distiller() and
#   geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
#   coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#                            panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
#   theme(legend.text=element_text(size=6), legend.title = element_text(size=8)) +
#   scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
#   scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
#   scale_fill_gradient2(name = "quantiles",space = "Lab", mid = "white", low = "#3288BD", high = "#D53E4F", na.value = "white", aesthetics = "fill",guide = "colourbar") + 
#   #labs(title = paste0("Copepod functional group ", fg, ", (", fg_text,"), #spp: ",number_spp,", annual mean")) +
#   theme(plot.title = element_text(size=8)) + theme(legend.text=element_text(size=5)) +
#   theme(legend.title = element_text(size=8))
# 
# ggsave(plot = g4, filename = paste0("Quantiles_binary_annual_ANN_avg.jpg"),dpi = 300, width = 7, height = 4)

