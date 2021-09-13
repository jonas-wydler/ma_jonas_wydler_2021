# -------------------------------------------------------------------------------------------------
#map functional diversity indices FEve, FDis, RaoQ, 12.01.2021, Jonas Wydler
# 
# - change GAM to ANN or GLM for all SDMs
#
# -------------------------------------------------------------------------------------------------
Sys.setenv(LANGUAGE= 'en')

library("tidyverse")
library("reshape2")
library("biomod2")
library("viridis")
#library("vegan")
library("FactoMineR")
library("dplyr")
library("parallel")
library("FD")
library("missMDA")
#global var
months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
world2 <- map_data("world2") # world coastline for maps 
world1 <- map_data("world")

wd_back <- ("/data/species_background_data/")
wd_ct_GAM <- ("/data/community_tables/ct_GAM")
wd_scores <- ("/data/Scores")
wd_omit <- ("/data/")
wd_clim <- ("/data/global_monthly_clims_1d/")
wd_script <- ("/data/")
wd_size <- ("/data/")

#load species to omit based on TSS score
setwd(wd_omit)
TSS_dat <- read.table("species_omit_TSS_v2.txt", h = T, sep = "\t")
TSS_means_omit <- subset(TSS_dat, TSS < 0.3)
TSS_means_omit_GAM <- subset(TSS_means_omit, SDM  == "GAM")


#projections
setwd(wd_ct_GAM)
files <- dir()[grep("table_mon_composition_baseline_GAM_",dir())]


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
load("df_GAM_v2.Rda")


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

# compfamd <- imputeFAMD(fct[,c(7:9,11:14,16:18)], npc = 4)
# FAMD <- FAMD(fct[,c(7:9,11:14,16:18)], tab.disj = compfamd$tab.disj)
# 
# eigFAMD<- data.frame(perc = FAMD$eig[,"percentage of variance"], nb = c(1:nrow(FAMD$eig)) )
# 


#f <- files[2]



f <- files[12]
res <- mclapply(files, function(f) {
  
  message(paste(f, sep = ""))
  name_month <- substr(f,nchar(f)-6,nchar(f)-4)
  setwd(wd_ct_GAM)
  comm <- read.table(file =f , sep = "\t", h = T)
  #comm <- get(load(f)) # dim(comm) ; colnames(comm)
  # Change colnames (remove brackets or points in the species names if necessary)
  colnames(comm) <-  gsub("(\\.|\\.)", "", colnames(comm))
  setwd(wd_script)
  
  comm <- comm %>% dplyr::select(-c(TSS_means_omit_GAM$species))
  
  # Colnames of copepods should match the species names in the table where you store the max body length data
  
  specieswithtrait<- unique(dat_traits_fd_v2$species) # vector f species names 
  common.spp <- intersect(specieswithtrait, colnames(comm[4:(length(comm)-1)])) # find the common species between those modelled (colnames in the monthly or annual community table and vector of names in the traits table
  
  #transform to 1/0
  species_tss_names <- unique(df_GAM_v2$species)
  common.spp2 <- intersect(species_tss_names, common.spp) 
  common.spp2 <- sort(common.spp2)
  # Subset the community table to retain only the species xic size information (can also adapt this code )
  subset <- comm[,c("cell_id","x","y",common.spp2)] # dim(subset)
  df_GAM_v2_subset <- subset(df_GAM_v2, is.element(df_GAM_v2$species,common.spp2))
  
  #sorting
  df_GAM_v2_subset$species <- as.character(df_GAM_v2_subset$species)
  df_GAM_v2_subset <- df_GAM_v2_subset[order(df_GAM_v2_subset$species),]
  
  df_sort <- subset[c(common.spp2)]
  df_sort <- df_sort[,order(colnames(df_sort))]
  subset[common.spp2] <- df_sort
  
  #carful, the order of the species names has to be the same, not a perfect solution
  #df_GAM_v2_subset$species == colnames(subset[4:length(subset)]) #this is to test
  
  
  df_GAM_v2_subset$tss_cutoff <- df_GAM_v2_subset$tss_cutoff/1000
  #subset2 <- subset
  subset[4:length(subset)] <- BinaryTransformation(subset[4:length(subset)],df_GAM_v2_subset$tss_cutoff )
  #abun <- subset2[4:length(subset2)]
  abun <- subset[4:length(subset)]
  # Melt to long format, to put species names as vector
  #m.sub <- melt(subset2, id.vars = c("cell_id","x","y"))
  m.sub <- melt(subset, id.vars = c("cell_id","x","y"))
  colnames(m.sub)[c(4,5)] <- c("species","presence")
  
  #df_fd <- data.frame(m.sub %>% group_by(cell_id) %>% 
  #summarize(x = unique(x), y = unique(y),sum_sp = sum(presence)))
  
  traits <- dat_traits_fd_v2[dat_traits_fd_v2$species %in% common.spp2, ]
  traits_v2 <- traits[,c("size", "Myelination","Spawning","Omnivore","Carnivore","Herbivore","Detritivore","Current","Cruise","Ambush")]
  rownames(traits_v2) <- traits$species
  
  #subset out cell ids that are completly 0 after NAs are beeing set to 0
  abun_na <- abun
  abun_na$na_count <- apply(abun, 1, function(x) sum(is.na(x)))
  
  #this part should later be done on the subset with cell ids so it can later be reconstructed
  abun_na <- abun_na[abun_na$na_count < 1 ,]
  abun_na <- abun_na[1:length(abun_na)-1]

  # -------------------------------------------------------------------------------------------------
  #FAMD
  #npfamd <- estim_ncpFAMD(traits_v2)
  pdf(file = NULL)
  compfamd <- imputeFAMD(traits_v2, npc = 4)
  FAMD <- FAMD(traits_v2, tab.disj = compfamd$tab.disj)
  dev.off()
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
  # -------------------------------------------------------------------------------------------------

 
  ex1 <- dbFD(x = famd_dist, a = abun_na, calc.FRic = T)

  #add xy coords to ex1
  subset_na <- subset
  subset_na$na_count <- apply(abun, 1, function(x) sum(is.na(x)))

  subset_na <- subset_na[subset_na$na_count < 1 ,]
  subset_na <- subset_na[1:length(subset_na)-1]

  df_ann<- as.data.frame(ex1)
  df_ann_v2 <- cbind(subset_na[1:3],df_ann)
  setwd(wd_script)
  save(df_ann_v2, file = paste0("dat_convex_binary_",name_month,"_GAM.Rda"))
  
  return(subset)
  
}, mc.cores = 12
)
# -------------------------------------------------------------------------------------------------
#test FD results
setwd(wd_script)

table_df <- bind_rows(res)
#test <- res[[1]]
save(table_df, file = paste0("dat_convex_annual_binary_GAM.Rda"))




