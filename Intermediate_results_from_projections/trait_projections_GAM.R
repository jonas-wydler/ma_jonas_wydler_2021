# -------------------------------------------------------------------------------------------------
#Global copepod Community weighted mean proportion of functional traits projections using SDMs and trait table, 
#Jonas Wydler 6.01.2021, updated 11.09.2021
# Replace GAM everywhere for GLM and ANN or use the provided scripts.
# -------------------------------------------------------------------------------------------------
Sys.setenv(LANGUAGE= 'en')

library("tidyverse")
library("reshape2")
library("biomod2")
library("viridis")
#library("vegan")
#library("FactoMineR")
library("dplyr")
library("parallel")
#global var
months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
world2 <- map_data("world2") # world coastline for maps 
world1 <- map_data("world")

#local
wd_back <- ("/data/species_background_data/")
wd_ct_GAM <- ("/data/community_tables/ct_GAM")
wd_omit <- ("/data/") #wd with the text file on TSS scores
wd_clim <- ("/data/global_monthly_clims_1d/")
wd_size <- ("/data/") #wd with tratis table
#wd_ensemble <- ("C:/Users/Jonas/ma/sdm/part3/data_for_ensemble/size/")
wd_script <- ("/data/")#where to save plots


#load species to omit based on TSS score
setwd(wd_omit)
TSS_dat <- read.table("species_omit_TSS_v2.txt", h = T, sep = "\t")
TSS_means_omit <- subset(TSS_dat, TSS < 0.3)
TSS_means_omit_GAM <- subset(TSS_means_omit, SDM  == "GAM")


#projections
setwd(wd_ct_GAM)
files <- dir()[grep("table_mon_composition_baseline_GAM_",dir())]


##############################################
#Read size
setwd(wd_size)
dat_traits <- read.csv("table_funct_traits_copepods_v2.csv", h = T, sep = ";", dec = ",")
#dat_size <- dat_traits[c("Species", "max_body_length")]
dat_traits_v2 <- dat_traits[c("Species", "n", "Myelination", "Spawning", "Trophism", "Feeding_mode", "max_body_length")]
colnames(dat_traits_v2) <- c("species", "n", "Myelination", "Spawning", "Trophism", "Feeding_mode", "size")
#colnames(dat_size) <- c("species", "max")

#fix names
for (species in dat_traits_v2){
  #message(paste0(species))
  speciesname <-dat_traits_v2$species
  #speciesname <- str_replace_all(unique(species_FG$species), "_", ".")
  speciesname <- gsub("\\(|\\)", "", speciesname)
  
  if( grepl(pattern = '-', x = speciesname, fixed = T) ) 
    speciesname <- str_replace_all(speciesname,'-','')
  
  dat_traits_v2$species <- speciesname
  
}
#dat_size <- na.omit(dat_size)
#dat_traits_v2 <- subset(dat_traits_v2, !is.na(dat_traits_v2$size))

weighted.proportion <- function (x, w){
 
  wones <- (x * w)/(sum(w))
  return(sum(wones))
}
files.base <- files

#Change the trait here and below for the other functional traits analysis.
trait <- "Trophism"

dat_traits_v3 <- subset(dat_traits_v2, !is.na(dat_traits_v2$Trophism))

#We decided to set Omnivore-Carnivore to Carnivore for this part of the analysis.
#Change the trait here and below for the other functional traits analysis
dat_traits_v3$Trophism[dat_traits_v3$Trophism == "Omnivore-Carnivore"]<-"Carnivore"

trophic_guilds <- unique(dat_traits_v3$Trophism)
#f <- files.base[1]

res <- mclapply(files.base, function(f) {
  
  # Message
  message(paste(f, sep = ""))
  name_month <- substr(f,nchar(f)-6,nchar(f)-4)
  setwd(wd_ct_GAM)
  comm <- read.table(file =f , sep = "\t", h = T)
  #comm <- get(load(f)) # dim(comm) ; colnames(comm)
  # Change colnames (remove brackets or points in the species names if necessary)
  colnames(comm) <-  gsub("(\\.|\\.)", "", colnames(comm))
  
  comm <- comm %>% dplyr::select(-c(TSS_means_omit_GAM$species))
  
  # Colnames of copepods should match the species names in the table where you store the max body length data
  
  specieswithtrait<- unique(dat_traits_v3$species) # vector f species names 
  common.spp <- intersect(specieswithtrait, colnames(comm[4:(length(comm)-1)])) # find the common species between those modelled (colnames in the monthly or annual community table and vector of names in the traits table
  # Subset the community table to retain only the species xic size information (can also adapt this code )
  subset <- comm[,c("cell_id","x","y",common.spp)] # dim(subset)
  
  
  rm(comm) ; gc()
  # Melt to long format, to put species names as vector
  m.sub <- melt(subset, id.vars = c("cell_id","x","y"))
  colnames(m.sub)[c(4,5)] <- c("species","HSI")
  # Provide species body size data to m.sub table
  m.sub$Trophism <- NA
  for(sp in unique(m.sub$species)) {
    trait <- dat_traits_v3[dat_traits_v3$species == sp,"Trophism"] #Change the trait here and below for the other functional traits analysis
    m.sub[m.sub$species == sp,"Trophism"] <- trait
  } # eo for loop
  # Check
  # summary(m.sub)
  # Use summarize to derive estimate of HSI weighted median body length (size structure)
  #require("Hmisc")  !!! https://stackoverflow.com/questions/33807624/understanding-ddply-error-message 
  require("matrixStats") ; require("dplyr")
  #Change the trait here for the other functional traits analysis.
  df_trait <- data.frame(m.sub %>% group_by(cell_id) %>% 
                       summarize(x = unique(x), y = unique(y),relfrq = weighted.proportion(x = ifelse(Trophism == "Omnivore-Herbivore", 1, 0),w = HSI))) #Change the trait here for the other functional traits analysis.
  # Quick map to check
  
  df4 <- df_trait
  df4$x.2 <- df4$x
  df4[df4$x < 0 ,"x.2"] <- (df4[df4$x < 0 ,"x"]) + 360
  #Here I rename the month to filler so that it works with changing names
  colnames(df4) <- c("cell_id","x","y","relfrq","x.2")

  g1 <- ggplot() + geom_raster(aes(x = x.2, y = y, fill = relfrq), data = df4) +
    # play around with palettes. You can also use scale_fill_distiller() and
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                             panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
    scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    scale_fill_viridis(na.value = "white") + labs(fill = "TODO")
  setwd(wd_script)
 
  
  
  return(df_trait)
  
}, mc.cores = 1
)
# eo mclapply - f in files.base
# Rbind
table_df_trait <- bind_rows(res) #; rm(res) ; gc()
# head(table.med.size) ; dim(table.med.size) ; summary(table.med.size)
# Derive ensemble estimate of median size structure
df_trait_2 <- data.frame(table_df_trait %>% group_by(cell_id) %>% 
                    summarize(x = unique(x), y = unique(y), relfrq = mean(relfrq, na.rm = T)))

df <- df_trait_2
df$x.2 <- df$x
df[df$x < 0 ,"x.2"] <- (df[df$x < 0 ,"x"]) + 360
#Here I rename the month to filler so that it works with changing names
colnames(df) <- c("cell_id","x","y","relfrq","x.2")

g2 <- ggplot() + geom_raster(aes(x = x.2, y = y, fill = relfrq), data = df) + 
  # play around with palettes. You can also use scale_fill_distiller() and
  geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
  coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                           panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") + 
  theme(legend.text=element_text(size=8)) +
  scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
  scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
  scale_fill_viridis(na.value = "white") + labs(fill = "Proportion of herbivorous species")


setwd(wd_script)
ggsave(plot = g2, filename = paste0("Copepod_proportion_herbivorous_map_.jpg"),dpi = 300, width = 7, height = 4)


# -------------------------------------------------------------------------------------------------
#save dfs for ensemble below
#
# df_carni <- df
# df_herbi <- df
# df_detri <- df
# df_omni <- df
# 
# test <- df_carni$relfrq + df_herbi$relfrq + df_detri$relfrq + df_omni$relfrq
# table(test)
# save(df_herbi,file="df_herbi_GAM.Rda")


#Calculate ensemble
setwd("C:/Users/Jonas/ma/sdm/part3/data_for_ensemble/")
load("df_omni_ANN.Rda")
omni_ANN <- df_omni
rm(df_omni)

load("df_omni_GLM.Rda")
omni_GLM <- df_omni
rm(df_omni)

load("df_omni_GAM.Rda")
omni_GAM <- df_omni
rm(df_omni)

omni_ensm <- bind_rows(omni_ANN,omni_GLM,omni_GAM)

df_trait_3 <- data.frame(omni_ensm %>% group_by(cell_id) %>% 
                           summarize(x = unique(x), y = unique(y), relfrq = mean(relfrq, na.rm = T)))

df <- df_trait_3
df$x.2 <- df$x
df[df$x < 0 ,"x.2"] <- (df[df$x < 0 ,"x"]) + 360
#Here I rename the month to filler so that it works with changing names
colnames(df) <- c("cell_id","x","y","relfrq","x.2")

g2 <- ggplot() + geom_raster(aes(x = x.2, y = y, fill = relfrq), data = df) + 
  # play around with palettes. You can also use scale_fill_distiller() and
  geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
  coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                           panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") + 
  theme(legend.text=element_text(size=8)) +
  scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
  scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
  scale_fill_viridis(na.value = "white") + labs(fill = "Proportion of omnivorous species")


setwd(wd_script)
ggsave(plot = g2, filename = paste0("Copepod_proportion_omnivorous_map_annual_ensemble.jpg"),dpi = 300, width = 7, height = 4)
