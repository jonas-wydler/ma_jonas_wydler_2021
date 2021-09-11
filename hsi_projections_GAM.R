#Global copepod habitat suitability projections using SDMs, Jonas Wydler 18.02.2020, updated 11.09.2021
# - Uses copepod community tables (from the SDMs) to calculate mean habitat suitability across species
# Replace ANN everywhere for GLM and ANN or use the provided scripts.
# -------------------------------------------------------------------------------------------------
Sys.setenv(LANGUAGE= 'en')

library("tidyverse")
library("reshape2")
library("biomod2")
library("viridis")
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
wd_script <- ("/data/")#where to save plots

wd_dat <- ("/data/") #wd where intermediate outputs are stored for ensemble
wd_map <- ("/data/") #wd where you want to save the map
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
dat_traits <- read.csv2("table_funct_traits_copepods_v2.csv")
#dat_traits<- dat_traits[!is.na(dat_traits$max_body_length),]
# dat_size <- dat_traits[c("Species", "max_body_length")]
# colnames(dat_size) <- c("species", "max")
colnames(dat_traits)[9] <- "size"
colnames(dat_traits)[3] <- "species"

for (species in dat_traits){
  message(paste0(species))
  speciesname <-dat_traits$species
  #speciesname <- str_replace_all(unique(species_FG$species), "_", ".")
  speciesname <- gsub("\\(|\\)", "", speciesname)
  
  if( grepl(pattern = '-', x = speciesname, fixed = T) ) 
    speciesname <- str_replace_all(speciesname,'-','')
  
  dat_traits$species <- speciesname
  
}


#f <- files[1] #to test
res <- mclapply(files, function(f) {
  
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
  #specieswithsizes <- unique(dat_size$species) # vector f species names 
  common.spp <- intersect(colnames(comm[4:(length(comm)-1)]),colnames(comm[4:(length(comm)-1)])) # find the common species between those modelled (colnames in the monthly or annual community table and vector of names in the traits table
  # Subset the community table to retain only the species xic size information (can also adapt this code )
  subset <- comm[,c("cell_id","x","y",common.spp)] # dim(subset)
  rm(comm) ; gc()

  # dplyr approach
  m.sub <- melt(subset, id.vars = c("cell_id","x","y"))
  colnames(m.sub)[c(4,5)] <- c("species","HSI")
  m.sub$HSI <- as.numeric(m.sub$HSI)
 
  # summary(m.sub)
  # Use summarize to derive estimate of HSI weighted median body length (size structure)
  #require("Hmisc")  !!! https://stackoverflow.com/questions/33807624/understanding-ddply-error-message 
  require("matrixStats") ; require("dplyr")
  
  #FROM THE HELP SECTION
  #If na.rm is FALSE an NA or NaN value in any of the arguments will cause a value of NA or NaN to be returned, otherwise NA and NaN values are ignored.
  
  df_fd <- data.frame(m.sub %>% group_by(cell_id) %>%
    summarize(x = unique(x), y = unique(y),sum_sp = sum(HSI, na.rm = FALSE)))

  return(df_fd)
  
}, mc.cores = 1
)
# eo mclapply - f in files.base
# Rbind
table.rich <- bind_rows(res) #; rm(res) ; gc()
# head(table.med.size) ; dim(table.med.size) ; summary(table.med.size)
an.rich <- data.frame(table.rich %>% group_by(cell_id) %>% summarize(x = unique(x), y = unique(y), an.rich = mean(sum_sp, na.rm = T)))
#an.rich <- data.frame(table.rich %>% group_by(cell_id) %>% summarize(x = unique(x), y = unique(y), an.rich = mean(rowsum, na.rm = T)))
# 
# setwd(wd_dat)
# save(an.rich, file = "rich_gam_dplyr_na_f.RData")

df <- an.rich
df$x.2 <- df$x
df[df$x < 0 ,"x.2"] <- (df[df$x < 0 ,"x"]) + 360
#Here I rename the month to filler so that it works with changing names
colnames(df) <- c("cell_id","x","y","an.rich","x.2")

g2 <- ggplot() + geom_raster(aes(x = x.2, y = y, fill = an.rich), data = df) + 
  # play around with palettes. You can also use scale_fill_distiller() and
  geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
  coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                           panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
  scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
  scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
  scale_fill_viridis(na.value = "white") + labs(fill = "Annual mean Richness [sum of HSI]")


setwd("C:/Users/Jonas/ma/sdm/part3/richness_maps/")
ggsave(plot = g2, filename = paste0("Copepod_annual_richness_GAM_rowsum_dplyr_rnarmF.jpg"),dpi = 300, width = 7, height = 4)

# -------------------------------------------------------------------------------------------------
#ensemble calculations
setwd(wd_dat)
load("rich_ANN_dplyr_na_f.RData")
dat_rich_ann <- an.rich;rm(an.rich);

load("rich_GAM_dplyr_na_f.RData")
dat_rich_gam <- an.rich;rm(an.rich);

load("rich_GLM_dplyr_na_f.RData")
dat_rich_glm <- an.rich;rm(an.rich);gc()

dat_rich_ensm <- rbind(dat_rich_ann, dat_rich_gam, dat_rich_glm)
dat_rich_ensm <- data.frame(dat_rich_ensm %>% group_by(cell_id) %>% summarize(x = unique(x), y = unique(y), an.rich = mean(an.rich, na.rm = T)))

df <- dat_rich_ensm
df$x.2 <- df$x
df[df$x < 0 ,"x.2"] <- (df[df$x < 0 ,"x"]) + 360
#Here I rename the month to filler so that it works with changing names

g2 <- ggplot() + geom_raster(aes(x = x.2, y = y, fill = an.rich), data = df) + 
  # play around with palettes. You can also use scale_fill_distiller() and
  geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
  coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                           panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
  scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
  scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
  scale_fill_viridis(na.value = "white") + labs(fill = "Annual mean Richness [sum of HSI]")

setwd(wd_map)
ggsave(plot = g2, filename = "Copepod_annual_richness_ensm_rowsum_dplyr_rnarmF.jpg", width = 7, height = 4, dpi = 300)

