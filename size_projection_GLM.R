# -------------------------------------------------------------------------------------------------
#Global copepod body size projections using SDMs and trait table, Jonas Wydler 14.12.2020, updated 11.09.2021
# - Uses copepod community tables (from the SDMs) and copepod body size information to project a CW median body size for the globe
# Replace GLM everywhere for GAM and ANN or use the provided scripts.
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
wd_ct_GLM <- ("/data/community_tables/ct_GLM")
wd_omit <- ("/data/") #wd with the text file on TSS scores
wd_clim <- ("/data/global_monthly_clims_1d/")
wd_size <- ("/data/") #wd with tratis table
#wd_ensemble <- ("C:/Users/Jonas/ma/sdm/part3/data_for_ensemble/size/")
wd_script <- ("/data/")#where to save plots

#load species to omit based on TSS score
setwd(wd_omit)
TSS_dat <- read.table("species_omit_TSS_v2.txt", h = T, sep = "\t")
TSS_means_omit <- subset(TSS_dat, TSS < 0.3)
TSS_means_omit_GLM <- subset(TSS_means_omit, SDM  == "GLM")


#projections
setwd(wd_ct_GLM)
files <- dir()[grep("table_mon_composition_baseline_GLM_",dir())]

##############################################
#Read size
setwd(wd_size)
dat_traits <- read.csv2("table_funct_traits_copepods_v2.csv")
dat_size <- dat_traits[c("Species", "max_body_length")]
colnames(dat_size) <- c("species", "max")

for (species in dat_size){
  message(paste0(species))
  speciesname <-dat_size$species
  #speciesname <- str_replace_all(unique(species_FG$species), "_", ".")
  speciesname <- gsub("\\(|\\)", "", speciesname)
  
  if( grepl(pattern = '-', x = speciesname, fixed = T) ) 
    speciesname <- str_replace_all(speciesname,'-','')
  
  dat_size$species <- speciesname
  
}
dat_size <- na.omit(dat_size)


files.base <- files
# files.base is a vector containing the names of the files containing the monthly/annual HSI community tables (you should have one per month and SDM, so 36?)
# each row is a grid cell and columns should be: species HSI, coordinates, cell_id etc.
res <- mclapply(files.base, function(f) {
  
  # Message
  message(paste(f, sep = ""))
  name_month <- substr(f,nchar(f)-6,nchar(f)-4)
  setwd(wd_ct_GLM)
  comm <- read.table(file =f , sep = "\t", h = T)
  #comm <- get(load(f)) # dim(comm) ; colnames(comm)
  # Change colnames (remove brackets or points in the species names if necessary)
  colnames(comm) <-  gsub("(\\.|\\.)", "", colnames(comm))
  
  comm <- comm %>% dplyr::select(-c(TSS_means_omit_GLM$species))
  
  # Colnames of copepods should match the species names in the table where you store the max body length data
  specieswithsizes <- unique(dat_size$species) # vector f species names 
  common.spp <- intersect(specieswithsizes,colnames(comm[4:(length(comm)-1)])) # find the common species between those modelled (colnames in the monthly or annual community table and vector of names in the traits table
  # Subset the community table to retain only the species xic size information (can also adapt this code )
  subset <- comm[,c("cell_id","x","y",common.spp)] # dim(subset)
  rm(comm) ; gc()
  # Melt to long format, to put species names as vector
  m.sub <- melt(subset, id.vars = c("cell_id","x","y"))
  colnames(m.sub)[c(4,5)] <- c("species","HSI")
   # Provide species body size data to m.sub table
            m.sub$size <- NA
            for(sp in unique(m.sub$species)) {
                s <- dat_size[dat_size$species == sp,"max"]
                m.sub[m.sub$species == sp,"size"] <- s
            } # eo for loop
            # Check
            # summary(m.sub)
            # Use summarize to derive estimate of HSI weighted median body length (size structure)
            #require("Hmisc")  !!! https://stackoverflow.com/questions/33807624/understanding-ddply-error-message 
            require("matrixStats") ; require("dplyr")
            # weightedMedian(x = m.sub$size, w = m.sub$HSI)
            # Use weightedMean if you prefer to an evarged instead of a median
            med.size <- data.frame(m.sub %>% group_by(cell_id) %>% summarize(x = unique(x), y = unique(y), med.size = weightedMedian(x = size, w = HSI, na.rm = F)))
            # head(med.size) ; summary(med.size)
            # Quick map to check
            
            df4 <- med.size
            df4$x.2 <- df4$x
            df4[df4$x < 0 ,"x.2"] <- (df4[df4$x < 0 ,"x"]) + 360
            #Here I rename the month to filler so that it works with changing names
            colnames(df4) <- c("cell_id","x","y","med.size","x.2")
            
            g1 <- ggplot() + geom_raster(aes(x = x.2, y = y, fill = med.size), data = df4) +
              # play around with palettes. You can also use scale_fill_distiller() and
              geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
              coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                                       panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
              scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
              scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
              scale_fill_viridis(na.value = "white") + labs(fill = "Copepod_nmedian size_n(mm)")
            setwd(wd_script)
            #ggsave(plot = g1, filename = paste0("Copepod_nmedian size_n(mm)_map_",name_month,".jpg"),dpi = 300, width = 7, height = 4)
                                                
           

            return(med.size)
    
        }, mc.cores = 1
)
# eo mclapply - f in files.base
# Rbind
table.med.size <- bind_rows(res) #; rm(res) ; gc()
# head(table.med.size) ; dim(table.med.size) ; summary(table.med.size)
# Derive ensemble estimate of median size structure
med.size2 <- data.frame(table.med.size %>% group_by(cell_id) %>% summarize(x = unique(x), y = unique(y), med.size = mean(med.size, na.rm = T)))

df <- med.size2
df$x.2 <- df$x
df[df$x < 0 ,"x.2"] <- (df[df$x < 0 ,"x"]) + 360
#Here I rename the month to filler so that it works with changing names
colnames(df) <- c("cell_id","x","y","med.size","x.2")

g2 <- ggplot() + geom_raster(aes(x = x.2, y = y, fill = med.size), data = df) + 
    # play around with palettes. You can also use scale_fill_distiller() and
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                             panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
    scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    scale_fill_viridis(na.value = "white") + labs(fill = "Copepod\nmedian size\n(mm)")

#ggsave(plot = g2, filename = paste0("Copepod_nmedian_size_n(mm)_map_annual_GLM_narmF.jpg"),dpi = 300, width = 7, height = 4)
setwd(wd_ensemble)
annual_size_GLM_naF <- med.size2
#save(annual_size_GLM_naF,file = "annual_size_GLM_naF")

