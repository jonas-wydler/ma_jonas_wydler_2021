# -------------------------------------------------------------------------------------------------
#15.01.2021, Jonas Wydler, 
#New implementation for copepod function group mapping
#UPDATE:27.01.2021, non-prop implementation again -> so sum of hsi per group divided by #group members
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

wd_omit <- ("C:/Users/Jonas/ma/sdm/background/low_TSS_species_projections/")
wd_clim <- ("C:/Users/Jonas/ma/sdm/global_monthly_clims_1d/")
wd_script <- ("C:/Users/Jonas/ma/sdm/part3/")
wd_fginfo <- ("C:/Users/Jonas/ma/sdm/")
wd_fgv2 <- ("C:/Users/Jonas/ma/sdm/part3/fg_maps_v2/")
wd_fd_data <- ("C:/Users/Jonas/ma/sdm/part3/data_for_ensemble/fg/GAM/")
sdm <- "GAM"

wd_back <- ("/data/species_background_data/")
wd_ct_GAM <- ("/data/community_tables/ct_GAM")
wd_omit <- ("/data/") #wd with the text file on TSS scores
wd_clim <- ("/data/global_monthly_clims_1d/")
wd_size <- ("/data/") #wd with tratis table
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
#Read FG
setwd(wd_fginfo)
load("species_FG_v2") #FG description with species that have no FG assigned being in their own "group".
colnames(species_FG_v2) <- c("species","FG_f11")

fg_desc <- read.csv2("FgNames_v2.txt")
species_FG_v_na <- subset(species_FG_v2, is.na(species_FG_v2$FG_f11))
species_FG_v_na$FG_f11 <- "NA"
levels(species_FG_v_na$FG_f11) <- "NA"

species_FG_n_na <- subset(species_FG_v2, !is.na(species_FG_v2$FG_f11))
levels(species_FG_n_na$FG_f11)

species_fg_df <- rbind(species_FG_n_na,species_FG_v_na)

#fix names


for (species in species_fg_df){
  
  speciesname <-species_fg_df$species
  #speciesname <- str_replace_all(unique(species_fg_df$species), "_", ".")
  speciesname <- gsub("\\(|\\)", "", speciesname)
  
  if( grepl(pattern = '-', x = speciesname, fixed = T) ) 
    speciesname <- str_replace_all(speciesname,'-','')
  
  species_fg_df$species <- speciesname
  
}

#fg <- "7" #to test

fgs <-  unique(species_fg_df$FG_f11)
#for loop across functional groups
for (fg in fgs){
  #f = files[1] #to test
  res <- mclapply(files, function(f) {
    message(paste(f, sep = ""))
    name_month <- substr(f,nchar(f)-6,nchar(f)-4)
    setwd(wd_ct_GAM) 
    
    
    comm <- read.table(file =f , sep = "\t", h = T)
    #comm <- get(load(f)) # dim(comm) ; colnames(comm)
    
    # Change colnames (remove brackets or points in the species names if necessary)
    colnames(comm) <-  gsub("(\\.|\\.)", "", colnames(comm))
    
    #drop species 
    comm <- comm %>% dplyr::select(-c(TSS_means_omit_GAM$species))
    specieswithfg <- unique(species_fg_df$species) 
    # find the common species between those modelled (colnames in the monthly or annual community table and vector of names in the traits table
    common.spp <- intersect(specieswithfg,colnames(comm[4:(length(comm)-1)])) 
    
    #common.spp <- setdiff(colnames(comm[4:(length(comm)-1)]),specieswithfg) 
    subset <- comm[,c("cell_id","x","y",common.spp)] # dim(subset)
    rm(comm) ; gc()
    
    
    #melt 
    m.sub <- melt(subset, id.vars = c("cell_id","x","y"))
    colnames(m.sub)[c(4,5)] <- c("species","HSI")
    
    #sp <- as.character(unique(m.sub$species)[1])
    m.sub$FG_f11 <- NA
    #m.sub$FG_f11 <- as.factor(m.sub$FG_f11)
    for(sp in unique(m.sub$species)) {
      grp <- as.character(species_fg_df[species_fg_df$species == sp,"FG_f11"])
      m.sub[m.sub$species == sp,"FG_f11"] <- as.character(grp)
    } # eo for loop
    m.sub$FG_f11 <- as.factor(m.sub$FG_f11)
    
    
   
    require("matrixStats") ; require("dplyr")
    number_spp <- length(intersect(species_fg_df[species_fg_df$FG_f11 == fg,]$species,common.spp))
    fg_text <- fg_desc[fg_desc$FG_f11 == fg,"FG_desc"][1]
    m.sub_fg <- subset(m.sub, m.sub$FG_f11 == fg)
    df_fg <- data.frame(m.sub_fg %>% group_by(cell_id) %>% 
                          summarize(x = unique(x), y = unique(y),relfrq = (sum(HSI, na.rm = F)/number_spp)))#mean (i.e. devided by #species in the group) hsi per fg 
    
    #this part is to make the map atlantic centred
    df <- df_fg
    df$x.2 <- df$x
    df[df$x < 0 ,"x.2"] <- (df[df$x < 0 ,"x"]) + 360
    #Here I rename the month to filler so that it works with changing names
    colnames(df) <- c("cell_id","x","y","relfrq","x.2")

    g1 <- ggplot() + geom_raster(aes(x = x.2, y = y, fill = relfrq), data = df) +
      # play around with palettes. You can also use scale_fill_distiller() and
      geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
      coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                               panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
      scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
      scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
      scale_fill_viridis(na.value = "white") + labs(fill = paste0("HSI:")) +
      labs(title = paste0("Copepod functional group ", fg, ", (", fg_text,"), #spp: ",number_spp,", in month ", name_month)) +
      theme(plot.title = element_text(size=8)) + theme(legend.text=element_text(size=5)) +
      theme(legend.title = element_text(size=8))
    
    #map for the na case
    # g1 <- ggplot() + geom_raster(aes(x = x.2, y = y, fill = relfrq), data = df) +
    #   # play around with palettes. You can also use scale_fill_distiller() and
    #   geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
    #   coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    #                            panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
    #   scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
    #   scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    #   scale_fill_viridis(na.value = "white") + labs(fill = paste0("Relative frequency:")) +
    #   labs(title = paste0("Copepod no group assigned, #spp: ",number_spp,", in month ", name_month)) +
    #   theme(plot.title = element_text(size=8)) + theme(legend.text=element_text(size=5)) +
    #   theme(legend.title = element_text(size=8))
    # 
    #setwd(paste0(wd_fgv2,sdm,"/",fg))
    #ggsave(plot = g1, filename = paste0("HSI_copepod_functional_group_",fg,"_",name_month,"_",sdm,".jpg"),dpi = 300, width = 7, height = 4)
    setwd(wd_script)
    
    
    
    return(df_fg)
    
  }, mc.cores = 1
  )
  # -------------------------------------------------------------------------------------------------
  #names
  fg_text <- fg_desc[fg_desc$FG_f11 == fg,"FG_desc"][1]
  
  #######################
  setwd(wd_ct_GAM)
  comm <- read.table(file = files[1] , sep = "\t", h = T)
  #comm <- get(load(f)) # dim(comm) ; colnames(comm)
  # Change colnames (remove brackets or points in the species names if necessary)
  colnames(comm) <-  gsub("(\\.|\\.)", "", colnames(comm))
  #drop species 
  comm <- comm %>% dplyr::select(-c(TSS_means_omit_GAM$species))
  specieswithfg <- unique(species_fg_df$species) 
  # find the common species between those modelled (colnames in the monthly or annual community table and vector of names in the traits table
  common.spp <- intersect(specieswithfg,colnames(comm[4:(length(comm)-1)])) 
  number_spp <- length(intersect(species_fg_df[species_fg_df$FG_f11 == fg,]$species,common.spp))
  #######################
  # -------------------------------------------------------------------------------------------------
  
  table_df_fg <- bind_rows(res) #; rm(res) ; gc()
  # head(table.med.size) ; dim(table.med.size) ; summary(table.med.size)
  table_df_fg_an <- data.frame(table_df_fg %>% group_by(cell_id) %>% 
                             summarize(x = unique(x), y = unique(y), relfrq = mean(relfrq, na.rm = T)))
  df2 <- table_df_fg_an
  df2$x.2 <- df2$x
  df2[df2$x < 0 ,"x.2"] <- (df2[df2$x < 0 ,"x"]) + 360
  #Here I rename the month to filler so that it works with changing names
  colnames(df2) <- c("cell_id","x","y","relfrq","x.2")

  g2 <- ggplot() + geom_raster(aes(x = x.2, y = y, fill = relfrq), data = df2) +
    # play around with palettes. You can also use scale_fill_distiller() and
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                             panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
    theme(legend.text=element_text(size=6), legend.title = element_text(size=8)) +
    scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    scale_fill_viridis(na.value = "white") + labs(fill = paste0("HSI:")) +
    labs(title = paste0("Copepod functional group ", fg, ", (", fg_text,"), #spp: ",number_spp,", annual mean")) +
    theme(plot.title = element_text(size=8)) + theme(legend.text=element_text(size=5)) +
    theme(legend.title = element_text(size=8))
   
   
  
  # g2 <- ggplot() + geom_raster(aes(x = x.2, y = y, fill = relfrq), data = df2) +
  #   # play around with palettes. You can also use scale_fill_distiller() and
  #   geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
  #   coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
  #                            panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
  #   scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
  #   scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
  #   scale_fill_viridis(na.value = "white") + labs(fill = paste0("Relative frequency:")) +
  #   labs(title = paste0("Copepod no group assigned, #spp: ",14,", annual")) +
  #   theme(plot.title = element_text(size=8)) + theme(legend.text=element_text(size=5)) +
  #   theme(legend.title = element_text(size=8))
  
  
  # setwd(paste0(wd_fgv2,sdm,"/",fg))
  # ggsave(plot = g2, filename = paste0("HSI_copepod_functional_group_",fg,"_annual_mean_",sdm,".jpg"),dpi = 300, width = 7, height = 4)
  # setwd(wd_fd_data)
  # save(table_df_fg_an, file = paste0("table_df_fg_an_GAM_hsi_",fg)) #save dataframe for further use
  
}
