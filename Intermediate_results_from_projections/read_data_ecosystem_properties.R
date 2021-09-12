
### 11/01/2021: R script to read and modify, when needed, the annual field of the selected marine ecosystem properties, whicb serve as proxies for marine ecosystem services provision

# - Primary production and CO2 removal from the atmosphere via the Biological Carbon Pump --> DeVries & Weber (2017)
# - Plankton size index (slope from a size spectrum derived from satellite products), shows where larger animals are/ could emerge --> Kostadinov et al. (2009,2016)
# - Small (< 30cm) pelagic fishes reported and unreported fish catches --> estimate of fisheries effort (NOT FISH BIOMASS or 'NATURAL' FISH PRODUCTION) and amount of planktivorous fishes caught by various fishing gears. Proxy of food provision for human societies --> Watson 2017
# - Marine biodiversity estimates: species richness estimates for various animal groups (tunas, sharks, cetaceans etc.) --> Tittensor et al. (2010) ; Marine Biodiversity is considered as an ecosysem service througth its links to recreational activities (tourism) but also how it influences human cultures and natural ecosystem functioning

#UPDATED: JONAS, 25.01.2021
Sys.setenv(LANGUAGE= 'en')

library("raster") #needed
library("tidyverse") #needed
library("maps") #needed
library("ncdf4") #needed


### -----------------------------------------------------------------------------------------------------------------------------------------------------------------
wd_script <- ("/data/")
wd_fish <- ("/data/Fish catches Watson 2017/")
wd_part <- ("/data/Carbon DeVries&Weber 2017/")
wd_npp <- ("/data/NPP ModisAqua/")
wd_planksize <- ("/data/Plankton size Kostadinov&al._2009/")
wd_biodiv <- ("/data/Marine Biodiversity Tittensor&al._2010/")
wd_dat <- ("/data/intermediate_results_from_projections/hsi/")

world2 <- map_data("world2") # world coastline for maps 
world1 <- map_data("world")

dat <- world2[c("long","lat")]
colnames(dat) <- c("x2","y") # these coords are not in 0.5 steps but smaller and somewhat in disorder
dat$cell_id <- factor(paste(round(dat$x2,1), round(dat$y,1), sep = "_"))
### -----------------------------------------------------------------------------------------------------------------------------------------------------------------
### 1°) Carbon related variables: net primary production (NPP), Flux of POC (at 100m depth or euphotic zone depth), e ratio (biological carbon pump efficiency) etc.
setwd(wd_part)

nc <- nc_open("Cexp_deVries_2017.nc")

setwd(wd_part)
lat <- raster::stack("Cexp_deVries_2017.nc", varname = "LAT")
lon <- raster::stack("Cexp_deVries_2017.nc", varname = "LON")
y <- as.data.frame(lat, xy = T) 
x <- as.data.frame(lon, xy = T)


vars <- c("NPP","FPOCex","FPOC100m","POCfast","POCslow","POCflux")
# v <- "FPOC100m"

### !!! You might need to adjust some bits of the lapply() below
clims <- lapply(vars, function(v) {
				# Get data from the nc file
				ras <- raster::stack("Cexp_deVries_2017.nc", varname = v)
				# Turn into ddf and provide coordinates
				d <- as.data.frame(ras, xy = T)
				d$x <- x[,3] 
				d$y <- y[,3]
				
				d$cell_id <- factor(paste(round(d$x,1), round(d$y,1), sep = "_"))
				d$avg <- rowMeans(as.matrix(d[,c(3:14)]), na.rm = T) 
				# Change colnames 
				colnames(d)[c(3:14)] <- paste("v",c(1:12), sep = "")
				
				colnames(d)[16] <- v
				if(v == "FPOCex") { 
					return(d[,c(1,2,15,16)]) 
				} else {
					return( d[,c(1,2,16)] ) 
				} # eo if else loop
		} # eo FUN
) # eo lapply
# Cbind
# str(clims)
poc <- do.call(cbind, clims)

poc <- poc[c("x","y","cell_id","NPP","FPOCex","FPOC100m","POCslow","POCflux")]
#colnames(poc) <- c("x.2","y","cell_id","NPP","FPOCex","FPOC100m","POCslow","POCflux")
poc$eratio <- (poc$FPOC100m) / (poc$NPP)
summary(poc)

# Remove values < 0
poc <- poc[which(poc$NPP > 300 & poc$FPOCex > 300 & poc$eratio >= 0 & poc$eratio <= 0.4),]

# Next, convert to raster and combine with div indices
spg <- poc[,c("x","y","FPOC100m")]
#colnames(spg)[1] <- "x.2"

#rename for coords
colnames(spg)[1] <- "x"
colnames(spg)[2] <- "y"

coordinates(spg) <- ~ x + y
gridded(spg) <- TRUE
ras <- raster(spg)
crs(ras) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

summary(poc$FPOC100m)

### ---------------------------------------------------------------------
### Alternate version for NPP: source satellite products from ModisAqua. Might need to play around with the spatial grid so it matches the 1°x1° cell grid we're using.
setwd(wd_npp)
nc <- nc_open("NPP_ModisAQUA_cbpm_2002-2015_mon.nc") # ; nc
npp <- raster::stack("NPP_ModisAQUA_cbpm_2002-2015_mon.nc", varname = "npp")
#plot(log1p(npp))

lat <- raster::stack("NPP_ModisAQUA_cbpm_2002-2015_mon.nc", varname = "Lon")
lon <- raster::stack("NPP_ModisAQUA_cbpm_2002-2015_mon.nc", varname = "Lat")
y <- as.data.frame(lat, xy = T)
x <- as.data.frame(lon, xy = T)

dNPP <- as.data.frame(npp, xy = T)
# dim(dNPP) ; head(dNPP) ; colnames(dNPP)
colnames(dNPP)[c(3:14)] <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
dNPP$x <- y$Lon
dNPP$y <- x$Lat
summary(dNPP)
# Convert to the 1x1 resolution
dNPP$x_1d <- round(dNPP$x, .1)
dNPP$y_1d <- round(dNPP$y, .1)

unique(dNPP$x_1d)  
# And flip x coordinates
dNPP$x2 <- dNPP$x_1d 
dNPP[dNPP$x_1d < 0 ,"x2"] <- (dNPP[dNPP$x_1d < 0 ,"x_1d"]) + 360
#unique(div$x)[order(unique(div$x))]
unique(dNPP$x2)[order(unique(dNPP$x2))]
# Need to add 0.5 to match div
dNPP$x2 <- (dNPP$x2)+0.5
#unique(div$y)[order(unique(div$y))] 
unique(dNPP$y_1d)[order(unique(dNPP$y_1d))] 
### Add +0.5
dNPP$y2 <- (dNPP$y_1d)+0.5
#unique(div$y)[order(unique(div$y))]
#unique(dNPP$y2)[order(unique(dNPP$y2))] 

# add an id and compute monthly means and then annual climatology
dNPP$id <- factor(paste(dNPP$x2, dNPP$y2, sep = "_"))
clims <- data.frame(dNPP %>% group_by(id) %>% summarize(x = unique(x2), y = unique(y2),
                                                        Jan = mean(Jan,na.rm = T), Feb = mean(Feb,na.rm = T), Mar = mean(Mar,na.rm = T),
                                                        Apr = mean(Apr,na.rm = T), May = mean(May,na.rm = T), Jun = mean(Jun,na.rm = T),
                                                        Jul = mean(Jul,na.rm = T), Aug = mean(Aug,na.rm = T), Sep = mean(Sep,na.rm = T),
                                                        Oct = mean(Oct,na.rm = T), Nov = mean(Nov,na.rm = T), Dec = mean(Dec,na.rm = T) )
) # eo ddf
summary(clims) ; dim(clims)

### Compute annual clim with rowMeans
clims$Annual <- rowMeans(as.matrix(clims[,c(4:15)]), na.rm = T)
summary(clims$Annual)
colnames(clims)[1] <- "cell_id"
colnames(clims)[16] <- "Annual_NPP_v2"
tmp_df <- clims[c("cell_id","x","y","Annual_NPP_v2")]
summary(tmp_df)

### extract fpoc100m
tmp_df$FPOC100m <- raster::extract(x = ras, y = tmp_df[,c("x","y")])
summary(tmp_df$FPOC100m)
#Why are here suddenly 64800 - 26920NA vlaues 
### Same with Cexp below euphotic zone
spg <- poc[,c("x","y","FPOCex")]
colnames(spg)[1] <- "x"
coordinates(spg) <- ~ x + y
gridded(spg) <- TRUE
ras <- raster(spg)
crs(ras) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
### extract and join with ddf
tmp_df$FPOCex <- raster::extract(x = ras, y = tmp_df[,c("x","y")])
summary(tmp_df$FPOCex)

### And NPP
spg <- poc[,c("x","y","NPP")]
colnames(spg)[1] <- "x"
coordinates(spg) <- ~ x + y
gridded(spg) <- TRUE
ras <- raster(spg)
crs(ras) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
### extract and join with ddf
tmp_df$NPP <- raster::extract(x = ras, y = tmp_df[,c("x","y")])
summary(tmp_df$NPP)

### Compute export efficency (FPOCex/ NPP)
tmp_df$e <- (tmp_df$FPOCex) / (tmp_df$NPP)
summary(tmp_df$e) # Remove values > 1 (not possible)
tmp_df <- tmp_df[tmp_df$e <= 0.4,]
summary(tmp_df)

### NPP, FPOCex and FPOC100m are expressed in mmolC/m^2/yr --> convert to the standard mgC/m2.day
### Simply multiply by 12.0107 (molar mass of Carbon) and divide by 365
tmp_df$FPOCex <- (tmp_df$FPOCex)*(12.0107/365)
tmp_df$NPP <- (tmp_df$NPP)*(12.0107/365)
tmp_df$FPOC100m <- (tmp_df$FPOC100m)*(12.0107/365)
summary(tmp_df)

data_particles <- tmp_df
setwd(wd_script)
save(data_particles, file = "data_particles.RData")
### Map !

g1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = sqrt(Annual_NPP_v2)), data = clims) +
 	scale_fill_viridis(name = "Annual NPP\nsqrt(mgC/m2/yr)", na.value = "white") +
 	geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(0,60,120,180,240,300,360),
               labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
 	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
 

g2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = e ), data = data_particles) +
  # play around with palettes. You can also use scale_fill_distiller() and
  geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
  coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                           panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
  theme(legend.text=element_text(size=6), legend.title = element_text(size=8)) +
  scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
  scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) + 
  scale_fill_viridis(na.value = "white")
       
# setwd(wd_script)
# ggsave(plot = g1, filename = paste0("Annual_NPP_ModisAqua.jpg"),dpi = 300, width = 7, height = 4)

### -------------------------------------------------------------------

### 2°) Plankton size index (shows where larger particles occur)
### Simple read.csv
setwd(wd_planksize)
plank_dat <- read.csv2("table_annual_PSDslope_Kostadinov_1d_29_01_20.txt", sep = "\t", h = T)

colnames(plank_dat) <- c("cell_id","x","y","slope_psi","x.2")
plank_dat$x.2 <- plank_dat$x.2 + 0.5
plank_dat$y <- plank_dat$y + 0.5
plank_dat$cell_id <- factor(paste(round(plank_dat$x.2,1), round(plank_dat$y,1), sep = "_"))
plank_dat$slope_psi <- as.numeric(plank_dat$slope_psi)

plank_dat_temp <- plank_dat[c("cell_id", "slope_psi")]
summary(plank_dat_temp) #around 35k

setwd(wd_script)
data_plank <- plank_dat[c("cell_id","x.2","y","slope_psi")]
colnames(data_plank)[2] <- "x"
save(data_plank, file = "data_plankton.RData")



lglgl <- merge(tmp_df,plank_dat_temp, by = "cell_id")
#lglgl <- lglgl[is.na(lglgl$cell_id) == F,]
### -------------------------------------------------------------------

### 3°) Small pelagic fish catches. Mean annual climatology derived from Watson 2017. 
### Simple get(load(filename.Rdata)) in R.
setwd(wd_fish)
load("clim_pelagic%3C30cm_catches_1990-2019_1d.Rdata")

df_fish <- clim2save
rm(clim2save);gc();

df_fish$x.2 <- df_fish$x 
df_fish[df_fish$x < 0 ,"x.2"] <- (df_fish[df_fish$x < 0 ,"x"]) + 360
colnames(df_fish)[1] <- "cell_id"
df_fish$cell_id <- factor(paste(round(df_fish$x.2,1), round(df_fish$y,1), sep = "_"))
data_fish <- df_fish[df_fish$logged >= 1,]
data_fish <- data_fish[c("cell_id", "x.2", "y", "logged")]
colnames(data_fish)[2] <- "x"

setwd(wd_script)
save(data_fish, file = "data_fish.RData")

temp <- data_fish[c("cell_id","logged")]
summary(temp)#10523 values


ggplot() + geom_tile(aes(x = x, y = y, fill = slope_psi), data = data_plank) +
  scale_fill_viridis(name = "") +
  geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
  coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
                                        labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
  scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
                     labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
  theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
        panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )

dff <- merge(lglgl,temp, by = "cell_id" )
### -------------------------------------------------------------------

### 4°) Marine biodiversity fields --> how to read shapefiles and convert to data.frames
### BEWARE: THESE ESTIMATES ARE AT A MUCH COARSER RESOLUTION, so you'll need to decrease the 1°x1° resolution of your variables for proper comparison
setwd(paste0(wd_biodiv,"01_Data"))
shape <- readOGR(dsn = "WCMC-019-PatternsBiodiversity2010-AcrossTaxa.shp")
# Shape2 for all taxa 
shape2 <- readOGR(dsn = "WCMC-019-PatternsBiodiversity2010-IndivTaxa.shp")
class(shape2) ; str(shape2) ; head(shape2)
# Convert to raster: use a standard 1x1 grid from WOA
setwd(wd_script)
ras <- raster::raster("woa13_decav_t01_01_reference_grd.nc")
ras ; plot(ras)
# ?rasterize the 3 fields
AllNorm <- rasterize(x = shape, y = ras, field = "AllNorm")
AllNorm
plot(AllNorm)

OceanNorm <- rasterize(x = shape, y = ras, field = "OceanNorm")
OceanNorm
plot(OceanNorm)

### Use euphausiids and forams for validation of zooplankton diversity fields? 
krill <- rasterize(x = shape2, y = ras, field = "Euphausiid")
forams <- rasterize(x = shape2, y = ras, field = "ForamCK")
plot(krill)
plot(forams)

# Convert to ddf
biodiv <- as.data.frame(shape, xy = T)
class(biodiv) ; dim(biodiv); summary(biodiv)
colnames(biodiv)[c(2:3)] <- c("x","y")

# Make longitudes match div
biodiv$x  <- round(biodiv$x, .1)
biodiv$y  <- round(biodiv$y, .1)

biodiv$x.2 <- biodiv$x 
biodiv[biodiv$x < 0 ,"x.2"] <- (biodiv[biodiv$x < 0 ,"x"]) + 360

biodiv$x.2 <- biodiv$x.2 + 0.5
biodiv$y <- biodiv$y + 0.5

biodiv$cell_id <- factor(paste(round(biodiv$x.2,1), round(biodiv$y,1), sep = "_"))
temp_biodiv <- biodiv[c("cell_id", "x.2", "y", "AllTaxa", "AllNorm", "CoastNorm", "OceanNorm")]
summary(temp_biodiv)#630 values
data_biodiv <- temp_biodiv

setwd(wd_script)
save(data_biodiv, file = "data_biodiv.RData")

dff2 <- merge(lglgl,temp_biodiv, by = "cell_id" )


# Test map
ggplot() + geom_tile(aes(x = x.2, y = y, fill = OceanNorm), data = data_biodiv) +
      scale_fill_viridis(name = "") +
      geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
      coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
                labels = c("0°W","60°W","120°W","180°W","-120°W","-60°W","0°W"), expand = c(0,0)) +
      scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
                labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
        theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
          panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )





