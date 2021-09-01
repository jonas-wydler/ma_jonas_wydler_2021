
##### ETHZ, Script from Fabio Benedetti, extended by Jonas Wydler on the 18.12.2020
#	 - Re-fitting the env predictors (updated version from 18/09/2020) to the species occurrences
#updated:1.9.2021 for legacy after leaving the group
Sys.setenv(LANGUAGE= 'en')
# -------------------------------------------------------------------------------------------
library("raster")
library("rgeos")
library("raster")
library("maptools")
library("rgdal")
library("tidyr")
library("ncdf4")
library("RColorBrewer") # nice color palettes (discrete and continuous)
library("viridis") # even better palettes for continuous data like SST etc.
library("cmocean")
# -------------------------------------------------------------------------------------------
#define working dirs
wd_data <- ("") #dir where the data

world2 <- map_data("world2") # world coastline for maps 
world1 <- map_data("world")

# 
# dat_test <- read.table("data_thinned_Centropages_furcatus_28_22_19_XX_XX_XX.txt", sep = ";", h = T)
# str(dat_test)
# 
# clims_test <- clims[[4]]
# clims_test_df <- as.data.frame(clims_test)
# 
# summary(clims_test_df$x.1)
# 
# #to plot the maps centered correctly
# clims_test_df$x.1_2 <- clims_test_df$x.1
# clims_test_df[clims_test_df$x.1 < 0 ,"x.1_2"] <- (clims_test_df[clims_test_df$x.1 < 0 ,"x.1"]) + 360
# 
# #world2$long2 <- world2$long[world2$long < 0 ,"long2"] <- (world2[world2$long < 0 ,"long"]) + 360
# 
# ### Example of a map using world2
# ggplot() + geom_raster(aes(x = x.1_2, y = y.1, fill = Sistar), data = clims_test_df) +
#     scale_fill_viridis(name = "") + # play around with palettes. You can also use scale_fill_distiller() and
#     geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
#     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
#     scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
#     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) 
# 
# WD <- getwd()

# --------------------------------------------------------------------------------------------------------------------------------

### Define the directory with the env climatologies
setwd(paste0("wd_data/global_monthly_clims_1d/"))
cl <- dir()[grep("18_09_20", dir())] ; cl
# Load each clim, turn into raster stack and store them in a list
# test with c = cl[[1]]
clims <- mclapply(cl, function(c) {
			# read table
			message(paste("Reading ",c, sep = ""))
			dat <- read.table(c, sep = ";", h = T)
			# convert to stack
			coordinates(dat) <- ~ x + y
			gridded(dat) <- TRUE
			# coerce to raster
			stak <- stack(dat)
			crs(stak) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
            # drop: Chla, logChla, PAR, MLD, MLPAR (newer fields available), indices are: 2,6,7,11,12,15 ### ADAPT TO YOUR LIKING HERE
            stk2 <- dropLayer(stak, i = c(2,6,7,11,12,15))
			# return
			return(stk2)
		}, mc.cores = 1 #fabio put 15 here, doesnt work on windows.
) # eo lapply
length(clims) # 12 elements, one per month, containing 24 fields 
for(i in c(1:12)) {
		message( paste0(names(clims[[i]]),sep = "__")) #had to change to paste0, not sure why
		message( paste("  ", sep = "") )
} # eo for loop
# --------------------------------------------------------------------------------------------------------------------------------
# #Chech for Environmental variable collinearity
# cilms_selection_df <- clims_test_df[,c(1,2,5:28)]
# cor(cilms_selection_df, use= "everything", method = "spearman")
# 
# cor(cilms_selection_df, use="complete.obs", method="kendall")
# --------------------------------------------------------------------------------------------------------------------------------
### In a mclapply, apply a function to fit the monthly env data to the occurrences 

### Got to the dir where species tables are stored
setwd(paste0("wd_data/data_species_v3.2v5.1_thinned/"))

files <- dir()[grep("data_",dir())]
files
#f <- files[101] # for testing function below before sending it to kryo with R CMD BATCH 

env.fitter <- function(f = files) {
    
            setwd(paste0("wd_data/data_species_v3.2v5.1_thinned/"))
            message(paste("Fitting the env data to ",f, sep = ""))
            message(paste("", sep = ""))
            message(paste("", sep = ""))
            data <- read.table(f, sep = "\t", h = T) # str(data) ; colnames(data)
            
            # Fit the env data from 'clims' first
    		env_annual <- clims[[1]]
    		# To extract the 1d coordinates corresponding to each occurrence
    		cells <- cellFromXY(subset(x = env_annual, subset = 1), as.matrix(data[,c("x","y")]))
    		coords <- data.frame(xyFromCell(subset(x = env_annual, subset = 1), cells),value = raster::extract(subset(x = env_annual, subset = 1), cells) ) 
    		
    		if( nrow(coords) == nrow(data) ) {
    			data$xbin_1d <- coords$x
    			data$ybin_1d <- coords$y
    		} # eo if else loop
            # summary(data)
    		data$Bathy <- raster::extract(x = subset(x = env_annual, subset = 'Bathy'), y = data[,c("xbin_1d","ybin_1d")], method = 'bilinear')
    		data$dSST <- raster::extract(x = subset(x = env_annual, subset = 'dSST'), y = data[,c("xbin_1d","ybin_1d")], method = 'bilinear')#change deltaT to DSST
            
    		# Create empty vectors for each predictor (20 in total since you already provided Bathy & dSST)#
    		vars <- names(clims[[2]])[c(5:19)] # vars ### MAKE SURE THESE MAKE SENSE
    		data[vars] <- NA

    		### And fill those (according to month) in a for loop
            layers <- c(1:12)
    		for(l in layers) {
				
    				# Get the monthly clims 
    				#message(paste("Matching month ", l, sep = ""))
    				envdata <- clims[[l]]
				
    				data[which(data$month == l),"logNO3"] <- raster::extract(x = subset(x = envdata, subset = 'logNO3'), 
    				y = data[which(data$month == l),c("xbin_1d","ybin_1d")], method = 'bilinear')
    				data[which(data$month == l),"logPO4"] <- raster::extract(x = subset(x = envdata, subset = 'logPO4'), 
    				y = data[which(data$month == l),c("xbin_1d","ybin_1d")], method = 'bilinear')
    				data[which(data$month == l),"logSiO2"] <- raster::extract(x = subset(x = envdata, subset = 'logSiO2'), 
    				y = data[which(data$month == l),c("xbin_1d","ybin_1d")], method = 'bilinear')
    				data[which(data$month == l),"Nstar"] <- raster::extract(x = subset(x = envdata, subset = 'Nstar'), 
    				y = data[which(data$month == l),c("xbin_1d","ybin_1d")], method = 'bilinear')
    				data[which(data$month == l),"NO3"] <- raster::extract(x = subset(x = envdata, subset = 'NO3'), 
    				y = data[which(data$month == l),c("xbin_1d","ybin_1d")], method = 'bilinear')
    				data[which(data$month == l),"PO4"] <- raster::extract(x = subset(x = envdata, subset = 'PO4'), 
    				y = data[which(data$month == l),c("xbin_1d","ybin_1d")], method = 'bilinear')
    				data[which(data$month == l),"Pstar"] <- raster::extract(x = subset(x = envdata, subset = 'Pstar'), 
    				y = data[which(data$month == l),c("xbin_1d","ybin_1d")], method = 'bilinear')
    				data[which(data$month == l),"SiO2"] <- raster::extract(x = subset(x = envdata, subset = 'SiO2'), 
    				y = data[which(data$month == l),c("xbin_1d","ybin_1d")], method = 'bilinear')
    				data[which(data$month == l),"Sistar"] <- raster::extract(x = subset(x = envdata, subset = 'Sistar'), 
    				y = data[which(data$month == l),c("xbin_1d","ybin_1d")], method = 'bilinear')
    				data[which(data$month == l),"SLA"] <- raster::extract(x = subset(x = envdata, subset = 'SLA'), 
    				y = data[which(data$month == l),c("xbin_1d","ybin_1d")], method = 'bilinear')
    				data[which(data$month == l),"SSS"] <- raster::extract(x = subset(x = envdata, subset = 'SSS'), 
    				y = data[which(data$month == l),c("xbin_1d","ybin_1d")], method = 'bilinear')
    				data[which(data$month == l),"SST"] <- raster::extract(x = subset(x = envdata, subset = 'SST'), 
    				y = data[which(data$month == l),c("xbin_1d","ybin_1d")], method = 'bilinear')
    				data[which(data$month == l),"Wind"] <- raster::extract(x = subset(x = envdata, subset = 'Wind'), 
    				y = data[which(data$month == l),c("xbin_1d","ybin_1d")], method = 'bilinear')
    				data[which(data$month == l),"dO2"] <- raster::extract(x = subset(x = envdata, subset = 'dO2'), 
    				y = data[which(data$month == l),c("xbin_1d","ybin_1d")], method = 'bilinear')
    				data[which(data$month == l),"EKE"] <- raster::extract(x = subset(x = envdata, subset = 'EKE'), 
    				y = data[which(data$month == l),c("xbin_1d","ybin_1d")], method = 'bilinear')
    				data[which(data$month == l),"logEKE"] <- raster::extract(x = subset(x = envdata, subset = 'logEKE'), 
    				y = data[which(data$month == l),c("xbin_1d","ybin_1d")], method = 'bilinear')

    		} # eo for loop - common WOA13 products
            
            # Fit the env data from the newer satellite products (Chlorophyll, PAR, MLPAR, MLD)
    		var2 <- c("Chla","logChla","MLD","PAR","MLPAR")
    		data[var2] <- NA
            
            for(v in var2) {
                    
                    #v <- 'logChla'
                    setwd(paste0("wd_data/global_monthly_clims_1d/"))
                    filename <- paste("clim_month",v,"27_11_19.txt", sep = "_")
    			    #message(paste("Reading ",filename, sep = ""))
    			    clim <- read.table(filename, sep = "\t", h = T)
    			    # convert to stack
    			    coordinates(clim) <- ~ x + y ;  gridded(clim) <- TRUE
    			    # coerce to raster
    			    stak <- stack(clim)
    			    crs(stak) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
                    # names(stak) ; plot(stak)
                    #plot(stak, band = 2)
                    #points(data[,c("x","y")], pch=16)
                    ### Need to rotate x coordinates for this actually ^^
                    data$x2bin_1d <- data$xbin_1d
                    data[data$xbin_1d < 0 ,"x2bin_1d"] <- (data[data$xbin_1d < 0 ,"xbin_1d"]) + 360
                                    
				    data[which(data$month == 1),v] <- raster::extract(x = subset(x = stak, subset = 'Jan'), 
				    y = data[which(data$month == 1),c("x2bin_1d","ybin_1d")], method = 'bilinear')
				    #
				    data[which(data$month == 2),v] <- raster::extract(x = subset(x = stak, subset = 'Feb'), 
				    y = data[which(data$month == 2),c("x2bin_1d","ybin_1d")], method = 'bilinear')
                    #
				    data[which(data$month == 3),v] <- raster::extract(x = subset(x = stak, subset = 'Mar'), 
				    y = data[which(data$month == 3),c("x2bin_1d","ybin_1d")], method = 'bilinear')
                    #
				    data[which(data$month == 4),v] <- raster::extract(x = subset(x = stak, subset = 'Apr'), 
				    y = data[which(data$month == 4),c("x2bin_1d","ybin_1d")], method = 'bilinear')
                    #
				    data[which(data$month == 5),v] <- raster::extract(x = subset(x = stak, subset = 'May'), 
				    y = data[which(data$month == 5),c("x2bin_1d","ybin_1d")], method = 'bilinear')
                    #
				    data[which(data$month == 6),v] <- raster::extract(x = subset(x = stak, subset = 'Jun'), 
				    y = data[which(data$month == 6),c("x2bin_1d","ybin_1d")], method = 'bilinear')
                    #
				    data[which(data$month == 7),v] <- raster::extract(x = subset(x = stak, subset = 'Jul'), 
				    y = data[which(data$month == 7),c("x2bin_1d","ybin_1d")], method = 'bilinear')
                    #
				    data[which(data$month == 8),v] <- raster::extract(x = subset(x = stak, subset = 'Aug'), 
				    y = data[which(data$month == 8),c("x2bin_1d","ybin_1d")], method = 'bilinear')
                    #
				    data[which(data$month == 9),v] <- raster::extract(x = subset(x = stak, subset = 'Sep'), 
				    y = data[which(data$month == 9),c("x2bin_1d","ybin_1d")], method = 'bilinear')
                    #
				    data[which(data$month == 10),v] <- raster::extract(x = subset(x = stak, subset = 'Oct'), 
				    y = data[which(data$month == 10),c("x2bin_1d","ybin_1d")], method = 'bilinear')
                    #
				    data[which(data$month == 11),v] <- raster::extract(x = subset(x = stak, subset = 'Nov'), 
				    y = data[which(data$month == 11),c("x2bin_1d","ybin_1d")], method = 'bilinear')
                    #
				    data[which(data$month == 12),v] <- raster::extract(x = subset(x = stak, subset = 'Dec'), 
				    y = data[which(data$month == 12),c("x2bin_1d","ybin_1d")], method = 'bilinear')	
                
            } # eo for loop - satellite products
            
            # Drop rotated x coord
            data <- subset(data, select = -c(x2bin_1d))
            
            # And fit new pCO2 product by P. Landschüzter (ESSD)
            pCO2 <- raster::stack("MPI-ULB-SOM_FFN_clim.nc")
            pCO2[pCO2 < 0] <- 50 # replacing wrong values
            # Disaggregate resolution from 1/4 to 1°
            pCO2v2 <- aggregate(pCO2, fact = 4) # pCO2v2
            rm(pCO2)
            
            data["pCO2"] <- NA
            
    		for(l in layers) {
    				data[which(data$month == l),"pCO2"] <- raster::extract(x = subset(x = pCO2v2, subset = l), 
    				y = data[which(data$month == l),c("xbin_1d","ybin_1d")], method = 'bilinear')
    		} # eo for loop - pCO2		
            
            # Print new table in approprate dir
            #setwd("/net/kryo/work/wydlerj/clim_and_species")
            setwd(paste0("wd_data/data_species_v3.2v5.1_thinned/"))
            # Change filename to save the table under
            filename <- str_replace_all(f,".txt","")
            filename <- str_replace_all(filename,"occurrence","fitted")
            filename <- paste(filename,"_fitted.txt", sep = "")
            
            write.table(x = data, filename, sep = ";")
    
    
} # eo env.fitter ---------------------------------

require("parallel")
mclapply(X = files, FUN = env.fitter, mc.cores = 32) # mc.cores speciifyies the nb of CPUs to run the fun in parallel on.
# On kryo, remember max CPU is = 48. I usually use 35 if it's not too crowded


