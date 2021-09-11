# -------------------------------------------------------------------------------------------
#ETHZ, script by Fabio Benedetti, adapted by Jonas Wydler.
# - Defining background data that is needed for SDMs.

#Updated for legacy on 1.09.2021
# -------------------------------------------------------------------------------------------

library("rgeos")
library("raster")
library("maptools")
library("rgdal")
library("tidyverse")
library("stringr")
library("reshape2")
library("geosphere")
library("ncdf4")
library("classInt")
library("parallel")

Sys.setenv(LANGUAGE= 'en')
world2 <- map_data("world2")
world <- map_data("world")

#load data that specifies the FG for each species
load("species_FG")
colnames(species_FG) <- c("species","FG_f11")
# -------------------------------------------------------------------------------------------

### Select dir where the combined clim/species tables are

setwd(paste0("wd_data/combined_species_clims/"))

files <- dir()[grep("data_thinned",dir())] # vector of filenames (species table with occurences and associated cliamtologies)

# We had to correct erroneous data for Calanus_finmarchicus. The issue was that there were wrong presences on the southeren hemisphere. 
# files_cf <- dir()[grep("data_thinned_Calanus_finmarchicus_28_22_19.txt",dir())] # vector of filenames (species table with occurrences and associated clims)
# cf <- read.table(files_cf, h = T, sep = "\t")
# cf_new <- subset(cf, y > 0)
# write.table(cf_new,"data_thinned_Calanus_finmarchicus_28_22_19.txt",sep="\t")

### Rbind all
# f <- files[5]
res <- mclapply(files, function(f) {
            message(paste("Retrieving ",f, sep = ""))
            message(paste("", sep = ""))
            data <- read.table(f, h = T, sep = ";")
            return(data)
        }, mc.cores = 1
) # eo mclapply
# Rbind
table <- bind_rows(res); rm(res); gc()

head(table); dim(table)
summary(table)

#Merge the combined clims/species table with information on FGs
table <- merge(table,species_FG, by  = "species", all=T)
# Set dir where you'll print the new species tables containing both the occurrences and the background data
setwd("")

# Specify the parameter(s) to stratify the sampled environment by 
vec.strat <- c("SST")

# Remove the brackets in some of the species names (e.g. Acartia species)
table$species <- gsub("\\(|\\)", "", table$species)
species <- unique(table$species)
species   # quickly check if these are correct


write.table(species,"speciestest.txt")
### Sites id
table$id2 <- paste(table$x, table$y, table$month, sep = "_") # 1 site = monthly-resolved 1°x1° grid cell (e.g. -0.5_0.5_1 or 65.5_-10.5_4 etc.)

speciestest <- species[c(1,2,3,4,5,6,14,17,46,66,78,99,101,220,333)]
### In a for loop across species (a mclapply is inside so we can't really parallel jobs across species datasets here, but it doesn't take long for most species)

#sp <- species[256]
for(sp in species2) {
	
		  message(paste("Drawing psAbs for ",sp, " ============================================ ", sep = ""))
		  FG <- unique(table[table$species == sp, "FG_f11"]) 
		  
		  # Get species data
		  subset <- table[table$species == sp,]
          n <- nrow(subset)
          # Identify sites where the species if found as present
          sites.with.pres <- unique(subset$id)
		      #sites.with.pres <- unique(subset$id2)
	      ### If, species belongs to FG 11, then draw pseudo absences only in FG11 sites (because missing southern ocean, according to target group approach)
	      if(!is.na(FG)&(FG == 11)) {
				rich <- data.frame(table[!((table$id2 %in% sites.with.pres) & table$FG_f11 == 11),] %>% group_by(id2) %>% summarise(n = n(), group.div = length(unique(species))))
	        
	        # summary(rich)
	          threshold <- 10 # Threshold of minimum rich per monthly site, could be 5 or 15...from the maps you showed me, 10 looks reasonable     
	          cells2keep <- rich[rich$group.div >= threshold,"id2"] # length(cells2keep)
              ### Define the sites that will constitute the background to draw pseudo-absences from
              bckgrnd <- table[table$id2 %in% cells2keep,]
              
	      } else {
	          
    		  # Define the background data sites depending on their minimum nb of 
              # Compute n species per sites WITHOUT THE SITES WHERE THE SPECIES OF INTEREST HAS BEEN FOUND
    	      rich <- data.frame(table[!(table$id2 %in% sites.with.pres),] %>% group_by(id2) %>% summarise(n = n(), group.div = length(unique(species))))
              # summary(rich)
              threshold <- 10 # Threshold of minimum rich per monthly site, could be 5 or 15...from the maps you showed me, 10 looks reasonable     
              cells2keep <- rich[rich$group.div >= threshold,"id2"] # length(cells2keep)
         
              ### Define the sites that will constitute the background to draw pseudo-absences from
              bckgrnd <- table[table$id2 %in% cells2keep,]
              
	      }
	
		  ### Quick sanity check							
          intersect(unique(bckgrnd$id2), sites.with.pres) # should be 0
								
		  ### Specify range of SST into which the values fall to drive sampling of absences proportionally to the overall presences points
		  x_envir <- bckgrnd[,c(vec.strat)] # vector of SST in the bckgrnd data
		  # Split ranges into 9  strata
		  breaks <- 9
		  ### Create a matrix that divides range into 9 equal parts; with two variables we get a maximum of 81 strata
		  x_breaks <- classIntervals(na.omit(x_envir), breaks, style = "equal")
		  x_matrix <- cbind(x_breaks$brks[1:breaks], x_breaks$brks[2:(breaks+1)], ID = 1:breaks)
		  colnames(x_matrix) <- c("low","up","ID")	   
          
          # Define vector of length of total points of environmental variable
		  x_reclass <- c(1:length(x_envir))
		  				
		  # Allocate points from full data to one of the nine environmental strata per variable
		  for(i in 1:breaks) {	
		      x_reclass[which(x_envir >= x_matrix[i,"low"] & x_envir <= x_matrix[i,"up"] )] <- x_matrix[i,"ID"]	
		  } # eo for loop

		  ### Create an ID indicating the stratum (unique combination of variables) into which each point falls in full data-frame
		  bckgrnd$x_rcls <- x_reclass
								
		  print( paste0(sp,", ",FG," | n = ",n, " | drawing psAbs")) # eo print
								
		  ### Extract frequencies by which points/sites of the target group fall into environmental strata. 
		  # Then, derive the number of desired absences for the focal model species per stratum. 
		  x_rcls_freq <- data.frame( table(bckgrnd$x_rcls) / length(bckgrnd$x_rcls) ) 		
		  # Give name to column
		  colnames(x_rcls_freq)[1] <- "x_rcls" 
		  # Convert to numeric
		  x_rcls_freq$x_rcls <- as.numeric(as.character(x_rcls_freq$x_rcls)) 
		  ### Add desired background points to be produced per stratum: generally 10 x more absences than presences
		  x_rcls_freq$prop_abs <- (n*10)*x_rcls_freq$Freq
		  # To round desired absences to integer: adds column difference between smaller closest integer and desired number
		  x_rcls_freq$prop_abs_0 <- ( x_rcls_freq$prop_abs - floor(x_rcls_freq$prop_abs) )
		  # To add column with random number between 0 and 1 (with steps of 0.01)
		  x_rcls_freq$prob <- sample(seq(0, 1, 0.01), nrow(x_rcls_freq), replace = T)
		  # To add column with "1"
		  x_rcls_freq$absences <- 1
		  # Round up absences for random subset
		  x_rcls_freq$absences[which(x_rcls_freq$prop_abs_0 > x_rcls_freq$prob)] <- ceiling(x_rcls_freq$prop_abs[which(x_rcls_freq$prop_abs_0 > x_rcls_freq$prob)])
		  # Round absences down for random subset
		  x_rcls_freq$absences[which(x_rcls_freq$prop_abs_0 < x_rcls_freq$prob)] <- floor(x_rcls_freq$prop_abs[which(x_rcls_freq$prop_abs_0 < x_rcls_freq$prob)]) 
		  
          # Skip strata without presences
		  absence_groups <- x_rcls_freq[x_rcls_freq$absences > 0,] 
		  # Select backround data, here including the points/sites of the focal species ('overlapping background')
		  absence_table <- bckgrnd ; gc()
			
		  # Randomly select background pts
		  # nnn <- nrow(absence_groups) 
          
          nnn <- 9
          
		  require("parallel")#parallel
		  psAbs <- lapply(X = c(1:nnn), FUN = function(i) {
	
						# Select available absences within stratum in question
						message(paste(i, sep = ""))
						grp_abs_table <- absence_table[absence_table$x_rcls == absence_groups[i,"x_rcls"],]

						# Define the max nb of absences that can be drawn 
						absence_num <- ifelse(
								# Test if the number of desired background pts is bigger than the available background points
								absence_groups[i,"absences"] > nrow( absence_table[absence_table$x_rcls == absence_groups[i,"x_rcls"],]),		
								# if TRUE the potential points are insufficient - however, save the number of available points as absence_num
								nrow(absence_table[absence_table$x_rcls == absence_groups[i,"x_rcls"],]),		
								# ELSE: save the number of desired background points as absence_num   		
								absence_groups[i,"absences"]
						) # eo if else loop
						
						# Randomly sample the background points from the table containing all possible absences for the stratum in question
						sampled_grp_abs_table <- grp_abs_table[sample(1:nrow(grp_abs_table), size = absence_num),]
						rm(absence_num, grp_abs_table)
				
						return(sampled_grp_abs_table)
				
					} # eo fun
			
		   ) # eo lapply
			
		   ### Merge presences (obs = 1) with absences (obs = 0)
		   pseudoabs <- data.frame(do.call("rbind", psAbs), obs = 0)
		   occ_table <- rbind(data.frame(subset, obs = 1),  pseudoabs[,c(colnames(subset),"obs")])
				
           rm(psAbs, nnn) ; gc()
           
		   ### Create column with weights = 1; weights are associated with presences and absences for modelling
		   occ_table$weights <- 1
		   # Compute ratio of presences to absences
		   abs_ratio <- nrow(occ_table[occ_table$obs == 1,]) / nrow(occ_table[occ_table$obs == 0,])
			
		   # Add the ratio as weight for the psAbs
		   occ_table$weights[occ_table$obs == 0] <- abs_ratio # For observation that are absences we add the ratio
		   row.names(occ_table) <- c(1:nrow(occ_table)) # 
                
           # Final checks
           # dim(occ_table) ; colnames(occ_table) ; str(occ_table)
           # summary(occ_table)
           if( length(unique(occ_table$obs)) == 2 ) {
               
    		   ### Save the data to train some ENMs later
    		   setwd(wd.background)
    		   message(paste("Saving species dataset for ",sp, " ============================================ ", sep = ""))
    		   write.table(occ_table, paste("backgrounddata_", sp,sep = ""), sep = ";")
			
    		   ### Clean some stuff 
    		   rm(occ_table, abs_ratio, pseudoabs, absence_table, x_rcls_freq, x_reclass, x_matrix, x_envir)
    		   gc()
				
               setwd(wd2)
               
           } else {
               
               message(paste("  ", sep = ""))
               message(paste(" |||  ISSUE, NO PSEUDO-ABSENCES WERE GENERATED  ||| ", sep = ""))
               message(paste("  ", sep = ""))
               
           }
            	
} # eo for loop


### Making sanity checks for the various datasets created: look at proprotions of 1 and 0 for each file            
setwd(wd.background)

files <- dir(wd.background)[grep("background",dir())] # 
#files2 <-  dir(wd.background)[grep("backgrounddata_Calanus_finmarchicus",dir())]   
for(f in files) {
           
    f <- files[5]
    d <- read.table(f, sep = ";", h = T)
    avg <- round(mean(d$obs, na.rm = T),4)
            
    if(avg < 0.1) {
        
        message(paste("", sep = ""))
                
    } else {
                
        message(paste("Mean proportion of 1:0 for ",unique(d[d$obs == 1,"species"])," ||  ",avg, sep = ""))
        
        require("ggthemes")
        map <- ggplot() + geom_point(aes(x = x, y = y, colour = factor(obs)), data = d, alpha = .5) +
                     geom_polygon(aes(x = long, y = lat, group = group), data = world,
                         fill = "grey75", colour = "black", size = 0.2) +
                     coord_quickmap() + theme_map()
        # save maps showing the background data
                setwd("")#dir where you want the maps
               
                ggsave(plot = map, filename = paste("mapbackground",unique(d[d$obs == 1,]$species),".png"),dpi = 600, width = 16, height = 9)
                setwd(wd.background)
            }
            
} # eo for loop - f in files


