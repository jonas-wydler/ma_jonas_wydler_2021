#Rscript by Fabio Benedetti, adapted by Jonas Wydler for Master thesis
# - Creates univariate response curves for copepod species
#27.11.2020
#updated:08.09.2021
# --------------------------------------------------------------------------------------------------------------------------------

#library("raster")
#library("sp")
#library("reshape2")
library("tidyverse")#needed
library("biomod2") #needed
#library("R.devices")

Sys.setenv(LANGUAGE= 'en')

#data.wd.back <- ("wd_data/species_background_data/")
data.wd.back2 <- ("C:/Users/Jonas/ma/sdm/part2/kryodw/") #delete
data.wd.back <- ("C:/Users/Jonas/ma/sdm/background")
# --------------------------------------------------------------------------------------------------------------------------------

# Vector of SDMs
SDMs <- c("ANN","GLM","GAM")
#SDMs <- "GAM"
eval_runs <- c("RUN1","RUN2","RUN3","RUN4","RUN5","RUN6","RUN7","RUN8","RUN9","RUN10") 

### Set modelling options
myBiomodOption <- BIOMOD_ModelingOptions(
	
						GLM = list( type = 'quadratic',
    								interaction.level = 0,
   			 						myFormula = NULL,
    								test = 'AIC',
    								family = binomial("logit"),
    								mustart = 0.5,
    								control = glm.control(epsilon = 1e-08, maxit = 50, trace = FALSE)),

						GAM = list(algo = 'GAM_mgcv',
								type = 's_smoother',
								k = 5,
								interaction.level = 0,
								myFormula = NULL,
								family = binomial("logit"),
								method = 'GCV.Cp',
								optimizer = c('outer','newton'),
								select = FALSE,
								knots = NULL,
								paraPen = NULL,
								control = list(nthreads = 1, irls.reg = 0, epsilon = 1e-07
								, maxit = 200, trace = FALSE, mgcv.tol = 1e-07, mgcv.half = 15
								, rank.tol = 1.49011611938477e-08
								, nlm = list(ndigit=7, gradtol=1e-06, stepmax=2, steptol=1e-04, iterlim=200, check.analyticals=0)
								, optim = list(factr=1e+07)
								, newton = list(conv.tol=1e-06, maxNstep=5, maxSstep=2, maxHalf=30, use.svd=0)
								, outerPIsteps = 0, idLinksBases = TRUE, scalePenalty = TRUE, keepData = FALSE))
								
) # eo modelling options

### Set the working directories

setwd(data.wd.back) # change accordingly

files <- dir()[grep("backgrounddata",dir())] # vector of species files
# Go to the niche.modelling dir, where your models will be stored    
setwd(paste(data.wd.back2,"/","niche.modelling_resp_curves", sep = ""))
second.wd <- getwd()
    
### In a for loop across species files, built the GAM from a set of chosen common predictors (same for all, like in Benedetti et al., 2018a), and derive the univariate response curves and their summarizing statistics (aka niche traits, so center and width)

vars <- c("SST","logChl","logNO3","MLD","Sistar", "PAR", "Nstar", "logEKE") # choose the top 6-8 variables you want to uyse for comapring the species in env niche space (Benedetti et al., 2018)

f <- files[225]
#files2nd <- files[327:385]
for(f in f) {
	
        	    # Get the data
        	    setwd(data.wd.back)
        	    message(paste("Modelling ",f, " =========================================================", sep = ""))
        	    data <- read.table(f, h = T, sep = ";")
                g <- unique(data[data$obs == 1,"FG_f11"]) # change "FG_f11" by your FG column
                            
                ### Change data colnames for: logChla --> logChl
              colnames(data)[grep("logChla", colnames(data))] <- "logChl"
              colnames(data)[grep("Chla", colnames(data))] <- "Chl"
                            
              data2 <- na.omit( data[,c("species","x","y","obs",vars)] )        	    
        	    
        	    n <- nrow( data2[data2$obs == 1,] )
        	    rm(data2)
							
        		# If n >= 50, continue, otherwise got to next species
        		if( n >= 50) {
								
        	  			### Initialisation: data formatting
        				myRespName <- str_replace_all(unique(data[data$obs == 1,"species"]), "_", ".")
        				myRespName <- gsub("\\(|\\)", "", myRespName)       
                        if( grepl(pattern = '-', x = myRespName, fixed = T) ) {
                            myRespName <- str_replace_all(myRespName,'-','')   
                        }
        				# the presence/absences data for our species 
        				myResp <- as.numeric(data$obs)
        				
        			# the XY coordinates of species data
        				myRespXY <- data[,c("x","y")]        
                myExpl <- data[,vars]
        				# weights vector
                weights <- na.omit(data[,c(vars,"weights")])[,"weights"]
                
                #weights <- data$weights
                
        				# data formating
                        setwd(second.wd)
                                
        	  			myBiomodData <- BIOMOD_FormatingData(resp.var = myResp, expl.var = myExpl, resp.xy = myRespXY, resp.name = myRespName)
							
        	 			### Species Distribution Modelling
                #nulldev() 
        				myBiomodModelOut <- BIOMOD_Modeling(myBiomodData, models = SDMs, models.options = myBiomodOption, NbRunEval = length(eval_runs), 
        								   DataSplit = 80, Yweights = weights, VarImport = 0, models.eval.meth = c('TSS'), SaveObj = F, do.full.models = F
        				) # eo modelling
        				#complains that weights doent match obs, weird because in other script no problem
        				
        				#Error in which(sapply(x, is.factor)) : argument to 'which' is not logical, w
        		
        			
        				 #dev.off()					
  
        				### Get evaluation scores
        				scores <- data.frame(get_evaluations(myBiomodModelOut))
        				# Use lapply to summarize all scores 
        				scores.ls <- lapply(SDMs, function(sdm) {
        							# Retrieve all scores for the sdm 's' 
        							tss <- scores["TSS",grep(paste("Testing.data.",sdm, sep = ""), colnames(scores))]
        							tss_cutoffs <- scores["TSS",grep(paste("Cutoff.",sdm, sep = ""), colnames(scores))]
        							colnames(tss) <- paste(sdm, eval_runs, sep = "_")
        							tss <- t(tss)
        							tss_cutoffs <- t(tss_cutoffs)
        							colnames(tss_cutoffs) <- "Cutoff_TSS"
        							return(cbind(tss, tss_cutoffs))
        					} 
                        ) # eo lapply
        				scores.tbl <- data.frame(do.call(rbind, scores.ls))
                scores.tbl$species <- unique(data[data$obs == 1,"species"])
                scores.tbl$group <- g
        				rm(scores,scores.ls)
                        
                        ### Get response curves
                        # sdm <- "GAM"
        				        # sdm <- "GLM"
        				        # sdm <- "ANN"
                for(sdm in SDMs) {
                            
            				response.plot <- BIOMOD_LoadModels(myBiomodModelOut, models = sdm)
                                
            				# Get response plot from biomod2's special function
                    nulldev() # allows not to display the automatic plots from response.plot2 (try without once though)
            				myRespPlot2D <- response.plot2(models = response.plot,
            								Data = get_formal_data(myBiomodModelOut,'expl.var'), plot = F,
            								show.variables = get_formal_data(myBiomodModelOut,'expl.var.names'),
            								do.bivariate = FALSE,fixed.var.metric= 'mean',legend = TRUE,
            								data_species = get_formal_data(myBiomodModelOut,'resp.var'),
            								save.file = "jpeg",
            								name = paste("Response Curves for ", myRespName, ", ",sdm),
            								size = 600
                            ) # eo myRespPlot2D
                    #dev.off()	

            				# Transform into a data.frame for ggploting							   
            				resp <- data.frame(myRespPlot2D)
            				
            				#casting the data into sep col but keep env data, super ugly and long but it works
            				resp$id <- as.integer(resp$id)
            				
            				
            				#logChl
            				subset_logChl <- subset(resp, expl.name == "logChl")
            				cast_logChl <- reshape2::dcast(subset_logChl, id + expl.val ~  expl.name + pred.name , value.var = c("pred.val"))
            				names(cast_logChl)[names(cast_logChl) == "expl.val"] <- "logChl"
            				#reorder columns
            				cast_logChl <- cast_logChl[, c(1, 2, 3,4, 5, 6, 7, 8 , 9, 10)]
            				
            				#logEKE
            				subset_logEKE <- subset(resp, expl.name == "logEKE")
            				cast_logEKE <- reshape2::dcast(subset_logEKE, id + expl.val ~  expl.name + pred.name , value.var = c("pred.val"))
            				names(cast_logEKE)[names(cast_logEKE) == "expl.val"] <- "logEKE"
            				#reorder columns
            				cast_logEKE <- cast_logEKE[, c(1, 2, 3,4, 5, 6, 7, 8 , 9, 10)]
            				
            				#MLD
            				subset_MLD <- subset(resp, expl.name == "MLD")
            				cast_MLD <- reshape2::dcast(subset_MLD, id + expl.val ~  expl.name + pred.name , value.var = c("pred.val"))
            				names(cast_MLD)[names(cast_MLD) == "expl.val"] <- "MLD"
            				#reorder columns
            				cast_MLD <- cast_MLD[, c(1, 2, 3,4, 5, 6, 7, 8 , 9, 10)]
            				
            				#Nstar
            				subset_Nstar <- subset(resp, expl.name == "Nstar")
            				cast_Nstar <- reshape2::dcast(subset_Nstar, id + expl.val ~  expl.name + pred.name , value.var = c("pred.val"))
            				names(cast_Nstar)[names(cast_Nstar) == "expl.val"] <- "Nstar"
            				#reorder columns
            				cast_Nstar <- cast_Nstar[, c(1, 2, 3, 4, 5, 6, 7, 8 , 9, 10)]
            				
            				#PAR
            				subset_PAR <- subset(resp, expl.name == "PAR")
            				cast_PAR <- reshape2::dcast(subset_PAR, id + expl.val ~  expl.name + pred.name , value.var = c("pred.val"))
            				names(cast_PAR)[names(cast_PAR) == "expl.val"] <- "PAR"
            				#reorder columns
            				cast_PAR <- cast_PAR[, c(1, 2, 3, 4, 5, 6, 7, 8 , 9, 10)]
            				
            				#Sistar
            				subset_Sistar <- subset(resp, expl.name == "Sistar")
            				cast_Sistar <- reshape2::dcast(subset_Sistar, id + expl.val ~  expl.name + pred.name , value.var = c("pred.val"))
            				names(cast_Sistar)[names(cast_Sistar) == "expl.val"] <- "Sistar"
            				#reorder columns
            				cast_Sistar <- cast_Sistar[, c(1, 2, 3, 4, 5, 6, 7, 8 , 9, 10)]
            				
            				#SST
            				subset_SST <- subset(resp, expl.name == "SST")
            				cast_SST <- reshape2::dcast(subset_SST, id + expl.val ~  expl.name + pred.name , value.var = c("pred.val"))
            				names(cast_SST)[names(cast_SST) == "expl.val"] <- "SST"
            				#reorder columns
            				cast_SST <- cast_SST[, c(1, 2, 3, 4, 5, 6, 7, 8 , 9, 10)]
            				
            				#logNO3
            				subset_logNO3 <- subset(resp, expl.name == "logNO3")
            				cast_logNO3 <- reshape2::dcast(subset_logNO3, id + expl.val ~  expl.name + pred.name , value.var = c("pred.val"))
            				names(cast_logNO3)[names(cast_logNO3) == "expl.val"] <- "logNO3"
            				#reorder columns
            				cast_logNO3 <- cast_logNO3[, c(1, 2, 3,4, 5, 6, 7, 8 , 9, 10)]
            				
            				casted_data <- cbind(cast_logChl, cast_logEKE[,2:10], cast_logNO3[,2:10], cast_MLD[,2:10], cast_Nstar[,2:10], cast_PAR[,2:10], cast_Sistar[,2:10], cast_SST[,2:10])
            				casted_data2 <- cbind(cast_SST[,2], cast_logChl[2], cast_logNO3[,2], cast_MLD[,2], cast_Sistar[,2], cast_PAR[,2], cast_Nstar[,2], cast_logEKE[,2], cast_SST[,3:10], cast_logChl[,3:10], cast_logNO3[,3:10], cast_MLD[,3:10], cast_Sistar[,3:10], cast_PAR[,3:10], cast_Nstar[,3:10], cast_logEKE[,3:10])
            				
            				
            				names(casted_data2)[names(casted_data2) == "cast_logEKE[, 2]"] <- "logEKE"
            				names(casted_data2)[names(casted_data2) == "cast_logNO3[, 2]"] <- "logNO3"
            				names(casted_data2)[names(casted_data2) == "cast_MLD[, 2]"] <- "MLD"
            				names(casted_data2)[names(casted_data2) == "cast_Nstar[, 2]"] <- "Nstar"
            				names(casted_data2)[names(casted_data2) == "cast_PAR[, 2]"] <- "PAR"
            				
            				names(casted_data2)[names(casted_data2) == "cast_Sistar[, 2]"] <- "Sistar"
            				names(casted_data2)[names(casted_data2) == "cast_SST[, 2]"] <- "SST"
            				names(casted_data2)[names(casted_data2) == "logChl[, 2]"] <- "logChl"
            				
            				colnames(casted_data2[,9:88]) <- paste(rep(vars, each = 10, times = 1)) 

            				resp  <- casted_data2
            				
                            avg.resp <- lapply(vars, function(v) {
                                        #v <- "SST" # for testing
                                        sub <- resp[,grep(v, colnames(resp))]
                                        sub[,paste(v,"_avg_resp",sep="")] <- rowMeans(sub[,c(2:11)], na.rm = F) # average across the 10 CV runs
                                        return(sub[,c(1,12)])	    
                                }
                            ) # eo lapply - v in vars     
                            # cbind
                            table.resp <- do.call(cbind,avg.resp)     
                            rm(avg.resp,myRespPlot2D,resp) ; gc()

                            ### Derive univariate niche traits ! (Brun et al., 2015 and Benedetti et al., 2018)
                            require("parallel")
                            niche.traits <- mclapply(vars, function(v) {
                                        
                                    #v <- "Sistar" # for testing
                                    sub <- table.resp[,grep(v, colnames(table.resp))]
                                    colnames(sub) <- c('x','y')
        								
                                    ### Niche center = median
                    								sub$area <- cumsum(sub$y)
                    								medarea <- (sub[length(sub$area),"area"] / 2)
                    								sub$dist1 <- abs(sub$area - medarea)
                    								center <- sub[which(sub$dist1 == min(sub$dist1)),'x']
                    								if( length(center) > 1) {
                    										center <- mean(center)
                    								} # eo if loop
                                        
                                                    ### Find x value where y is maximal 
                                                    center2 <- sub[which(sub$y == max(sub$y)),'x']
                    								if( length(center2) > 1) {
                    										center2 <- mean(center2)
                    								} # eo if loop
                                        
                    								### Niche breadth = interquantile range between the 10th and the 90th
                    								q1 <- (sub[length(sub$area),"area"] / 10)
                    								q3 <- (9*(sub[length(sub$area),"area"]) / 10)
                    								sub$dist2 <- abs(sub$area - q1)
                    								sub$dist3 <- abs(sub$area - q3)
                    								lower <- sub[which(sub$dist2 == min(sub$dist2)),'x']
                    								if( length(lower) > 1) {
                    										lower <- mean(lower)
                    								} # eo if loop
                    								upper <- sub[which(sub$dist3 == min(sub$dist3)),'x']
                    								if( length(upper) > 1) {
                    										upper <- mean(upper)
                    								} # eo if loop
                    								breadth <- upper - lower
                                        
                                                    # And try also weighted mean and sd to estimate center and breadth
                                                    require("Hmisc")
                                                    wtd.quantiles <- wtd.quantile(x = sub$x, weights = ((sub$y)/max(sub$y)))
                                                    center3 <- wtd.quantiles[3]  # weighted median
                                                    lower2 <- wtd.quantile(x = sub$x, weights = ((sub$y)/max(sub$y)), probs = seq(from = 0, to = 1, by = 0.1) )[2]
                                                    upper2 <- wtd.quantile(x = sub$x, weights = ((sub$y)/max(sub$y)), probs = seq(from = 0, to = 1, by = 0.1) )[9]
                                                    breadth2 <- upper2 - lower2
                                        
                                                    # And compute the range of 'y' to idetify which covariates have flat responses or not
                                                    # (could be used to weight in the PCA, or whatever analysis of niche traits variations)
                                                    range <- max(sub$y) - min(sub$y)
                                        
                                                    ### Summarize in a data.frame that you will return
                                                    # center = niche center based on median derived from cumsums, like in Benedetti et al., 2018
                                                    # center2 = value where 'y' (mean HSI) is the highest (rename 'optim') --> USELESS
                                                    # center3 = niche center based on weighted median (weights = normalized HSI) (rename center2)
                                                    # lower = 10th quantile based on cumsum area
                                                    # upper = 90th quantile based on cumsum area
                                                    # breadth = upper - lower
                                                    # lower2 = weighted quantiles
                                                    # upper2 = weighted quantiles
                                                    # breadth2 = upper2 - lower2
                                                                                
                                                    traits <- data.frame(species = unique(data[data$obs == 1,"species"]), group = g,
                                                        var = v, center = center, width = breadth, lower = lower, upper = upper, 
                                                        optim = center2, center2 = center3, width2 = breadth2, lower2 = lower2, 
                                                        upper2 = upper2, prob.range = range)
                                            
                                                    # Return niche traits ; not plot yet 
                                                    return(traits)
        								
                                                }, mc.cores =  1
                                    
                                ) # eo lapply
                                
                                ### Bind species straits
                                table.niche.traits <- bind_rows(niche.traits)
                                # Idea: if you plan to choose the niche traits in a PCA-like analysis, you might want to weight in their contribution according to how much HSI varies (so prob.range)
                                # In env niche space the variables that present stronger resp curves should have higher weight; other a niceh center derived from a rather flat response curve
                                # will have the same weight a niche center that actually has a strong responses and thus a more meaningful niche center
                                table.niche.traits$weight.pca <- (table.niche.traits$prob.range)/ max(table.niche.traits$prob.range)    
                                         
                                rm(niche.traits) ; gc()	
                                
                                ### Make a panel plot with all of species resp curves (n = 8 per species, change accordingly to your choice) and highlighting their niche traits
                                ### Make use of 'table.resp' and 'table.niche.traits' in a lapply and return a 'list' that will contain the 8 plots, 
                                ### arrange them in a panel with ggarrange
                                niche.plots <- lapply(vars, function(v) { 
                                                #v = "logNO3"
                                                sub <- table.resp[,grep(v, colnames(table.resp))]
                                                colnames(sub) <- c('x','y')
                                                sub.traits <- table.niche.traits[table.niche.traits$var == v,]
                                            
                                                plot <- ggplot() + geom_path(aes(x = x, y = y), linetype = "solid", data = sub, group = 1) +
                                                            geom_vline(xintercept = sub.traits$optim, colour = "#a50026", linetype = "solid") +
                                                            geom_vline(xintercept = sub.traits$center, colour = "#d73027", linetype = "solid") +
                                                            geom_vline(xintercept = sub.traits$center2, colour = "#313695", linetype = "solid") +
                                                            geom_vline(xintercept = sub.traits$lower, colour = "#f46d43", linetype = "dashed") +
                                                            geom_vline(xintercept = sub.traits$upper, colour = "#f46d43", linetype = "dashed") +
                                                            geom_vline(xintercept = sub.traits$lower2, colour = "#74add1", linetype = "dashed") +
                                                            geom_vline(xintercept = sub.traits$upper2, colour = "#74add1", linetype = "dashed") +
                                                            xlab(v) + ylab("Mean HSI") + scale_y_continuous(limits = c(0,1)) +
                                                            theme_classic()
                                                      
                                                # Return 
                                                return(plot)
                                
                                        } # eo FUN 
                                    
                                ) # eo lapply - v in vars
                                
                                require("ggpubr")
                                panel <- ggarrange(niche.plots[[1]],niche.plots[[2]],niche.plots[[3]],niche.plots[[4]],
                                            niche.plots[[5]],niche.plots[[6]],niche.plots[[7]],niche.plots[[8]],
                                            ncol = 4, nrow = 2, align = "hv")
                                         
                                message(paste("", sep = ""))
                                message(paste("Saving plots and tables for ",myRespName, sep = ""))
                                message(paste("", sep = ""))
                                         
                                ### Print panel in appropriate dir
                                setwd(paste(data.wd.back,"/resp_curves/plots/", sep = "")) # change accordingly
                                ggsave(plot = panel, filename = paste("panel_respcurv_",sdm,"_group_",unique(data[data$obs == 1,"species"]),".jpg", sep = ""), dpi = 300, height = 4, width = 10)
                                
                                ### Save niche traits table
                                setwd(paste(data.wd.back,"/resp_curves/niche_traits/", sep = "")) # change accordingly
            					          save(table.niche.traits, file = paste("table_niche_traits_",sdm,"_group_",unique(data[data$obs == 1,"species"]),".Rdata", sep = ""))
 
                                rm(niche.plots, table.niche.traits, table.resp)
                                gc()
                                setwd(second.wd)		
                            
                        } # eo for loop sdm in SDMs
                        
    					### Save evaluation scores again, should be nearly the same as before  
    					setwd(paste(data.wd.back,"/","Scores", sep = ""))
    					save(scores.tbl, file = paste("eval_scores_for_respcurv_",sdm,"_",unique(data[data$obs == 1,"species"]),".Rdata", sep = ""))
                        rm(scores.tbl)
                        gc()
                        setwd(second.wd)
								
        			} else {
							
        				    message(paste(" NOT ENOUGH PRESENCES WITH ENV PREDICTORS VALUES || n = ", n, sep = ""))
								
        			}
        
    } # eo for loop - f in files
    
    
# } # eo for loop - cat in categories
    
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
