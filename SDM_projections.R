#Rscript by Fabio Benedetti, adapted by Jonas Wydler for Master thesis
# - SDM projections, creates community tables for each SDM
#note: species "Scaphocalanus.echinatus" fails for GLM
#30.11.2020
#updated:04.09.2021
### ================================================================================================================

#library("raster")
#library("sp")
#library("stringr")
#library("reshape2")
library("tidyverse")
library("biomod2") 
#library("R.devices")



wd.back <- ("wd_data/species_background_data/")
wd.clim <- ("wd_data/global_monthly_clims_1d/")
wd.output.GLM <- ("wd_data/species_background_data/niche.modelling_GLM/")
wd.output.GAM <- ("wd_data/species_background_data/niche.modelling_GAM/")
wd.output.ANN <- ("wd_data/species_background_data/niche.modelling_ANN/")

Sys.setenv(LANGUAGE= 'en')

### ================================================================================================================

### 1°) Run SDM projections by adjusting the predictors accordng to SDMs and functional group.

# Vector of SDMs
SDMs <- c('GAM') # you can play around with the SDms of course, just make sure to change the script below accordingly
# Vector of eval_runs :
eval_runs <- c("RUN1","RUN2","RUN3","RUN4","RUN5","RUN6","RUN7","RUN8","RUN9","RUN10") 
# vector of months, for 2°)
months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")

### Load the env stacks (winter  and summer, like January & August) for projection
setwd(wd.clim) # the dir where the 12 monthly climatologies used to fit the SDMs are stored
apr <- read.table("glob_stack_month_apr_18_09_20.txt", h = T, sep = ";")
jul <- read.table("glob_stack_month_jul_18_09_20.txt", h = T, sep = ";")
oct <- read.table("glob_stack_month_oct_18_09_20.txt", h = T, sep = ";")
jan <- read.table("glob_stack_month_jan_18_09_20.txt", h = T, sep = ";")

apr <- apr[-which(apr$SSS < 20),] 
apr <- apr[-which(apr$Bathy > -175),] 
jul <- jul[-which(jul$SSS < 20),] 
jul <- jul[-which(jul$Bathy > -175),] 
oct <- oct[-which(oct$SSS < 20),] 
oct <- oct[-which(oct$Bathy > -175),] 
jan <- jan[-which(jan$SSS < 20),] 
jan <- jan[-which(jan$Bathy > -175),] 

### Load the remaining 8 climatologies for annual proj of diversity
annual <- TRUE
# annual <- FALSE

if(annual == TRUE) {
		
		feb <- read.table("glob_stack_month_feb_18_09_20.txt", h = T, sep = ";")
		mar <- read.table("glob_stack_month_mar_18_09_20.txt", h = T, sep = ";")
		may <- read.table("glob_stack_month_may_18_09_20.txt", h = T, sep = ";")
		jun <- read.table("glob_stack_month_jun_18_09_20.txt", h = T, sep = ";")
		aug <- read.table("glob_stack_month_aug_18_09_20.txt", h = T, sep = ";")
		sep <- read.table("glob_stack_month_sep_18_09_20.txt", h = T, sep = ";")
		nov <- read.table("glob_stack_month_nov_18_09_20.txt", h = T, sep = ";")
		dec <- read.table("glob_stack_month_dec_18_09_20.txt", h = T, sep = ";")
	
		feb <- feb[-which(feb$SSS < 20),] 
		feb <- feb[-which(feb$Bathy > -175),] 
		mar <- mar[-which(mar$SSS < 20),] 
		mar <- mar[-which(mar$Bathy > -175),] 
		may <- may[-which(may$SSS < 20),] 
		may <- may[-which(may$Bathy > -175),] 
		jun <- jun[-which(jun$SSS < 20),] 
		jun <- jun[-which(jun$Bathy > -175),] 
		aug <- aug[-which(aug$SSS < 20),] 
		aug <- aug[-which(aug$Bathy > -175),] 
		sep <- sep[-which(sep$SSS < 20),] 
		sep <- sep[-which(sep$Bathy > -175),] 
		nov <- nov[-which(nov$SSS < 20),] 
		nov <- nov[-which(nov$Bathy > -175),] 
		dec <- dec[-which(dec$SSS < 20),] 
		dec <- dec[-which(dec$Bathy > -175),]

} # eo if loop : annual switch


setwd(wd.back)


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
								, outerPIsteps = 0, idLinksBases = TRUE, scalePenalty = TRUE, keepData = FALSE) )
                        
) # eo modelling options



SpeciesNicheModelling <- function(f = files) {
	
	                        # Get the data
	                        setwd(wd.back)
	                        message(paste("Modelling ",f, " =========================================================", sep = ""))
	                        
	                        data <- read.table(f, h = T, sep = ";")
	                        
                            g <- unique(data[data$obs == 1,"FG_f11"]) # replace 'group2' by your FG column
                            
                            ### Change data colnames for: logChla --> logChl  ### CHECK IF THIS USEFUL IN YOUR CASE
                            colnames(data)[grep("logChla", colnames(data))] <- "logChl"
                            
                            colnames(data)[grep("Chla", colnames(data))] <- "Chl"
                            
                            ### Define vector of predictors as a fucntion of the species' FG and the SDM according to ranking
                            
                            #test if the value for "group" is empty or NA
                            {if ((length(g) == 0) || (is.na(g) == T)){
                              message(paste("NO FG assigned, general vars assinged"))
                              #in this case: assign best general set of predictors
                              vars <- c("SST","PAR", "MLD", "logNO3", "logChl")
                            }else if(g == "1" & sdm == "ANN") {
                                vars <- c("SST","PAR","MLD","Nstar","Sistar")
                            } else if(g == "1" & sdm == "GAM") {
                                vars <- c("SST","logNO3","logChl","PAR","logEKE") 
                            } else if(g == "1" & sdm == "GLM") {
                                vars <- c( "SST","logNO3","logChl","Nstar","logEKE")
                            } else if(g == "2" & sdm == "ANN" ) {
                                vars <- c("SST","MLD","PAR","Nstar","Sistar")
                            } else if(g == "2" & sdm == "GAM" ) {
                                vars <- c("SST","logNO3","logChl","PAR","logEKE" )
                            } else if(g == "2" & sdm == "GLM" ) {
                                vars <- c("SST", "logNO3", "logChl", "logEKE", "Nstar") 
                            } else if(g == "3" & sdm == "ANN" ) {
                              vars <- c("SST", "PAR", "MLD", "logEKE", "Nstar")
                            } else if(g == "3" & sdm == "GAM" ) {
                              vars <- c("SST", "logEKE", "logNO3", "logChl", "Sistar")
                            } else if(g == "3" & sdm == "GLM" ) {
                              vars <- c("logEKE", "SST", "logNO3", "logChl", "Nstar") 
                            } else if(g == "4" & sdm == "ANN" ) {
                              vars <- c("SST", "PAR", "MLD", "Nstar", "Sistar")
                            } else if(g == "4" & sdm == "GAM" ) {
                              vars <- c("SST", "logChl", "PAR", "logEKE", "logNO3")
                            } else if(g == "4" & sdm == "GLM" ) {
                              vars <- c("logChl", "logEKE", "SST", "logNO3", "PAR") 
                            } else if(g == "5" & sdm == "ANN" ) {
                              vars <- c("SST", "PAR", "MLD", "Sistar", "logEKE")
                            } else if(g == "5" & sdm == "GAM" ) {
                              vars <- c("SST", "logEKE", "logChl", "Sistar", "Wind")
                            } else if(g == "5" & sdm == "GLM" ) {
                              vars <- c("SST", "logEKE", "logChl", "Wind", "logNO3") 
                            } else if(g == "6" & sdm == "ANN" ) {
                              vars <- c("SST", "MLD", "PAR", "Sistar", "Wind")
                            } else if(g == "6" & sdm == "GAM" ) {
                              vars <- c("logChl", "logEKE", "SST", "pCO2", "PAR")
                            } else if(g == "6" & sdm == "GLM" ) {
                              vars <- c("logChl", "logEKE", "logNO3", "pCO2", "SST")#SST was forced here instead of sistar 
                            } else if(g == "7" & sdm == "ANN" ) {
                              vars <- c("SST", "MLD", "PAR", "Sistar", "Nstar")
                            } else if(g == "7" & sdm == "GAM" ) {
                              vars <- c("SST", "logNO3", "logChl", "PAR", "logEKE")
                            } else if(g == "7" & sdm == "GLM" ) {
                              vars <- c("SST", "logNO3", "logChl", "logEKE", "Nstar") 
                            } else if(g == "8" & sdm == "ANN" ) {
                              vars <- c("SST", "MLD", "PAR", "Nstar", "Sistar")
                            } else if(g == "8" & sdm == "GAM" ) {
                              vars <- c("SST", "logNO3", "logEKE", "logChl", "MLD")
                            } else if(g == "8" & sdm == "GLM" ) {
                              vars <- c("logNO3", "SST", "logChl", "PAR", "logEKE") 
                            } else if(g == "9" & sdm == "ANN" ) {
                              vars <- c("SST", "PAR", "MLD", "Sistar", "logEKE")
                            } else if(g == "9" & sdm == "GAM" ) {
                              vars <- c("SST", "logNO3", "logChl", "logEKE", "MLD")
                            } else if(g == "9" & sdm == "GLM" ) {
                              vars <- c("SST", "logNO3", "logChl", "logEKE", "Nstar") 
                            } else if(g == "10" & sdm == "ANN" ) {
                              vars <- c("SST", "logNO3", "Sistar", "MLD", "logChl")
                            } else if(g == "10" & sdm == "GAM" ) {
                              vars <- c("SST", "logNO3", "Sistar", "logChl", "Nstar")
                            } else if(g == "10" & sdm == "GLM" ) {
                              vars <- c("SST", "logNO3", "logChl", "Sistar", "Nstar") 
                            } else if(g == "11" & sdm == "ANN" ) {
                              vars <- c("SST", "PAR", "MLD", "Sistar", "logEKE")
                            } else if(g == "11" & sdm == "GAM" ) {
                              vars <- c("SST", "logChl", "logEKE", "logNO3", "Sistar")
                            } else if(g == "11" & sdm == "GLM" ) {
                              vars <- c("SST", "logChl", "logNO3", "logEKE", "Nstar") 
                               }else {
                            
                            } 
                            }
                            
                            
                            
	                        data2 <- na.omit(data[,c("species","x","y","obs",vars)] )
	                        n <- nrow( data2[data2$obs == 1,] )
	                        rm(data2)
                            
							# If n >= 50, continue, otherwise got to next species
							if( n >= 50 ) {
								
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
		
                                
	  			  				myBiomodData <- BIOMOD_FormatingData(resp.var = myResp, expl.var = myExpl, resp.xy = myRespXY, resp.name = myRespName)
							
	 			   				### Modelling
                                nulldev()
								myBiomodModelOut <- BIOMOD_Modeling(myBiomodData, models = sdm, models.options = myBiomodOption, NbRunEval = length(eval_runs), 
								   DataSplit = 80, Yweights = weights, VarImport = 0, models.eval.meth = c('ROC','TSS'), SaveObj = F, do.full.models = F
								) # eo modelling
                                dev.off()					
  
								### Get evaluation scores
								scores <- data.frame(get_evaluations(myBiomodModelOut))
								# Use lapply to summarize all scores (10 scores per 6 SDMs) AND PROBABILITY THRESH
								scores.ls <- lapply(sdm, function(sdms) {
												
                                                # Retrieve all scores for the sdm 's' 
												tss <- scores["TSS",grep(paste("Testing.data.",sdms, sep = ""), colnames(scores))]
												tss_cutoffs <- scores["TSS",grep(paste("Cutoff.",sdms, sep = ""), colnames(scores))]
												colnames(tss) <- paste(sdms, eval_runs, sep = "_")
												colnames(tss) <- paste(sdms, eval_runs, sep = "_")
												tss <- t(tss)
												tss_cutoffs <- t(tss_cutoffs)
												colnames(tss_cutoffs) <- "Cutoff_TSS"
                        
                                                
                        #retrieve auc scores
												auc <- scores["ROC",grep(paste("Testing.data.",sdm, sep = ""), colnames(scores))]
												auc_cutoffs <- scores["ROC",grep(paste("Cutoff.",sdm, sep = ""), colnames(scores))]
												colnames(auc) <- paste(sdm, eval_runs, sep = "_")
												colnames(auc_cutoffs) <- paste(sdm, eval_runs, sep = "_")
												auc <- t(auc)
												auc_cutoffs <- t(auc_cutoffs)
												colnames(auc_cutoffs) <- "Cutoff_AUC"
												                        
												return(cbind(tss, tss_cutoffs, auc, auc_cutoffs))
												
								} ) # eo lapply
								scores.tbl <- data.frame(do.call(rbind, scores.ls))
                                scores.tbl$species <- unique(data[data$obs == 1,"species"])
                                scores.tbl$group <- g
								rm(scores,scores.ls)
							
	  			  				### Make projections (they will be printed in the species folder)
								message(paste("Projecting ",myRespName, " =========================================================", sep = ""))
                                
                                # dim(apr[,vars]) ; summary(apr[,vars])                       
                                 myBiomodProj <- BIOMOD_Projection(
                                                     modeling.output = myBiomodModelOut,
                                                     new.env = apr[,vars],
                                                     proj.name = paste("projection",myRespName,"apr",sdm, sep = "_"),
                                                     selected.models = 'all',
                                                     binary.meth = 'TSS',
                                                     compress = 'xz',
                                                     clamping.mask = FALSE
                                 ) # eo projection
    
                                 # Project niches in July conditions
                                 myBiomodProj <- BIOMOD_Projection(
                                                     modeling.output = myBiomodModelOut,
                                                     new.env = jul[,vars],
                                                     proj.name = paste("projection",myRespName,"jul",sdm, sep = "_"),
                                                     selected.models = 'all',
                                                     binary.meth = 'TSS',
                                                     compress = 'xz',
                                                     clamping.mask = FALSE
                                 ) # eo projection
    
                                 # Project niches in October conditions
                                 myBiomodProj <- BIOMOD_Projection(
                                                     modeling.output = myBiomodModelOut,
                                                     new.env = oct[,vars],
                                                     proj.name = paste("projection",myRespName,"oct",sdm, sep = "_"),
                                                     selected.models = 'all',
                                                     binary.meth = 'TSS',
                                                     compress = 'xz',
                                                     clamping.mask = FALSE
                                 ) # eo projection
    
                                 # Project niches in January conditions
                                 myBiomodProj <- BIOMOD_Projection(
                                                     modeling.output = myBiomodModelOut,
                                                     new.env = jan[,vars],
                                                     proj.name = paste("projection",myRespName,"jan",sdm, sep = "_"),
                                                     selected.models = 'all',
                                                     binary.meth = 'TSS',
                                                     compress = 'xz',
                                                     clamping.mask = FALSE
                                 ) # eo projection
							
								### Create a switch for annual projections (annual = T/F)
								if( annual == TRUE ) {
                                    
									# Feb
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = feb[,vars],
														proj.name = paste("projection",myRespName,"feb",sdm, sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									) # eo projection
									# Mar
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = mar[,vars],
														proj.name = paste("projection",myRespName,"mar",sdm, sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									) # eo projection
									# May
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = may[,vars],
														proj.name = paste("projection",myRespName,"may",sdm, sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									) # eo projection
									# Jun
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = jun[,vars],
														proj.name = paste("projection",myRespName,"jun",sdm, sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									) # eo projection
									# Aug
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = aug[,vars],
														proj.name = paste("projection",myRespName,"aug",sdm, sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									) # eo projection	
									# Sep
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = sep[,vars],
														proj.name = paste("projection",myRespName,"sep",sdm, sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									)  # eo projection
									# Nov
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = nov[,vars],
														proj.name = paste("projection",myRespName,"nov",sdm, sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									)  # eo projection
									# Dec
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = dec[,vars],
														proj.name = paste("projection",myRespName,"dec",sdm, sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									)  # eo projection
									
								} else { 
                                    
									message(paste("============== NOT PERFORMING ANNUAL PROJECTIONS ==============", sep = ""))
                                    
								} # eo if else loop

								### Save evaluation scores
								setwd(paste(wd.back,"/","Scores", sep = ""))
								save(scores.tbl, file = paste("eval_scores_for_projections_test_",unique(data[data$obs == 1,"species"]),"_",sdm,".Rdata", sep = ""))
								setwd(wd.back)			
								
							} else {
							
								message(paste("WRONG GROUP ('Other_') OR JUST NOT ENOUGH PRESENCES WITH ENV PREDICTORS VALUES || n = ", n, sep = ""))
								
							}
  
} # eo NicheModelling


### For loop to apply the function above in mclapply per SDMs

#shorter version without check for failed runs. 
setwd(wd.back)

files <- dir()[grep("backgrounddata_",dir())]
#"Scaphocalanus.echinatus" fails 
# sdm <- SDMs[1]
for (sdm in SDMs){
  setwd(wd.back)
  message(paste(" ", sep = ""))	
  message(paste("RUNNING ",sdm," ", sep = ""))
  message(paste(" ", sep = ""))
  
  require("parallel")
  mclapply(X = files, SpeciesNicheModelling, mc.cores = 1)
  
  setwd(wd.back)
}


### ==============================================================================================================================


### 2°) Mclapply that allows yout ot retrieve the mean SDM species projection and retunr a community matrix (one for each month and SDM)

baseline.community.extracter <- function(m = months) {

              message(paste(" ", sep = ""))
              message(paste("Retrieving baseline copepod FG probabilities for ", m, sep = ""))
              message(paste(" ", sep = ""))

              # Load env variables (for coordinates, they are not saved by biomod in the projections .Rdata files, unfortunately)
              setwd(wd.clim)
              env <- read.table(paste("glob_stack_month_",m,"_18_09_20.txt", sep = ""), h = T, sep = ";")
              env <- env[-which(env$SSS < 20),]
              env <- env[-which(env$Bathy > -175),]

              # And for each SDM
              for(sdm in SDMs) {

                      message(paste("SDM == ",sdm, sep = ""))

                      ### Go to the niche.modelling directory of SDM == sdm
                      setwd(paste(wd.back,"niche.modelling_",sdm, sep = ""))  # change accordingly

                      ### Retrieve the vector of species folder (one per successfully projected species)
                      spp <- dir()
                      spp <- gsub("\\(|\\)", "", spp)

                      probas <- lapply(X = spp, FUN = function(sp) {

                                  # Got to species projection dir
                                  setwd(paste(wd.back,"niche.modelling_",sdm,"/",sp,"/", sep = ""))
                                  message(paste("Loading projections for ", sp, "  ================================  ", sep = ""))

                                  # If the 12 monthly projections are done
                                  if(sum(grepl("proj_projection_", dir())) == 12 ) {
                                        setwd(paste(wd.back,"/niche.modelling_",sdm,"/",sp,"/", sep = "")) 
                                        setwd(paste(paste("proj_projection_",sp,"_",m,"_",sdm, sep = ""),"/", sep = "")) # check if this is correct
                                        d <- get(load(paste("proj_projection_",sp,"_",m,"_",sdm,"_",sp, ".RData", sep = "")))  # this too
                                        resModel <- d[,sdm,,]
                                        resModel <- apply(resModel, 1, mean, na.rm = F) # retrieve the average of the 10 CV runs, dont remove NAs (cells with missing vars or land)
                                        resModel <- (resModel/1000) # re-scale to 0-1

                                         # Return
                                         return( data.frame(cell_id = paste(env$x, env$y, sep = "_"), x = env$x, y = env$y, species = gsub("\\.","_",sp), HSI = resModel) )

                                      } else {

                                          message(paste("Skipping because not projection (yet) for ", sp, "  ================================", sep = ""))

                                      } # eo if else loop

                              } # eo FUN

                      ) # eo lapply
                      tbl <- dplyr::bind_rows(probas)
                      rm(probas); gc()

                      ### Dcast to put species' average habitat suitability index (HSI) as columns
                      d_tbl <- dcast(tbl, cell_id + x + y ~ species, fun.aggregate = mean, na.rm = T, value.var = "HSI")
                      d_tbl$SDM <- sdm

                      # Save !
                      message(paste("Saving community table for | ",sdm," | ",m, sep = ""))
                      message(paste(" ", sep = ""))

                      ### Go to the dir where you want to save community tables
                      setwd(paste0(wd.back,"ct_",sdm,"/"))
                      write.table(d_tbl, file = paste("table_mon_composition_baseline_",sdm,"_",m,".txt", sep = ""), sep = "\t")

                      setwd(wd.back)
                      rm(d_tbl,tbl) ; gc()

              } # eo for loop - sdm in SDMs

} # eo FUN - baseline.community.extracter

### mclapply
require("parallel")

mclapply(X = months[1], FUN = baseline.community.extracter, mc.cores = 1)


### ==============================================================================================================================
### ==============================================================================================================================
### ==============================================================================================================================
