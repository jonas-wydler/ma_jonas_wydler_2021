# ETH - Script by Fabio Benedetti, adapated by Jonas Wydler 
#  - Rank environmental predictor
#Updated: For legacy on 01.09.2021
# --------------------------------------------------------------------------------------------------------------------------------
library("biomod2")
#library("rJava")
library("tidyverse")
library("stringr")
library("reshape2")
library("geosphere")
library("ncdf4")
library("classInt")
library("parallel")
#library("corrgram")
library("biomod2")
#library("rJava") # necessary in case you use maxent. Warning: there are often conflicts between libraraies and the Java version you jave installed on your computer

 
Sys.setenv(LANGUAGE= 'en')
setwd(paste0("wd_data/species_background_data/")) # go to dir with background data
# --------------------------------------------------------------------------------------------------------------------------------

myBiomodOption <- BIOMOD_ModelingOptions(
	
						GLM = list(type = 'quadratic', interaction.level = 0, myFormula = NULL, test = 'AIC', family = binomial("logit"),
    								mustart = 0.5, control = glm.control(epsilon = 1e-08, maxit = 50, trace = FALSE)),

						GAM = list(algo = 'GAM_mgcv', type = 's_smoother',
								k = 5,interaction.level = 0, myFormula = NULL,
								family = binomial("logit"),method = 'GCV.Cp',
								optimizer = c('outer','newton'), select = FALSE,
								knots = NULL, paraPen = NULL,
								control = list(nthreads = 1, irls.reg = 0, epsilon = 1e-07
								, maxit = 200, trace = FALSE, mgcv.tol = 1e-07, mgcv.half = 15
								, rank.tol = 1.49011611938477e-08
								, nlm = list(ndigit=7, gradtol=1e-06, stepmax=2, steptol=1e-04, iterlim=200, check.analyticals=0)
								, optim = list(factr=1e+07)
								, newton = list(conv.tol=1e-06, maxNstep=5, maxSstep=2, maxHalf=30, use.svd=0)
 								, outerPIsteps = 0, idLinksBases = TRUE, scalePenalty = TRUE, keepData = FALSE) )#, 

            				# MAXENT = list( path_to_maxent.jar = ("/net/kryo/work/wydlerj/background_data/"), # should not have ANY special char or space
            				# 		maximumiterations = 200,
            				# 		visible = FALSE,
            				# 		linear = TRUE,
            				# 		quadratic = TRUE,
            				# 		product = TRUE,
            				# 		threshold = FALSE,  # Disable threshold features. You wat that to avid fitting step-like resp curves
            				# 		hinge = TRUE,
            				# 		lq2lqptthreshold = 80,
            				# 		l2lqthreshold = 10,
            				# 		hingethreshold = 15,
            				# 		beta_threshold = -1,
            				# 		beta_categorical = -1,
            				# 		beta_lqp = -1,
            				# 		beta_hinge = -1,
            				# 		defaultprevalence = 0.5)
								
) # eo modelling options


VarsRanker <- function(f = files) {
	
							# Get the data
              setwd(paste0("wd_data/species_background_data/"))
              data <- read.table(f, h = T, sep = ";")
							data2 <- na.omit( data[,c(vars,"obs","species")] )
							n <- nrow( data2[data2$obs == 1,] )
              sp <- unique(data2[data2$obs == 1,"species"]) # species name as it's recorded in the file (not myRespName)
							rm(data2)
                            
              message(paste("Modelling ",sp, " =========================================================", sep = ""))
							
							# If n >= 50, continue, otherwise got to next species
							if(n >= 50) {
								
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

								# Environmental variables
								myExpl <- data[,vars]

								# weights vector
								weights <- na.omit(data[,c(vars,"weights")])[,"weights"]
		
								# formatage des données
	  			  		myBiomodData <- BIOMOD_FormatingData(resp.var = myResp, expl.var = myExpl, resp.xy = myRespXY, resp.name = myRespName)
                                
                                # Go to niche modelling directory (a folder called 'niche.modelling.ranks.test' within the species data directory)
                                ### WARNING: biomod2 really does not function if you have any special character or space in that character string
	  			  		                setwd(paste0("wd_data/species_background_data/"))
	  			  				
								
                                myBiomodModelOut <- BIOMOD_Modeling(myBiomodData, models = SDMs, models.options = myBiomodOption, 
								                         NbRunEval = length(eval_runs), DataSplit = 80, 
								                         Yweights = weights, VarImport = 30, models.eval.meth = c('ROC','TSS'),
								                         SaveObj = FALSE, do.full.models = FALSE
								) # eo modelling
                                
								### Get evaluation scores
								scores <- data.frame(get_evaluations(myBiomodModelOut))
								# Use lapply to summarize all scores (10 scores per 6 SDMs) AND PROBABILITY THRESH
								scores.ls <- lapply(SDMs, function(sdm) {
                                    
												# Retrieve all TSS scores for the sdm 's' 
												tss <- scores["TSS",grep(paste("Testing.data.",sdm, sep = ""), colnames(scores))]
												tss_cutoffs <- scores["TSS",grep(paste("Cutoff.",sdm, sep = ""), colnames(scores))]
												colnames(tss) <- paste(sdm, eval_runs, sep = "_")
												colnames(tss) <- paste(sdm, eval_runs, sep = "_")
												tss <- t(tss)
												tss_cutoffs <- t(tss_cutoffs)
												colnames(tss_cutoffs) <- "Cutoff_TSS"
                                                
												# Same with AUC/ROC
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
                scores.tbl$group <- unique(data[data$obs == 1,"group2"])
								rm(scores,scores.ls)
							
	  			  				### ?variables_importance
                                ranks <- data.frame(get_variables_importance(myBiomodModelOut))
                                ranks$predictor <- rownames(ranks)
                                m.ranks <- melt(ranks, id.var = "predictor")
                                colnames(m.ranks) <- c("predictor","model","rank")
                                m.ranks$model <- gsub(".", "_", m.ranks$model, fixed = T)
                                m.ranks$SDM <- data.frame(do.call(rbind,strsplit(as.character(m.ranks$model), split = "_")))[,1]
                                m.ranks$RUN <- data.frame(do.call(rbind,strsplit(as.character(m.ranks$model), split = "_")))[,2]
                                m.ranks$species <- unique(data[data$obs == 1,"species"])
                                m.ranks$group <- unique(data[data$obs == 1,"FG_f11"]) # FG here
                                                                    
								### Save predictors ranks in dedicated dir
                setwd(paste0("wd_data/Ranks/"))
                                
                save(m.ranks, file = paste("table_ranks_", sp,".Rdata", sep = ""))
                                
                ### Save evaluation scores in dedicated dir
                setwd(paste0("wd_data/Scores/"))
								save(scores.tbl, file = paste("table_scores_", sp,".Rdata", sep = ""))
                                
								setwd(paste0("wd_data/species_background_data/")) 			
								
								message(paste("Species ", sp, " Done.", sep = ""))
							} else {
							
								message(paste("NOT ENOUGH PRESENCES WITH ENV PREDICTORS VALUES || n = ", n, sep = ""))
								
							}
  
} # eo NicheModelling

require("parallel")


### SDM types to run
SDMs <- c('GLM','GAM','ANN')
#SDMs <- c('GAM', 'GLM')
# Avoid using RF and GBM for sure (overfitting), you can try MARS and FDA (though I can tell you both provide similar response curves as GAMs and GLMs) and MAXENT if you want. But these 3 should suffice. 

### Nb of cross-evaluations runs
eval_runs <- c("RUN1","RUN2","RUN3","RUN4","RUN5") # No need to run 10 CV runs here 
#eval_runs <- c("RUN1","RUN2")
vars <- c("SST","MLD","PAR","pCO2","logEKE","Nstar","logChla","Wind","Sistar","logNO3") # vector of environmental predictors names, make sure they match colnames in species files


setwd(paste0("wd_data/species_background_data/")) 	
files <- dir()[grep("backgrounddata_",dir())] # vector of filenames
#f <- files[32]
mclapply(X = files, FUN = VarsRanker, mc.cores = 32)


# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
