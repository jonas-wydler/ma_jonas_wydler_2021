
#Investigating SST scores to merge with distribution maps, 13.12.2020


library("tidyverse")
library("stringr")
library("reshape2")
library("geosphere")
library("ncdf4")
library("classInt")
library("parallel")
library("biomod2")
library("rJava")
library("forcats")
library("dplyr")


#===================================================================
wd_scores <- ("C:/Users/Jonas/ma/sdm/background/Scores/Scores_all_from_projections/")
wd_omit <- ("C:/Users/Jonas/ma/sdm/background/low_TSS_species_projections/")
wd_script <- ("C:/Users/Jonas/ma/sdm/part3/")


setwd(wd_scores)

scores <- dir()[grep("", dir())] ; scores
res <- lapply(scores,function(c){
  setwd(wd_scores)
  message(paste("Reading Scores ",c, sep = ""))
  data <- load(c)#carful with name, called scores.tb1 for me
  
  return(scores.tbl)
})
table_scores <- bind_rows(res); rm(res); rm(scores.tbl); gc()
table_scores <- tibble::rownames_to_column(table_scores, "model")

#this part is to get SDM and RUN out of rowname, it's not ideal and these variables should be defined in the model directly
table_scores$temp <- substr(x = as.character(table_scores$model), start = 1, stop =8)
table_scores$SDM <- data.frame(do.call(rbind,strsplit(as.character(table_scores$temp), split = "_")))[,1]
table_scores$RUN <- data.frame(do.call(rbind,strsplit(as.character(table_scores$temp), split = "_")))[,2]
table_scores$temp <-NULL
#===================================================================
df_GAM <- subset(table_scores, SDM == "GAM")
df_GLM <- subset(table_scores, SDM == "GLM")
df_ANN <- subset(table_scores, SDM == "ANN")

df_GAM_v2 <- df_GAM %>%
  group_by(species) %>%
  summarise(TSS = mean(TSS), tss_cutoff = mean(Cutoff_TSS))

df_GLM_v2 <- df_GLM %>%
  group_by(species) %>%
  summarise(TSS = mean(TSS), tss_cutoff = mean(Cutoff_TSS))

df_ANN_v2 <- df_ANN %>%
  group_by(species) %>%
  summarise(TSS = mean(TSS), tss_cutoff = mean(Cutoff_TSS))

df_GAM_v2$SDM <- "GAM"
df_GLM_v2$SDM <- "GLM"
df_ANN_v2$SDM <- "ANN"

species_omit_TSS_v2 <- rbind(df_GAM_v2, df_GLM_v2, df_ANN_v2)

setwd(wd_omit)

write.table(species_omit_TSS_v2, file = "species_omit_TSS_v2.txt", sep = "\t",
            row.names = FALSE)

setwd(wd_script)
save(df_GAM_v2,file="df_GAM_v2.Rda")
save(df_GLM_v2,file="df_GLM_v2.Rda")
save(df_ANN_v2,file="df_ANN_v2.Rda")

