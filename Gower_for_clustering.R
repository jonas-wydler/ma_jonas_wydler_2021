#Clustering of functional traits using Gower Distance, 12.10.2020

setwd("C:/Users/Jonas/Documents/ETH/master/ma/r/codev2/")
Sys.setenv(LANG = "en")

# --------------------------------------------------------------------------------------------------------------------
library("tidyverse")
library("FD")

# --------------------------------------------------------------------------------------------------------------------
### Load the functional traits table and explore
fct <- read.csv("table_funct_traits_copepods_v2.csv", h = T, sep = ";", dec = ",")
dim(fct) ; head(fct);  colnames(fct)
str(fct)
fct <- fct[,c(3:20)]

### Count NA
names <- colnames(fct)[c(7,8,9,10,15)] ; names
fct$na_count <- apply(fct[,names], 1, function(x) sum(is.na(x)))
# Drop species with missing body length and more than two missing traits
fct <- fct[!is.na(fct$max_body_length),]
fct <- fct[fct$na_count < 2,]

# --------------------------------------------------------------------------------------------------------------------
### Gower Distance matrix

fct <- fct[]
#convert to 1 and 0
fct$Myelination <- as.integer(as.logical(fct$Myelination))
fct$Omnivore <- as.integer(as.logical(fct$Omnivore))
fct$Carnivore <- as.integer(as.logical(fct$Carnivore))
fct$Herbivore <- as.integer(as.logical(fct$Herbivore))
fct$Detritivore <- as.integer(as.logical(fct$Detritivore))
fct$Current <- as.integer(as.logical(fct$Current))
fct$Cruise <- as.integer(as.logical(fct$Cruise))
fct$Ambush <- as.integer(as.logical(fct$Ambush))

# Change size colname
colnames(fct)[7] <- "size"

# Compute Gower's distance matrix, with all species having 0 or just 1 NA and then just 0 NA
gow <- gowdis(fct[,c(7:9,11:14,16:18)])

dat_gow <- gow

save(gow, file = "gow_withNA.RData")
