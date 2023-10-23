rm(list = ls())

setwd("~/Documents/lab/araucaria_niches")

library(readr)
library(dplyr)
library(CoordinateCleaner)
library(raster)
library(ape)
library(stringi)

#envar <- raster::getData("worldclim", var = "tmin", res = 2.5)
#envar <- raster::getData("worldclim", var = "bio", res = 2.5)

tree <- read.tree("data/stull&al2021_supermatrix_dated.tre")
tax <- read.csv("data/classification.csv", sep = "\t") #WFO 2023

path_envs <- list.files(path = "wc2-5", pattern='.bil$', recursive = TRUE,
                        full.names = TRUE)
envsdt <- stack(path_envs)

cyc_data <- readr::read_tsv("data/Cycadopsida.txt")
gin_data <- readr::read_tsv("data/Ginkgoopsida.txt")
gne_data <- readr::read_tsv("data/Gnetopsida.txt")
pin_data <- readr::read_tsv("/Volumes/Personal/lab/araucaria_niches/data/Pinopsida.txt")

cyc_clean <- cyc_data  %>%
  setNames(tolower(names(.))) %>%
  filter(coordinateprecision < 0.01 | is.na(coordinateprecision)) %>% 
  #filter(!basisofrecord %in% c("FOSSIL_SPECIMEN","LIVING_SPECIMEN", "PRESERVED_SPECIMEN")) %>%
  filter(basisofrecord %in% c("HUMAN_OBSERVATION")) %>%
  filter(!coordinateuncertaintyinmeters %in% c(301,3036,999,9999)) %>% 
  filter(decimallatitude != 0 | decimallongitude != 0) %>%
  #filter(is.na(as.factor(nwords(locality)))) %>%
  cc_cen(buffer = 2000) %>% # remove country centroids within 2km 
  cc_cap(buffer = 2000) %>% # remove capitals centroids within 2km
  cc_inst(buffer = 2000) %>% # remove zoo and herbaria within 2km 
  cc_sea() %>% # remove from ocean
  distinct(decimallongitude, decimallatitude, specieskey,datasetkey, .keep_all = TRUE) %>%
  glimpse() 

cyc_clean$species <- str_replace_all(cyc_clean$species, " ", "_")

gin_clean <- gin_data  %>%
  setNames(tolower(names(.))) %>%
  filter(coordinateprecision < 0.01 | is.na(coordinateprecision)) %>% 
  #filter(!basisofrecord %in% c("FOSSIL_SPECIMEN","LIVING_SPECIMEN", "PRESERVED_SPECIMEN")) %>%
  filter(basisofrecord %in% c("HUMAN_OBSERVATION")) %>%
  filter(!coordinateuncertaintyinmeters %in% c(301,3036,999,9999)) %>% 
  filter(decimallatitude != 0 | decimallongitude != 0) %>%
  #filter(is.na(as.factor(nwords(locality)))) %>%
  cc_cen(buffer = 2000) %>% # remove country centroids within 2km 
  cc_cap(buffer = 2000) %>% # remove capitals centroids within 2km
  cc_inst(buffer = 2000) %>% # remove zoo and herbaria within 2km 
  cc_sea() %>% # remove from ocean
  distinct(decimallongitude, decimallatitude, specieskey,datasetkey, .keep_all = TRUE) %>%
  glimpse() 

gin_clean$species <- stri_replace_all_fixed(gin_clean$species, " ", "_")

gne_clean <- gne_data  %>%
  setNames(tolower(names(.))) %>%
  filter(coordinateprecision < 0.01 | is.na(coordinateprecision)) %>% 
  #filter(!basisofrecord %in% c("FOSSIL_SPECIMEN","LIVING_SPECIMEN", "PRESERVED_SPECIMEN")) %>%
  filter(basisofrecord %in% c("HUMAN_OBSERVATION")) %>%
  filter(!coordinateuncertaintyinmeters %in% c(301,3036,999,9999)) %>% 
  filter(decimallatitude != 0 | decimallongitude != 0) %>%
  #filter(is.na(as.factor(nwords(locality)))) %>%
  cc_cen(buffer = 2000) %>% # remove country centroids within 2km 
  cc_cap(buffer = 2000) %>% # remove capitals centroids within 2km
  cc_inst(buffer = 2000) %>% # remove zoo and herbaria within 2km 
  cc_sea() %>% # remove from ocean
  distinct(decimallongitude, decimallatitude, specieskey,datasetkey, .keep_all = TRUE) %>%
  glimpse() 

gne_clean$species <- stri_replace_all_fixed(gne_clean$species, " ", "_")

pin_clean <- pin_data  %>%
  setNames(tolower(names(.))) %>%
  filter(coordinateprecision < 0.01 | is.na(coordinateprecision)) %>% 
  #filter(!basisofrecord %in% c("FOSSIL_SPECIMEN","LIVING_SPECIMEN", "PRESERVED_SPECIMEN")) %>%
  filter(basisofrecord %in% c("HUMAN_OBSERVATION")) %>%
  filter(!coordinateuncertaintyinmeters %in% c(301,3036,999,9999)) %>% 
  filter(decimallatitude != 0 | decimallongitude != 0) %>%
  #filter(is.na(as.factor(nwords(locality)))) %>%
  cc_cen(buffer = 2000) %>% # remove country centroids within 2km 
  cc_cap(buffer = 2000) %>% # remove capitals centroids within 2km
  cc_inst(buffer = 2000) %>% # remove zoo and herbaria within 2km 
  cc_sea() %>% # remove from ocean
  distinct(decimallongitude, decimallatitude, specieskey,datasetkey, .keep_all = TRUE) %>%
  glimpse() 

pin_clean$species <- stri_replace_all_fixed(pin_clean$species, " ", "_")

### creating a SP object to subset env data
cyc_clean_sp <- cyc_clean
gin_clean_sp <- gin_clean
gne_clean_sp <- gne_clean
pin_clean_sp <- pin_clean

coordinates(cyc_clean_sp) <- ~decimallongitude+decimallatitude
coordinates(gin_clean_sp) <- ~decimallongitude+decimallatitude
coordinates(gne_clean_sp) <- ~decimallongitude+decimallatitude
coordinates(pin_clean_sp) <- ~decimallongitude+decimallatitude

crs(cyc_clean_sp) <- crs(envsdt)
crs(gin_clean_sp) <- crs(envsdt)
crs(gne_clean_sp) <- crs(envsdt)
crs(pin_clean_sp) <- crs(envsdt)

# extracting env data
env_cyc <- raster::extract(envsdt, cyc_clean_sp, fun = mean, na.rm = T, sp = T)
env_gin <- raster::extract(envsdt, gin_clean_sp, fun = mean, na.rm = T, sp = T)
env_gne <- raster::extract(envsdt, gne_clean_sp, fun = mean, na.rm = T, sp = T)
env_pin <- raster::extract(envsdt, pin_clean_sp, fun = mean, na.rm = T, sp = T)

env_cyc_df <- data.frame(env_cyc)
env_gin_df <- data.frame(env_gin)
env_gne_df <- data.frame(env_gne)
env_pin_df <- data.frame(env_pin)

## cleaning the data from extra stuff
env_cyc_df_clean <- env_cyc_df %>% 
  dplyr::select(species, decimallongitude, decimallatitude, bio1:tmin9)

env_gin_df_clean <- env_gin_df %>% 
  dplyr::select(species, decimallongitude, decimallatitude, bio1:tmin9)

env_gne_df_clean <- env_gne_df %>% 
  dplyr::select(species, decimallongitude, decimallatitude, bio1:tmin9)

env_pin_df_clean <- env_pin_df %>% 
  dplyr::select(species, decimallongitude, decimallatitude, bio1:tmin9)

#save(env_cyc_df_clean, file = "data/env_cyc_df_clean.RData")
#save(env_gin_df_clean, file = "data/env_gin_df_clean.RData")
#save(env_gne_df_clean, file = "data/env_gne_df_clean.RData")
#save(env_pin_df_clean, file = "data/env_pin_df_clean.RData")

#load(file = "data/env_cyc_df_clean.RData")
#load(file = "data/env_gin_df_clean.RData")
#load(file = "data/env_gne_df_clean.RData")
#load(file = "data/env_pin_df_clean.RData")

#checking taxon sampling

tax$scientificName <- stri_replace_all_fixed(tax$scientificName, " ", "_")

tax_gym <- tax[tax$majorGroup == "G" & 
               tax$taxonRank == "species" & 
               tax$taxonomicStatus == "Accepted", ] # 1204 spp

spp_cyc <- tax$scientificName[tax$majorGroup == "G" & 
                              tax$taxonRank == "species" &
                              tax$taxonomicStatus == "Accepted" &
                              tax$family %in% c("Cycadaceae", "Zamiaceae")] # 356 spp

spp_gin <- tax$scientificName[tax$majorGroup == "G" & 
                              tax$taxonRank == "species" &
                              tax$taxonomicStatus == "Accepted" &
                              tax$family == "Ginkgoaceae"] # 1 spp

spp_gne <- tax$scientificName[tax$majorGroup == "G" & 
                              tax$taxonRank == "species" &
                              tax$taxonomicStatus == "Accepted" &
                              tax$family %in% c("Ephedraceae",
                                                "Welwitschiaceae",
                                                "Gnetaceae")] # 1 spp

spp_pin <- tax$scientificName[tax$majorGroup == "G" & 
                              tax$taxonRank == "species" &
                              tax$taxonomicStatus == "Accepted" &
                              tax$family %in% c("Araucariaceae",
                                                "Cupressaceae",
                                                "Pinaceae",
                                                "Podocarpaceae",
                                                "Sciadopityaceae",
                                                "Taxaceae")] # 847 spp

cyc_spp_acc <- unique(env_cyc_df_clean$species[env_cyc_df_clean$species %in% 
                                                 spp_cyc])
cyc_spp_acc <- cyc_spp_acc[complete.cases(cyc_spp_acc)]

dif_cyc <- setdiff(unique(env_cyc_df_clean$species), cyc_spp_acc)
dif_cyc <- dif_cyc[complete.cases(dif_cyc)]
dif_cyc <- cbind(dif_cyc, NA)

for (i in 1:nrow(dif_cyc)) {
  ID <- tax_gym$acceptedNameUsageID[tax_gym$scientificName == dif_cyc[i, 1]]
  syn <- tax_gym$scientificName[tax_gym$taxonID == ID & tax_gym$taxonomicStatus == "Accepted"]
  dif_cyc[i, 2] <- ifelse(length(syn) == 0, NA, 
                          paste(strsplit(syn, "_")[[1]][c(1, 2)], collapse = "_"))
}
dif_cyc <- dif_cyc[complete.cases(dif_cyc[, 2]), ]

for (i in 1:length(env_cyc_df_clean$species)) {
  if (env_cyc_df_clean$species[i] %in% dif_cyc[, 1]) {
    env_cyc_df_clean$species[i] <- dif_cyc[, 2][dif_cyc[, 1] == env_cyc_df_clean$species[i]]
  } 
}
cyc_spp_acc <- unique(env_cyc_df_clean$species[env_cyc_df_clean$species %in% 
                                                 spp_cyc])
cyc_spp_acc <- cyc_spp_acc[complete.cases(cyc_spp_acc)]

gen_cyc <- unique(stri_split_fixed(spp_cyc, "_", simplify = T)[, 1])

pro_cyc <- matrix(ncol = 3, nrow = length(gen_cyc))
for (i in 1:length(gen_cyc)) {
  spp <- tax_gym$scientificName[tax_gym$genus == gen_cyc[i]]
  spp_dat <- cyc_spp_acc[cyc_spp_acc %in% spp]
  pro_cyc[i, 1] <- length(spp)
  pro_cyc[i, 2] <- length(spp_dat)
  pro_cyc[i, 3] <- (length(spp_dat) * 100)/length(spp)
}
colnames(pro_cyc) <- c("tax", "data", "prop")
rownames(pro_cyc) <- gen_cyc

gne_spp_acc <- unique(env_gne_df_clean$species[env_gne_df_clean$species %in% 
                                                 spp_gne])
gne_spp_acc <- gne_spp_acc[complete.cases(gne_spp_acc)]

dif_gne <- setdiff(unique(env_gne_df_clean$species), gne_spp_acc)
dif_gne <- dif_gne[complete.cases(dif_gne)]
dif_gne <- cbind(dif_gne, NA)

for (i in 1:nrow(dif_gne)) {
  ID <- tax_gym$acceptedNameUsageID[tax_gym$scientificName == dif_gne[i, 1]]
  syn <- tax_gym$scientificName[tax_gym$taxonID == ID & tax_gym$taxonomicStatus == "Accepted"]
  dif_gne[i, 2] <- ifelse(length(syn) == 0, NA, 
                          paste(strsplit(syn, "_")[[1]][c(1, 2)], collapse = "_"))
}
dif_gne <- dif_gne[complete.cases(dif_gne[, 2]), ]

for (i in 1:length(env_gne_df_clean$species)) {
  if (env_gne_df_clean$species[i] %in% dif_gne[, 1]) {
    env_gne_df_clean$species[i] <- dif_gne[, 2][dif_gne[, 1] == env_gne_df_clean$species[i]]
  } 
}
gne_spp_acc <- unique(env_gne_df_clean$species[env_gne_df_clean$species %in% 
                                                 spp_gne])
gne_spp_acc <- gne_spp_acc[complete.cases(gne_spp_acc)]

gen_gne <- unique(stri_split_fixed(spp_gne, "_", simplify = T)[, 1])

pro_gne <- matrix(ncol = 3, nrow = length(gen_gne))
for (i in 1:length(gen_gne)) {
  spp <- tax_gym$scientificName[tax_gym$genus == gen_gne[i]]
  spp_dat <- gne_spp_acc[gne_spp_acc %in% spp]
  pro_gne[i, 1] <- length(spp)
  pro_gne[i, 2] <- length(spp_dat)
  pro_gne[i, 3] <- (length(spp_dat) * 100)/length(spp)
}
colnames(pro_gne) <- c("tax", "data", "prop")
rownames(pro_gne) <- gen_gne

pin_spp_acc <- unique(env_pin_df_clean$species[env_pin_df_clean$species %in% 
                                                 spp_pin])
pin_spp_acc <- pin_spp_acc[complete.cases(pin_spp_acc)]

dif_pin <- setdiff(unique(env_pin_df_clean$species), pin_spp_acc)
dif_pin <- dif_pin[complete.cases(dif_pin)]
dif_pin <- cbind(dif_pin, NA)

for (i in 1:nrow(dif_pin)) {
  ID <- tax_gym$acceptedNameUsageID[tax_gym$scientificName == dif_pin[i, 1]]
  syn <- tax_gym$scientificName[tax_gym$taxonID == ID & tax_gym$taxonomicStatus == "Accepted"]
  dif_pin[i, 2] <- ifelse(length(syn) == 0, NA, 
                          paste(strsplit(syn, "_")[[1]][c(1, 2)], collapse = "_"))
}
dif_pin <- dif_pin[complete.cases(dif_pin[, 2]), ]

for (i in 1:length(env_pin_df_clean$species)) {
  if (env_pin_df_clean$species[i] %in% dif_pin[, 1]) {
    env_pin_df_clean$species[i] <- dif_pin[, 2][dif_pin[, 1] == env_pin_df_clean$species[i]]
  } 
}
pin_spp_acc <- unique(env_pin_df_clean$species[env_pin_df_clean$species %in% 
                                                 spp_pin])
pin_spp_acc <- pin_spp_acc[complete.cases(pin_spp_acc)]

gen_pin <- unique(stri_split_fixed(spp_pin, "_", simplify = T)[, 1])

pro_pin <- matrix(ncol = 3, nrow = length(gen_pin))
for (i in 1:length(gen_pin)) {
  spp <- tax_gym$scientificName[tax_gym$genus == gen_pin[i]]
  spp_dat <- pin_spp_acc[pin_spp_acc %in% spp]
  pro_pin[i, 1] <- length(spp)
  pro_pin[i, 2] <- length(spp_dat)
  pro_pin[i, 3] <- (length(spp_dat) * 100)/length(spp)
}
colnames(pro_pin) <- c("tax", "data", "prop")
rownames(pro_pin) <- gen_pin

## maintaing in the dataset only those species also present in the ncbi data

all_spp <- c(unique(env_cyc_df_clean$species),
             unique(env_gin_df_clean$species),
             unique(env_gne_df_clean$species),
             unique(env_pin_df_clean$species))

spp_tr_cyc <- env_cyc_df_clean$species %in% intersect(env_cyc_df_clean$species,
                                                      tree$tip.label)
env_cyc_df_clean2 <- env_cyc_df_clean[spp_tr_cyc,] %>% na.omit() %>% droplevels()

spp_tr_gin <- env_gin_df_clean$species %in% intersect(env_gin_df_clean$species,
                                                      tree$tip.label) 
env_gin_df_clean2 <- env_gin_df_clean[spp_tr_gin,] %>% na.omit() %>% droplevels()

spp_tr_pin <- env_pin_df_clean$species %in% intersect(env_pin_df_clean$species,
                                                      tree$tip.label) 
env_pin_df_clean2 <- env_pin_df_clean[spp_tr_pin,] %>% na.omit() %>% droplevels()

#write.csv(env_cyc_df_clean2, "data/Environmental_Data_Cyc.csv")
#write.csv(env_gin_df_clean2, "data/Environmental_Data_Gin.csv")
#write.csv(env_pin_df_clean2, "data/Environmental_Data_Pin.csv")

env_cyc_df_clean2 <- read.csv("data/Environmental_Data_Cyc.csv",
                              stringsAsFactors = TRUE, header = TRUE) %>%
  dplyr::select(-X) %>% na.omit() %>% droplevels()

env_gin_df_clean2 <- read.csv("data/Environmental_Data_Gin.csv",
                              stringsAsFactors = TRUE, header = TRUE) %>%
  dplyr::select(-X) %>% na.omit() %>% droplevels()

env_pin_df_clean2 <- read.csv("data/Environmental_Data_Pin.csv",
                              stringsAsFactors = TRUE, header = TRUE) %>%
  dplyr::select(-X) %>% na.omit() %>% droplevels()

env_df_clean2 <- rbind(env_cyc_df_clean2, env_gin_df_clean2, env_pin_df_clean2)

#save(env_df_clean2, file = "data/env_df_clean2.RData")
#load(file = "data/env_df_clean2.RData")

#### ~~~~ Calculating Hypervolumes and holes therein using persistence homology ~~~~ ####
#### ~~~~~   Hypervolumes 
library(hypervolume)
library(tidyr)
library(purrr)
library(TDAstats)
library(ggplot2)

env_df_clean_pca2 <- prcomp(env_df_clean2[, 4:34])
env_df_clean_pca2 <- as.data.frame(
                            cbind(as.character(env_df_clean2$species),
                                  env_df_clean2$decimallongitude,
                                  env_df_clean2$decimallatitude,
                                  env_df_clean_pca2$x[, 1:3]))
colnames(env_df_clean_pca2) <- c("species", "lon", "lat", "PC1", "PC2", "PC3")
env_df_clean_pca2[, 2:6] <- apply(env_df_clean_pca2[, 2:6], 2, as.numeric)
keep_cyc <- names(table(env_df_clean_pca2$species)[table(env_df_clean_pca2$species) > 2])
keep_cyc <- keep_cyc[keep_cyc != "Zamia_purpurea"]
keep_cyc <- keep_cyc[keep_cyc != "Zamia_chigua"]
env_df_clean_pca2 <- subset(env_df_clean_pca2, subset = species %in% keep_cyc)

## Calculating hypervolumes ##
#8:59
env_df_clean_hypervolume <- env_df_clean_pca2 %>%
  nest(-species) %>%
  mutate(hypervol = map(data, function(.x){
    hyper <- hypervolume_gaussian(.x[3:5])
    hyper
  }))
b <- Sys.time()
#### ~~ Persistence Homology calculations ~~~ ####
## Functions 
subsample.distance_demo <- function (x, size, d, d.max = NULL, 
                                     replacement = FALSE, latlong = FALSE,
                                     echo = FALSE) {
  
  if (missing(x))
    
    stop("Must define a spatial object")
  
  if (missing(d))
    
    stop("Must define minimum separation distance")
  
  if (!is.null(d.max)) {
    
    if (d.max <= d)
      
      stop("Maximum distance must be larger than minimum")
    
  }
  
  if (!any(class(x)[1] == c("SpatialPointsDataFrame", "SpatialPolygonsDataFrame")))
    
    stop("x must be sp class polygons or points")
  
  if (latlong == TRUE) {
    
    message("geographic projection distances must be in kilometers")
    
  }
  
  if (size >= nrow(x))
    
    stop("subsample size must be smaller than population")
  
  rs <- sample(1:nrow(x), 1)
  
  s <- x[rs, ]
  
  if (replacement == FALSE) {
    
    x <- x[-rs, ]
    
  }
  
  deval = TRUE
  
  for (i in 2:size) {
    
    nsamp = 0
    
    while (deval == TRUE) {
      
      rs <- sample(1:nrow(x), 1)
      
      pts.dist = sp::spDists(s, x[rs, ], longlat = latlong)
      
      if (is.null(d.max)) {
        
        deval <- any(pts.dist < d, na.rm = TRUE)
        
      }
      
      else {
        
        deval <- any(pts.dist < d, na.rm = TRUE) | any(pts.dist >
                                                         
                                                         d.max, na.rm = TRUE)
        
      }
      
      nsamp = nsamp + 1
      
      if (echo)
        
        cat("Sample iteration=", nsamp, "\n")
      
      if (nsamp == nrow(x))
        
        break
      
    }
    
    if (echo) {
      
      cat("\n", "Min distance for", i, "=", min(pts.dist,
                                                
                                                na.rm = TRUE), "\n")
      
      cat(" Max distance for", i, "=", max(pts.dist, na.rm = TRUE),
          
          "\n")
      
    }
    
    if (nsamp == nrow(x)) {
      
      message(paste0("Warning: sampling cannot converge at n=",
                     
                     size, " returning n=", nrow(s)))
      
      return(s)
      
    }
    
    deval = TRUE
    
    s <- rbind(s, x[rs, ])
    
    if (replacement == FALSE) {
      
      x <- x[-rs, ]
      
    }
    
  }
  
  return(s)
  
}

subsample_distance <- function(df, size, d, ...) {
  
  sub.meuse <- c()
  
  for(i in 1:ncol(df)){
    
    for(j in 1:ncol(df)){
      
      if(i == j | j > i){
        
        next
        
      } else {
        
        testsp <- df
        
        sp::coordinates(testsp) <- c(i, j)
        
        sub <- subsample.distance_demo(testsp, size = size/ncol(df), d = d) ##d is for the minimum distance between the points sampled
        
        sub.meuse <- bind_rows(sub.meuse, as.data.frame(sub))
        
      }
      
    }
    
  }
  
  return(sub.meuse)
  
}

####
###### ~~~~~ automating the calculations for PH ~~~~~~~ #####
env_cyc_df_clean_hypervolume2 <- env_cyc_df_clean_hypervolume %>%
  mutate(randomHyper = map(hypervol, function(.x){
    out <- apply(.x@RandomPoints, 2, base::scale)
  })) %>%
  mutate(persistenceHomology = map(randomHyper, function(.x){
    set.seed(145) ## setting seed for reproducibility
    subset_hyper <- as.data.frame(.x) 
    subset_hyper <- subsample_distance(subset_hyper,
                                       size = 300,
                                       d = 0.2) ## subseting  
    PH_calc <- calculate_homology(subset_hyper, dim = 2) # calculating ph
    return(list(hyper_subset = subset_hyper,
                persHomol = as.data.frame(PH_calc)))
  }))

env_cyc_df_clean_hypervolume2 <- env_cyc_df_clean_hypervolume %>%
  mutate(randomHyper = map(hypervol, function(.x){
    out <- apply(.x@RandomPoints, 2, base::scale)
  })) %>%
  mutate(persistenceHomology = map(randomHyper, function(.x){
    set.seed(145) ## setting seed for reproducibility
    subset_hyper <- as.data.frame(.x) 
    subset_hyper <- subsample_distance(subset_hyper,
                                       size = 300,
                                       d = 0.2) ## subseting  
    PH_calc <- calculate_homology(subset_hyper, dim = 2) # calculating ph
    return(list(hyper_subset = subset_hyper,
                persHomol = as.data.frame(PH_calc)))
  }))

env_cyc_df_clean_hypervolume3 <- env_cyc_df_clean_hypervolume2 %>%
  mutate(hull = map(randomHyper, function(.x){
    hullout <- with(as.data.frame(.x), chull(PC1, PC2))
    hullout
  }),
  out = map2(randomHyper, hull, ~ .x[.y,,drop=FALSE]))

#save(env_cyc_df_clean_hypervolume3, file = "data/env_araucaria_df_clean_hypervolume3.RData")

tree_cyc <- read.nexus("data/Cycadales_Condamine_2015.tre")
index <- which(env_cyc_df_clean_hypervolume3$species %in% tree_cyc$tip.label)
env_cyc_df_clean_hypervolume3 <- env_cyc_df_clean_hypervolume3[index,]

## Plots of the point clouds
plot_points <- list()
for(i in 1:nrow(env_cyc_df_clean_hypervolume3)){
  plot_points[[i]] <- ggplot(as.data.frame(env_cyc_df_clean_hypervolume3$out[[i]]), 
                             aes(PC1, PC2, col = PC3)) +
    geom_polygon(alpha = 0.05, fill = "royalblue2", col = NA) +
    geom_point(data = as.data.frame(env_cyc_df_clean_hypervolume3$randomHyper[[i]]),
               aes(PC1, PC2, col = PC3), inherit.aes = FALSE, alpha = 0.07, size = 0.2) +
    geom_point(data = env_cyc_df_clean_hypervolume3$persistenceHomology[[i]]$hyper_subset,
               aes(PC1, PC2, col = PC3), 
               inherit.aes = FALSE, alpha = 1, size = 1.5) +
    ggtitle(paste0(env_cyc_df_clean_hypervolume3$species[i])) +
    scale_color_viridis_c() +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 9),
          plot.title = element_text(size = 9, face = "bold.italic", hjust = 0.5))
  
}

## plots of the persistence diagrams ##
plot_persistence <- list()
for(i in 1:nrow(env_cyc_df_clean_hypervolume3)){
  plot_persistence[[i]] <- ggplot(data = env_cyc_df_clean_hypervolume3$persistenceHomology[[i]]$persHomol,
                                  aes(x = birth, 
                                      y = death, 
                                      shape = as.factor(dimension), 
                                      color = as.factor(dimension))) +
    geom_point(alpha = 0.9) +
    geom_abline(slope = 1, size = 0.2) +
    xlab('Birth') +
    ylab('Death') +
    theme_linedraw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 9),
          axis.title = element_text(size = 11),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 9)) +
    scale_color_manual('Dimension', values = c("darkseagreen3", "royalblue2", "hotpink")) +
    scale_shape_manual('Dimension', values=c(15,17,19))
  
}

library(patchwork)
final_plot <- c(rbind(plot_points, plot_persistence))
wrap_plots(final_plot) + plot_layout(ncol = 4)



### ~~~~ Pairwise bottleneck distances of the persistence homology ~~~~~ ###
### Estimating pairwise bottleneck distances from the PHs of the species ###
## Function to estimate the distance between two persistence diagrams
homology_distance <- function(x, y, dim = 2, ...){
  output <- c()
  for(i in 0:dim){
    dist <- TDA::bottleneck(as.matrix(x), as.matrix(y), i)
    output <- cbind(output, dist)
    colnames(output)[i + 1] <- paste(i)
  }
  return(as.data.frame(output))
}


## Pairwise distance of the persistence (bottleneck distances)
pairwise_persistence <- function(df){
  output <- c()
  ## for loop for pairwise interactions ##
  for(i in 1:nrow(df)){
    for(j in 1:nrow(df)){
      if(i <= j){ #excluding self-comparison and double-comparison
        next
      }else{
        pairwise_dist_out <- homology_distance(df$persistenceHomology[[i]]$persHomol,
                                               df$persistenceHomology[[j]]$persHomol)
        sp1 <- df$species[j]
        sp2 <- df$species[i]
        
        outdata <- data.frame(species1 = sp1, species2 = sp2, dist = pairwise_dist_out)
        output <- as.data.frame(bind_rows(output, outdata))
      }
    }
  }
  
  output_final <- data.frame(species1 = output$species2, 
                             species2 = output$species1) %>%
    bind_cols(., output[3:5]) %>%
    bind_rows(., output)
  
  return(as.data.frame(output_final))
}




homology_distance(env_araucaria_df_clean_hypervolume3$persistenceHomology[[1]]$persHomol,
                  env_araucaria_df_clean_hypervolume3$persistenceHomology[[2]]$persHomol)


### This works but needs to change the name of the object from "test" to someting permanent
test <- pairwise_persistence(env_araucaria_df_clean_hypervolume3)
colnames(test) <- c("species1", "species2", "dist.0", "dist.1", "dist.2")



plot_dist0 <- ggplot(test, aes(species1, species2, fill= dist.0)) + 
  geom_tile() +
  scale_fill_viridis_c('Distance\n(Dim 0)', option = "mako") + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = -0.0)) +
  xlab('') +  scale_x_discrete(drop = FALSE) +
  ylab('Species 2')


plot_dist1 <- ggplot(test, aes(species1, species2, fill= dist.1)) + 
  geom_tile() +
  scale_fill_viridis_c('Distance\n(Dim 1)', option = "mako") + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = -0.0),
        axis.text.y = element_blank()) +
  xlab('Species 1') + 
  ylab('')


plot_dist2 <- ggplot(test, aes(species1, species2, fill= dist.2)) + 
  geom_tile() +
  scale_fill_viridis_c('Distance\n(Dim 2)', option = "mako") + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = -0.0),
        axis.text.y = element_blank()) +
  xlab('') + 
  ylab('')

library(patchwork)
final_plot_distance <-wrap_plots(plot_dist0, plot_dist1, plot_dist2)
final_plot_distance






## Summary dimensions
summary_persistence <- function(df){
  output <- c()
  ## for loop for pairwise interactions ##
  for(i in 1:nrow(df)){
    outdata <- df$persistenceHomology[[i]]$persHomol %>%
      mutate(surv = death - birth) %>%
      group_by(dimension) %>%
      summarise(avg_surv = mean(surv),
                sd_surv = sd(surv),
                max_surv = max(surv)) %>%
      mutate(species = df$species[i])
    output <- as.data.frame(bind_rows(output, outdata))
  }
  output <- output[c(5,1,2,3,4)]
  return(as.data.frame(output))
}


## calculating average survival of connected components
persistence_info_dt <- summary_persistence(env_araucaria_df_clean_hypervolume3) 
#

volume_hyper <-  env_araucaria_df_clean_hypervolume3 %>%
  mutate(hv_vol = purrr::map(hypervol, function(.x){
    .x@Volume
  })) %>%
  dplyr::select(species, hv_vol) %>%
  tidyr::unnest(hv_vol) %>%
  right_join(., persistence_info_dt) 

## subsetting only species in the tree
index_vol <- which(volume_hyper$species %in% araucariaTree_rooted$tip.label)
volume_hyper <- volume_hyper[index_vol,]


outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

araucariaTree_cleaned <- drop.tip(treeRatchet, as.vector(outersect(treeRatchet$tip.label, volume_hyper$species)))
araucariaTree_cleaned$edge.length <- ifelse(araucariaTree_cleaned$edge.length == 0,
                                            2000,
                                            2000+araucariaTree_cleaned$edge.length) ## I cannot have branches with zero length or polytomies, so using this to increase branch sizes
plot(araucariaTree_cleaned)

### ±±± Plotting the hypervolume for the main text
## plotting phylogeny with continuous traits
volume_hyper_dim0 <- volume_hyper %>%
  filter(dimension == 2) %>%
  droplevels() %>%
  tibble::column_to_rownames(., "species")

## Volume of hypervolume
volume_hyper_vol_plot <- setNames(volume_hyper_dim0[,1], rownames(volume_hyper_dim0))
obj1<-contMap(araucariaTree_cleaned,volume_hyper_vol_plot,plot=FALSE)
obj1<-setMap(obj1,invert=TRUE)
plot(obj1,fsize=c(0.9,0.5),outline=TRUE,lwd=c(4,7),leg.txt="Volume")

## dimension (Plot of the dimensions, done one by one (note: needs to change the "filter" dimension to obtain the others, as above))
volume_hyper_dim0_plot <- setNames(volume_hyper_dim0[,5], rownames(volume_hyper_dim0))
obj2<-contMap(araucariaTree_cleaned,volume_hyper_dim0_plot,plot=FALSE)
obj2<-setMap(obj2,invert=TRUE)
plot(obj2,fsize=c(0.9,0.5),outline=TRUE,lwd=c(4,7),leg.txt="Dim 2")


fancydt <- volume_hyper %>% dplyr::select(-sd_surv) %>%
  pivot_wider(., values_from = avg_surv, names_from = dimension) %>%
  tibble::column_to_rownames(., "species")

fancyPlot<-fancyTree(araucariaTree_cleaned,type="scattergram",X=fancydt[,c(1:4)])

signalphy <- phylosig(araucariaTree_cleaned, volume_hyper_dim0_plot, method = "K",
                      test = TRUE, nsim = 999)
signalphy$sim.K
str(signalphy)
plot(signalphy)


calculate_physig <- function(tree, 
                             data,
                             names = NULL,
                             method = "K",
                             test = TRUE, 
                             nsim = 999, 
                             dim = 2, ...){
  output <- c()
  inner <- c()
  for(j in 0:dim){
    volume_hyper <- data %>%
      filter(dimension == j) %>%
      droplevels() %>%
      tibble::column_to_rownames(., "species") ## subsetting the data and dim
    
    volume_hyper <- setNames(volume_hyper[,1], rownames(volume_hyper))
    signalphy <- phylosig(tree, 
                          volume_hyper, method = method,
                          test = test, nsim = nsim)# physignal calculation
    inner <-  bind_rows(inner, 
                        data.frame(dim = j,
                                   simK = signalphy$sim.K,
                                   pval = signalphy$P))
  }
  output <- bind_rows(output, inner)
  return(output)
}



hypervolume_sig <- calculate_physig(araucariaTree_cleaned, volume_hyper[c(1,2,3)]) %>%
  mutate(type = "Hypervolume")

avgdist_sig <- calculate_physig(araucariaTree_cleaned, volume_hyper[c(1,4,3)]) %>%
  mutate(type = "Mean dist")

sddist_sig <- calculate_physig(araucariaTree_cleaned, volume_hyper[c(1,5,3)]) %>%
  mutate(type = "SD dist")

maxdist_sig <- calculate_physig(araucariaTree_cleaned, volume_hyper[c(1,6,3)]) %>%
  mutate(type = "Max dist")

final_signal <- bind_rows(hypervolume_sig,
                          avgdist_sig,
                          sddist_sig,
                          maxdist_sig)
tail(avgdist_sig)

dummy <- final_signal  %>%
  group_by(dim, type) %>%
  summarize(mean = mean(simK),
            pvalue = mean(pval))


signal_plot <- ggplot(data =final_signal, 
                      aes(x = simK, col = type, fill = type)) +
  geom_density(alpha = 0.4)+
  facet_grid(dim ~ type) +
  geom_vline(data = dummy, aes(xintercept = mean, color = type), linetype = "dashed") + 
  geom_text(data = dummy, aes(label=paste0('p=', round(pvalue, 3))), 
            x = 2.5, y = 3, hjust=1, vjust=1,size = 3.5,
            inherit.aes = FALSE) + 
  theme_linedraw() + 
  scale_color_manual('Variable',
                     breaks = c("Hypervolume", "Mean dist", "SD dist", "Max dist"),
                     values = c("Purple4", "Firebrick3", "tomato4", "royalblue2"))+
  scale_fill_manual('Variable',
                    breaks = c("Hypervolume", "Mean dist", "SD dist", "Max dist"),
                    values = c("Purple4", "Firebrick3", "tomato4","royalblue2")) + 
  theme(panel.grid = element_blank(),
        strip.background.x = element_rect(fill = 'grey80'),
        strip.background.y = element_rect(fill = 'grey95'),
        strip.text = element_text(colour = "black")) +
  xlab('Blomberg\'s K\n(Simulated)') + 
  ylab('Density')
signal_plot






## ~~~ Plotting heatmap with phylogeny ~~~~ ###
# build a ggplot with a geom_tree
araucariaTree2 <- root(treeRatchet, outgroup = "Acmopyle pancheri", resolve.root = TRUE)
final_tree <- ggtree::ggtree(araucariaTree2, 
                             branch.length = "none", 
                             col = "black") +
  ggtree::geom_tippoint(size = 0.5, alpha = 1, col = "black") + 
  #geom_text(aes(label=node), hjust=-.3) +
  ggtree:: geom_rootpoint(pch = 10, size= 0.4) +
  ggtree::geom_tiplab(size = 4, col = "black", size = 4) +
  ggtree::geom_nodelab(nudge_x = -0.6, nudge_y = 0.5, size = 2.1, col = "grey30") + 
  theme(plot.margin = unit(c(14,8,14,8), "mm")) + 
  xlim(0,18) 
#geom_strip(14, 15, barsize=1, color='gold2', 
#           label = "Root", offset=6.5, offset.text = 0.2)
final_tree
## Getting the order of the tips in the phylogeny
is_tip <- araucariaTree$edge[, 2] <= length(araucariaTree$tip.label)
ordered_tips <- araucariaTree$edge[is_tip, 2]
araucariaTree$tip.label[ordered_tips]


## reordering the plot of the heatmaps
head(test)
test$species1 <- reorder(test$species1, araucariaTree$tip.label[ordered_tips])


test$species1[match(araucariaTree$tip.label[ordered_tips], test$species1)]
test$species1[match(test$species1, araucariaTree$tip.label[ordered_tips])]

##below works but is not permanent
test2 <- test[order(match(test$species1, araucariaTree$tip.label[ordered_tips])),]
rownames(test2 ) <- 1:nrow(test2)
