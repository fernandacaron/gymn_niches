rm(list = ls())

setwd("~/Documents/lab/gymnosperm_niches")

library(readr)
library(dplyr)
library(CoordinateCleaner)
library(raster)
library(ape)
library(stringi)
library(hypervolume)
library(tidyr)
library(purrr)
library(TDAstats)
library(ggplot2)
library(patchwork)
library(phytools)

#envar <- raster::getData("worldclim", var = "tmin", res = 2.5)
#envar <- raster::getData("worldclim", var = "bio", res = 2.5)

tree <- read.tree("data/stull&al2021_supermatrix_dated.tre")
tax <- read.csv("data/classification.csv", sep = "\t") #WFO 2023
tax <- tax[tax$taxonRank == "species" & tax$taxonomicStatus == "Accepted", ]
tax$scientificName <- stri_replace_all_fixed(tax$scientificName, " ", "_")

path_envs <- list.files(path = "wc2-5", pattern='.bil$', recursive = TRUE,
                        full.names = TRUE)
envsdt <- stack(path_envs)

# DO NOT RUN! (too long) - skip to line 150
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

# creating a SP object to subset env data
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

# cleaning the data from extra stuff
env_cyc_df_clean <- env_cyc_df %>% 
  dplyr::select(species, decimallongitude, decimallatitude, bio1:tmin9)

env_gin_df_clean <- env_gin_df %>% 
  dplyr::select(species, decimallongitude, decimallatitude, bio1:tmin9)

env_gne_df_clean <- env_gne_df %>% 
  dplyr::select(species, decimallongitude, decimallatitude, bio1:tmin9)

env_pin_df_clean <- env_pin_df %>% 
  dplyr::select(species, decimallongitude, decimallatitude, bio1:tmin9)

#write.csv(env_cyc_df_clean, "data/Environmental_Data_Cyc.csv")
#write.csv(env_gin_df_clean, "data/Environmental_Data_Gin.csv")
#write.csv(env_gne_df_clean, "data/Environmental_Data_Gne.csv")
#write.csv(env_pin_df_clean, "data/Environmental_Data_Pin.csv")

env_cyc_df_clean2 <- read.csv("data/Environmental_Data_Cyc.csv",
                              stringsAsFactors = TRUE, header = TRUE) %>%
  dplyr::select(-X) %>% na.omit() %>% droplevels()

env_gin_df_clean2 <- read.csv("data/Environmental_Data_Gin.csv",
                              stringsAsFactors = TRUE, header = TRUE) %>%
  dplyr::select(-X) %>% na.omit() %>% droplevels()

env_gne_df_clean2 <- read.csv("data/Environmental_Data_Gne.csv",
                              stringsAsFactors = TRUE, header = TRUE) %>%
  dplyr::select(-X) %>% na.omit() %>% droplevels()

env_pin_df_clean2 <- read.csv("data/Environmental_Data_Pin.csv",
                              stringsAsFactors = TRUE, header = TRUE) %>%
  dplyr::select(-X) %>% na.omit() %>% droplevels()

env_df_clean2 <- rbind(env_cyc_df_clean2, env_gin_df_clean2, env_gne_df_clean2,
                       env_pin_df_clean2)

spp_tree <- env_df_clean2$species %in% intersect(env_df_clean2$species,
                                                 tree$tip.label)
env_df_clean3 <- env_df_clean2[spp_tree,] %>% na.omit() %>% droplevels()

fam_tax <- tax$family
names(fam_tax) <- tax$scientificName

gen_tax <- tax$genus
names(gen_tax) <- tax$scientificName

env_df_gen <- env_df_clean3 %>% 
  mutate(species = ifelse(species %in% names(gen_tax), 
                          gen_tax[as.character(species)],
                          NA)) %>% 
  na.omit() %>% 
  droplevels()
#save(env_df_gen, file = "data/env_df_gen.RData")
#load(file = "data/env_df_gen.RData")

env_df_fam <- env_df_clean3 %>% 
  mutate(species = ifelse(species %in% names(fam_tax), 
                          fam_tax[as.character(species)],
                          NA)) %>% 
  na.omit() %>% 
  droplevels()
#save(env_df_fam, file = "data/env_df_fam.RData")
#load(file = "data/env_df_fam.RData")

tr_gen <- tree
for (i in 1:length(tr_gen$tip.label)) {
  tr_gen$tip.label[i] <- gen_tax[tr_gen$tip.label[i]]
}
gen <- unique(env_df_gen$species)
tr_gen <- keep.tip(tr_gen, gen)

tr_fam <- tree
for (i in 1:length(tr_fam$tip.label)) {
  tr_fam$tip.label[i] <- fam_tax[tr_fam$tip.label[i]]
}
fam <- unique(env_df_fam$species)
tr_fam <- keep.tip(tr_fam, fam)

# calculating hypervolumes and holes therein using persistence homology

# hypervolumes 
env_df_gen_pca2 <- prcomp(env_df_gen[, 4:34])
env_df_gen_pca2 <- as.data.frame(
                            cbind(as.character(env_df_gen$species),
                                  env_df_gen$decimallongitude,
                                  env_df_gen$decimallatitude,
                                  env_df_gen_pca2$x[, 1:3]))
colnames(env_df_gen_pca2) <- c("genus", "lon", "lat", "PC1", "PC2", "PC3")
env_df_gen_pca2[, 2:6] <- apply(env_df_gen_pca2[, 2:6], 2, as.numeric)

env_df_fam_pca2 <- prcomp(env_df_fam[, 4:34])
env_df_fam_pca2 <- as.data.frame(
                            cbind(as.character(env_df_fam$species),
                                  env_df_fam$decimallongitude,
                                  env_df_fam$decimallatitude,
                                  env_df_fam_pca2$x[, 1:3]))
colnames(env_df_fam_pca2) <- c("family", "lon", "lat", "PC1", "PC2", "PC3")
env_df_fam_pca2[, 2:6] <- apply(env_df_fam_pca2[, 2:6], 2, as.numeric)

# calculating hypervolumes

## FALTOU: Larix, Abies, Picea, Juniperus, Pinus 

env_df_gen_hypervolume9 <- env_df_gen_pca2 %>%
  filter(!genus %in% c("Austrotaxus", "Falcatifolium", "Nothotsuga")) %>%
  filter(genus == "Larix") %>%
  nest(data = -genus) %>%
  mutate(hypervol = purrr::map(data, function(.x) {
    hyper <- hypervolume_gaussian(.x[3:5])
    hyper
  }))
save(env_df_gen_hypervolume9, file = "data/env_df_gen_hypervolume9.RData")

# DO NOT RUN, too long - skip to line 335 to load the finished dataset

# dividing per family
env_df_sci_hypervolume <- env_df_fam_pca2 %>%
  nest(data = -family) %>% filter(family == "Sciadopityaceae") %>%
  mutate(hypervol = map(data, function(.x) {
    hyper <- hypervolume_gaussian(.x[3:5])
    hyper
    }))
#save(env_df_sci_hypervolume, file = "data/env_df_sci_hypervolume.RData")
#load(file = "data/env_df_sci_hypervolume.RData")

env_df_gne_hypervolume <- env_df_fam_pca2 %>%
  nest(data = -family) %>% filter(family == "Gnetaceae") %>%
  mutate(hypervol = map(data, function(.x) {
    hyper <- hypervolume_gaussian(.x[3:5])
    hyper
  }))
#save(env_df_gne_hypervolume, file = "data/env_df_gne_hypervolume.RData")
#load(file = "data/env_df_gne_hypervolume.RData")

env_df_wel_hypervolume <- env_df_fam_pca2 %>%
  nest(data = -family) %>% filter(family == "Welwitschiaceae") %>%
  mutate(hypervol = map(data, function(.x) {
    hyper <- hypervolume_gaussian(.x[3:5])
    hyper
  }))
#save(env_df_wel_hypervolume, file = "data/env_df_wel_hypervolume.RData")
#load(file = "data/env_df_wel_hypervolume.RData")

env_df_cyc_hypervolume <- env_df_fam_pca2 %>%
  nest(data = -family) %>% filter(family == "Cycadaceae") %>%
  mutate(hypervol = map(data, function(.x) {
    hyper <- hypervolume_gaussian(.x[3:5])
    hyper
  }))
#save(env_df_cyc_hypervolume, file = "data/env_df_cyc_hypervolume.RData")
#load(file = "data/env_df_cyc_hypervolume.RData")

env_df_ara_hypervolume <- env_df_fam_pca2 %>%
  nest(data = -family) %>% filter(family == "Araucariaceae") %>%
  mutate(hypervol = map(data, function(.x) {
    hyper <- hypervolume_gaussian(.x[3:5])
    hyper
  }))
#save(env_df_ara_hypervolume, file = "data/env_df_ara_hypervolume.RData")
#load(file = "data/env_df_ara_hypervolume.RData")

env_df_zam_hypervolume <- env_df_fam_pca2 %>%
  nest(data = -family) %>% filter(family == "Zamiaceae") %>%
  mutate(hypervol = map(data, function(.x) {
    hyper <- hypervolume_gaussian(.x[3:5])
    hyper
  }))
#save(env_df_zam_hypervolume, file = "data/env_df_zam_hypervolume.RData")
#load(file = "data/env_df_zam_hypervolume.RData")

env_df_eph_hypervolume <- env_df_fam_pca2 %>%
  nest(data = -family) %>% filter(family == "Ephedraceae") %>%
  mutate(hypervol = map(data, function(.x) {
    hyper <- hypervolume_gaussian(.x[3:5])
    hyper
  }))
#save(env_df_eph_hypervolume, file = "data/env_df_eph_hypervolume.RData")
#load(file = "data/env_df_eph_hypervolume.RData")

env_df_gin_hypervolume <- env_df_fam_pca2 %>%
  nest(data = -family) %>% filter(family == "Ginkgoaceae") %>%
  mutate(hypervol = map(data, function(.x) {
    hyper <- hypervolume_gaussian(.x[3:5])
    hyper
  }))
#save(env_df_gin_hypervolume, file = "data/env_df_gin_hypervolume.RData")
#load(file = "data/env_df_gin_hypervolume.RData")

env_df_pod_hypervolume <- env_df_fam_pca2 %>%
  nest(data = -family) %>% filter(family == "Podocarpaceae") %>%
  mutate(hypervol = map(data, function(.x) {
    hyper <- hypervolume_gaussian(.x[3:5])
    hyper
  }))
#save(env_df_pod_hypervolume, file = "data/env_df_pod_hypervolume.RData")
#load(file = "data/env_df_pod_hypervolume.RData")

env_df_tax_hypervolume <- env_df_fam_pca2 %>%
  nest(data = -family) %>% filter(family == "Taxaceae") %>%
  mutate(hypervol = map(data, function(.x) {
    hyper <- hypervolume_gaussian(.x[3:5])
    hyper
  }))
#save(env_df_tax_hypervolume, file = "data/env_df_tax_hypervolume.RData")
#load(file = "data/env_df_tax_hypervolume.RData")

# TOO HEAVY - try to run on another place
env_df_cup_hypervolume <- env_df_fam_pca2 %>%
  nest(data = -family) %>% filter(family == "Cupressaceae") %>%
  mutate(hypervol = map(data, function(.x) {
    hyper <- hypervolume_gaussian(.x[3:5])
    hyper
  }))
#save(env_df_cup_hypervolume, file = "data/env_df_cup_hypervolume.RData")
#load(file = "data/env_df_cup_hypervolume.RData")

# TOO HEAVY - try to run on another place
env_df_pin_hypervolume <- env_df_fam_pca2 %>%
  nest(data = -family) %>% filter(family == "Pinaceae") %>%
  mutate(hypervol = map(data, function(.x) {
    hyper <- hypervolume_gaussian(.x[3:5])
    hyper
  }))
#save(env_df_pin_hypervolume, file = "data/env_df_pin_hypervolume.RData")
#load(file = "data/env_df_pin_hypervolume.RData")

env_df_hypervolume <- bind_rows(env_df_sci_hypervolume,
                                env_df_gne_hypervolume,
                                env_df_wel_hypervolume,
                                env_df_cyc_hypervolume,
                                env_df_ara_hypervolume,
                                env_df_zam_hypervolume,
                                env_df_eph_hypervolume,
                                env_df_gin_hypervolume,
                                env_df_pod_hypervolume,
                                env_df_tax_hypervolume,
                                env_df_cup_hypervolume,
                                env_df_pin_hypervolume)

#save(env_df_hypervolume, file = "env_df_hypervolume.RData")
#load(file = "env_df_hypervolume.RData")

# persistence homology (PH) calculations
subsample.distance_demo <- function (x, size, d, d.max = NULL, 
                                     replacement = FALSE, latlong = FALSE,
                                     echo = FALSE) {
  
  if (missing(x))
    stop("Must define a spatial object")
  
  if (missing(d))
    stop("Must define minimum separation distance")
  
  if (!is.null(d.max)) {
    if (d.max <= d) stop("Maximum distance must be larger than minimum")
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
      } else {
        deval <- any(pts.dist < d, na.rm = TRUE) | any(pts.dist > d.max, 
                                                       na.rm = TRUE)
      }
      
      nsamp = nsamp + 1
      
      if (echo)
        cat("Sample iteration=", nsamp, "\n")
      
      if (nsamp == nrow(x))
        break
    }
    
    if (echo) {
      cat("\n", "Min distance for", i, "=", min(pts.dist, na.rm = TRUE), "\n")
      cat(" Max distance for", i, "=", max(pts.dist, na.rm = TRUE), "\n")
    }
    
    if (nsamp == nrow(x)) {
      message(paste0("Warning: sampling cannot converge at n=", size, 
                     " returning n=", nrow(s)))
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
    
    for(j in 1:ncol(df)) {
      
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

# automating the calculations for PH 
env_df_hypervolume2 <- env_df_hypervolume %>%
  mutate(randomHyper = purrr::map(hypervol, function(.x) {
    out <- apply(.x@RandomPoints, 2, base::scale)
  })) %>%
  mutate(persistenceHomology = purrr::map(randomHyper, function(.x) {
    set.seed(145) # setting seed for reproducibility
    subset_hyper <- as.data.frame(.x) 
    subset_hyper <- subsample_distance(subset_hyper, size = 300, 
                                       d = 0.2) # subseting  
    PH_calc <- calculate_homology(subset_hyper, dim = 2) # calculating PH
    return(list(hyper_subset = subset_hyper,
                persHomol = as.data.frame(PH_calc)))
  }))

env_df_hypervolume3 <- env_df_hypervolume2 %>%
  mutate(hull = purrr::map(randomHyper, function(.x) {
    hullout <- with(as.data.frame(.x), chull(PC1, PC2))
    hullout
  }),
  out = map2(randomHyper, hull, ~ .x[.y,,drop = FALSE]))

index <- which(env_df_hypervolume3$family %in% tr_fam$tip.label)
env_df_hypervolume3 <- env_df_hypervolume3[index,]

#save(env_df_hypervolume3, file = "data/env_df_hypervolume3.RData")

# plots of the point clouds

plot_points <- list()
for (i in 1:nrow(env_df_hypervolume3)) {
  plot_points[[i]] <- ggplot(as.data.frame(env_df_hypervolume3$out[[i]]), 
                             aes(PC1, PC2, col = PC3)) +
    geom_polygon(alpha = 0.05, fill = "royalblue2", col = NA) +
    geom_point(data = as.data.frame(env_df_hypervolume3$randomHyper[[i]]),
               aes(PC1, PC2, col = PC3), inherit.aes = FALSE, alpha = 0.07, 
               size = 0.2) +
    geom_point(data = env_df_hypervolume3$persistenceHomology[[i]]$hyper_subset,
               aes(PC1, PC2, col = PC3), 
               inherit.aes = FALSE, alpha = 1, size = 1.5) +
    ggtitle(paste0(env_df_hypervolume3$family[i])) +
    scale_color_viridis_c() +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 9),
          plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
  
}

# plots of the persistence diagrams
plot_persistence <- list()
for (i in 1:nrow(env_df_hypervolume3)) {
  plot_persistence[[i]] <- ggplot(data = env_df_hypervolume3$persistenceHomology[[i]]$persHomol,
                                  aes(x = birth, 
                                      y = death, 
                                      shape = as.factor(dimension), 
                                      color = as.factor(dimension))) +
    geom_point(alpha = 0.9) +
    geom_abline(slope = 1, size = 0.2) +
    xlab("Birth") +
    ylab("Death") +
    theme_linedraw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 9),
          axis.title = element_text(size = 11),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 9)) +
    scale_color_manual("Dimension", 
                       values = c("darkseagreen3", "royalblue2", "hotpink")) +
    scale_shape_manual("Dimension", values = c(15, 17, 19))
}

final_plot <- c(rbind(plot_points, plot_persistence))

pdf("figures/Fig1_F.pdf", height = 9, width = 10)
wrap_plots(final_plot) + plot_layout(ncol = 4)
dev.off()

# pairwise bottleneck distances of the persistence homology 
# estimating pairwise bottleneck distances from the PHs of the species

# function to estimate the distance between two persistence diagrams
homology_distance <- function(x, y, dim = 2, ...) {
  output <- c()
  for (i in 0:dim) {
    dist <- TDA::bottleneck(as.matrix(x), as.matrix(y), i)
    output <- cbind(output, dist)
    colnames(output)[i + 1] <- paste(i)
  }
  return(as.data.frame(output))
}

# pairwise distance of the persistence (bottleneck distances)
pairwise_persistence <- function(df) {
  output <- c()
  # for loop for pairwise interactions
  for (i in 1:nrow(df)) {
    for (j in 1:nrow(df)) {
      if (i <= j) { # excluding self-comparison and double-comparison
        next
      } else {
        pairwise_dist_out <- homology_distance(df$persistenceHomology[[i]]$persHomol,
                                               df$persistenceHomology[[j]]$persHomol)
        fam1 <- df$family[j]
        fam2 <- df$family[i]
        
        outdata <- data.frame(family1 = fam1, family2 = fam2, dist = pairwise_dist_out)
        output <- as.data.frame(bind_rows(output, outdata))
      }
    }
  }
  
  output_final <- data.frame(family1 = output$family2, 
                             family2 = output$family1) %>%
    bind_cols(., output[3:5]) %>%
    bind_rows(., output)
  
  return(as.data.frame(output_final))
}

homology_distance(env_df_hypervolume3$persistenceHomology[[1]]$persHomol,
                  env_df_hypervolume3$persistenceHomology[[2]]$persHomol)

pairw_persist_hypervol <- pairwise_persistence(env_df_hypervolume3)
colnames(pairw_persist_hypervol) <- c("family1", "family2", "dist.0", "dist.1", "dist.2")

plot_dist0 <- ggplot(pairw_persist_hypervol, 
                     aes(family1, family2, fill = dist.0)) + 
  geom_tile() +
  scale_fill_viridis_c("Distance\n(Dim 0)", option = "mako") + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = -0.0)) +
  xlab("") +  scale_x_discrete(drop = FALSE) +
  ylab("Family 2")


plot_dist1 <- ggplot(pairw_persist_hypervol, 
                     aes(family1, family2, fill = dist.1)) + 
  geom_tile() +
  scale_fill_viridis_c("Distance\n(Dim 1)", option = "mako") + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = -0.0),
        axis.text.y = element_blank()) +
  xlab("Family 1") + 
  ylab("")


plot_dist2 <- ggplot(pairw_persist_hypervol, 
                     aes(family1, family2, fill = dist.2)) + 
  geom_tile() +
  scale_fill_viridis_c("Distance\n(Dim 2)", option = "mako") + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = -0.0),
        axis.text.y = element_blank()) +
  xlab("") + 
  ylab("")

final_plot_distance <- wrap_plots(plot_dist0, plot_dist1, plot_dist2)

pdf("figures/Fig2_F.pdf", width = 12, height = 4)
final_plot_distance
dev.off()

# summary dimensions
summary_persistence <- function(df) {
  output <- c()
  # for loop for pairwise interactions
  for (i in 1:nrow(df)) {
    outdata <- df$persistenceHomology[[i]]$persHomol %>%
      mutate(surv = death - birth) %>%
      group_by(dimension) %>%
      summarise(avg_surv = mean(surv),
                sd_surv = sd(surv),
                max_surv = max(surv)) %>%
      mutate(family = df$family[i])
    output <- as.data.frame(bind_rows(output, outdata))
  }
  output <- output[c(5, 1, 2, 3, 4)]
  return(as.data.frame(output))
}

# calculating average survival of connected components
persistence_info_dt <- summary_persistence(env_df_hypervolume3) 

volume_hyper <-  env_df_hypervolume3 %>%
  mutate(hv_vol = purrr::map(hypervol, function(.x) {
    .x@Volume
  })) %>%
  dplyr::select(family, hv_vol) %>%
  tidyr::unnest(hv_vol) %>%
  right_join(., persistence_info_dt) 

# subsetting only species in the tree
index_vol <- which(volume_hyper$family %in% tr_fam$tip.label)
volume_hyper <- volume_hyper[index_vol,]

outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

tr_fam_cleaned <- drop.tip(tr_fam, as.vector(outersect(tr_fam$tip.label, 
                                                              volume_hyper$family)))
plot(tr_fam_cleaned)

# plotting the hypervolume for the main text
  # plotting phylogeny with continuous traits

layout(matrix(1:4, ncol = 2, byrow = T))

volume_hyper_dim0 <- volume_hyper %>%
  filter(dimension == 0) %>%
  droplevels() %>%
  tibble::column_to_rownames(., "family")

# volume of hypervolume
volume_hyper_vol_plot <- setNames(volume_hyper_dim0[, 1], 
                                  rownames(volume_hyper_dim0))
obj1 <- contMap(tr_fam_cleaned, volume_hyper_vol_plot, plot = FALSE)
obj1 <- setMap(obj1, invert = TRUE)
plot(obj1, fsize = c(0.9, 0.5), outline = TRUE, lwd = c(4, 7), 
     leg.txt = "Volume")

# dimension
volume_hyper_dim0_plot <- setNames(volume_hyper_dim0[, 3], 
                                   rownames(volume_hyper_dim0))
obj2 <- contMap(tr_fam_cleaned, volume_hyper_dim0_plot, plot = FALSE)
obj2 <- setMap(obj2, invert = TRUE)
plot(obj2, fsize = c(0.9, 0.5), outline = TRUE, lwd = c(4, 7), 
     leg.txt = "Dim 0")

# dimension
volume_hyper_dim1 <- volume_hyper %>%
  filter(dimension == 1) %>%
  droplevels() %>%
  tibble::column_to_rownames(., "family")

volume_hyper_dim1_plot <- setNames(volume_hyper_dim1[, 3], 
                                   rownames(volume_hyper_dim1))
obj2 <- contMap(tr_fam_cleaned, volume_hyper_dim1_plot, plot = FALSE)
obj2 <- setMap(obj2, invert = TRUE)
plot(obj2, fsize = c(0.9, 0.5), outline = TRUE, lwd = c(4, 7), 
     leg.txt = "Dim 1")

# dimension
volume_hyper_dim2 <- volume_hyper %>%
  filter(dimension == 2) %>%
  droplevels() %>%
  tibble::column_to_rownames(., "family")

volume_hyper_dim2 <- setNames(volume_hyper_dim2[, 3], 
                              rownames(volume_hyper_dim2))
obj2 <- contMap(tr_fam_cleaned, volume_hyper_dim2, plot = FALSE)
obj2 <- setMap(obj2, invert = TRUE)
plot(obj2, fsize = c(0.9, 0.5), outline = TRUE, lwd = c(4, 7), 
     leg.txt = "Dim 2")

#
#
#fancydt <- volume_hyper %>% dplyr::select(-sd_surv) %>%
#  pivot_wider(., values_from = avg_surv, names_from = dimension) %>%
#  tibble::column_to_rownames(., "species")
#
#fancyPlot<-fancyTree(tr_fam_cleaned,type="scattergram",X=fancydt[,c(1:4)])
#
#

calculate_physig <- function(tree, 
                             data,
                             names = NULL,
                             method = "K",
                             test = TRUE, 
                             nsim = 999, 
                             niter = 999,
                             dim = 2, ...){
  output <- c()
  inner <- c()
  for (j in 0:dim) {
    volume_hyper <- data %>%
      filter(dimension == j) %>%
      droplevels() %>%
      tibble::column_to_rownames(., "family") # subsetting the data and dim
    
    volume_hyper <- setNames(volume_hyper[,1], rownames(volume_hyper))
    signalphy <- phylosig(tree, 
                          volume_hyper, method = method,
                          test = test, nsim = nsim, niter = niter) # physignal calculation
    if (method == "K") {
      inner <- bind_rows(inner, data.frame(dim = j,
                                           estimate = signalphy$K,
                                           sim = signalphy$sim.K
                                           pval = signalphy$P))
    } else {
      inner <- bind_rows(inner, data.frame(dim = j,
                                           estimate = signalphy$lambda,
                                           pval = signalphy$P))
    }

  }
  output <- bind_rows(output, inner)
  return(output)
}

hypervolume_sig_K <- calculate_physig(tr_fam_cleaned, 
                                      volume_hyper[c(1, 2, 3)]) %>%
  mutate(type = "Hypervolume")

avgdist_sig_K <- calculate_physig(tr_fam_cleaned, 
                                  volume_hyper[c(1, 4, 3)]) %>%
  mutate(type = "Mean dist")

sddist_sig_K <- calculate_physig(tr_fam_cleaned, 
                                 volume_hyper[c(1, 5, 3)]) %>%
  mutate(type = "SD dist")

maxdist_sig_K <- calculate_physig(tr_fam_cleaned, 
                                  volume_hyper[c(1, 6, 3)]) %>%
  mutate(type = "Max dist")

hypervolume_sig_L <- calculate_physig(tr_fam_cleaned, 
                                      volume_hyper[c(1, 2, 3)],
                                      method = "lambda") %>%
  mutate(type = "Hypervolume")

avgdist_sig_L <- calculate_physig(tr_fam_cleaned, 
                                  volume_hyper[c(1, 4, 3)],
                                  method = "lambda") %>%
  mutate(type = "Mean dist")

sddist_sig_L <- calculate_physig(tr_fam_cleaned, 
                                 volume_hyper[c(1, 5, 3)],
                                 method = "lambda") %>%
  mutate(type = "SD dist")

maxdist_sig_L <- calculate_physig(tr_fam_cleaned, 
                                  volume_hyper[c(1, 6, 3)],
                                  method = "lambda") %>%
  mutate(type = "Max dist")

final_signal <- bind_rows(cbind(method = "K", hypervolume_sig_K),
                          cbind(method = "L", hypervolume_sig_L),
                          cbind(method = "K", avgdist_sig_K),
                          cbind(method = "L", avgdist_sig_L),
                          cbind(method = "K", sddist_sig_K),
                          cbind(method = "L", sddist_sig_L),
                          cbind(method = "K", maxdist_sig_K),
                          cbind(method = "L", maxdist_sig_L))

#
#
#dummy <- final_signal  %>%
#  group_by(dim, method, type) %>%
#  summarize(mean = mean(estimate),
#            pvalue = mean(pval))
#
#

signal_plot <- final_signal %>%
  filter(method == "K") %>%
  ggplot(aes(x = estimate, col = type, fill = type)) +
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
