#### Araucaria Final analysis (Ecology Letters)
### Before revisions (i.e. version prior to rejection)
### Original script folder and name: Niche_hole>Alignment_Tree_Building_Jan2023_Fixed1.R
## June 2023
# Note that custom-made functions are given as they are used.




### 1) Constructing the phylogeny using 2 barcoding data 

## Functions
nwords <- function(string, pseudo=F){
  ifelse( pseudo, 
          pattern <- "\\S+", 
          pattern <- "[[:alpha:]]+" 
  )
  str_count(string, pattern)
} 



duplicateToConsensus <- function(sequences){
  ##getting the names in the sequence
  metanames <- stringr::word(as.character(names(sequences)), 
                             start = 2, end = 3)
  original_length <- length(sequences) #original length of the sequences
  # loop to find the duplicates, create the consensus, replace in the original sequence file and fix their names
  for(i in 1:original_length){
    if(sum(str_detect(metanames, metanames[i])) == 1){
      names(sequences)[i] <- metanames[i]
      metanames[i] <- "DONE"
    } else {
      if(metanames[i] == "DONE"){
        next
      } else {
        consensus <- msaConsensusSequence(msa(sequences[which(str_detect(metanames, metanames[i])),]))
        sequences[length(sequences) + 1,] <- consensus
        names(sequences)[length(sequences)] <- paste(metanames[i])
        metanames[which(str_detect(metanames, metanames[i]))] <- "DONE"
      }
    }
  }
  sequences[which(nwords(names(sequences)) == 2),]
}


### Multiple alignment tutorial ##
# https://bioconductor.org/packages/devel/bioc/vignettes/msa/inst/doc/msa.pdf

#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("msa")
#BiocManager::install("phyloseq")
#BiocManager::install("ggtree")

### Loading packages

library(msa)
library(stringr)
library(stringi)
library(alakazam)
library(seqinr)
library(ape)
library(phyloseq)
library(ggplot2)
library(scales)
library(phangorn)
library(adegenet)
library(phytools)



## Loading 
mySequences <- readAAStringSet("Araucaria_matKrbcL_combined.txt")
cleaned_mySequences <- duplicateToConsensus(mySequences)
myFirstAlignment <- msa(cleaned_mySequences)

# detecting which ones are of the same species
conMat <- consensusMatrix(myFirstAlignment)

#Converting to distance matrix of alignment
MyConvertedAlign <- msaConvert(myFirstAlignment, 
                               type="seqinr::alignment")
dist_matrix_align <- seqinr::dist.alignment(MyConvertedAlign, "identity")

## Plotting the tree
araucariaTree_rooted <- nj(dist_matrix_align)
araucariaTree_rooted$root.edge <- "Prumnopitys ferruginoides"
araucariaTree_rooted  <- ape::root(araucariaTree_rooted, outgroup = "Prumnopitys ferruginoides")
plot(araucariaTree_rooted, main="Araucariae rbcL tree")
araucariaTree_rooted$tip.label

dna2 <- as.phyDat(myFirstAlignment)
fun <- function(x) upgma(dist.ml(x))
bs_upgma <- bootstrap.phyDat(dna2,  fun)
plotBS(araucariaTree_rooted, bs_upgma, main="UPGMA")

#Parsimony
parsimony(araucariaTree_rooted, dna2)

treeRatchet  <- pratchet(dna2, trace = 0, minit=200)
parsimony(treeRatchet, dna2)

treeRatchet  <- acctran(treeRatchet, dna2)
treeRatchet  <- multi2di(treeRatchet) ## must use multi2di instead of di2 function to resolve polytomies

if(inherits(treeRatchet, "multiPhylo")){
  treeRatchet <- unique(treeRatchet)
}

treeRatchet$edge.length <- treeRatchet$edge.length + 4000
plotBS(midpoint(treeRatchet), type="phylogram",
       digit = 2)

treeRatchet$root.edge <- "Prumnopitys ferruginoides"
treeRatchet <- root(treeRatchet, outgroup =  "Prumnopitys ferruginoides")





### 2) Obtaining the environmental variables given the coordinates of each of the species.
###.   Note that this can be quite tedious. You can skip this step by loading the "saved" data.


### ** ~~~~ To replicate Ecol Letters analysis, can Skip from here ~~~~ ** ###
### loading the environmental data
path_envs <- list.files(path = '/Users/s23jm9/Dropbox/University of Aberdeen/Research-Projects/1.Research Projects/Triatomidae-project/wc2-5', 
                        pattern='.bil$', recursive=TRUE, full=TRUE)
envsdt <- stack(path_envs) ## I have these saved and can share. You can also get these from worldclim and other precipitation elevation databases

library(sp)
library(raster)
library(dplyr)
library(CoordinateCleaner)

levels(as.factor(araucaria_data$species))
### At the moment need to download and then save as .txt file for import (manually)
araucaria_data <- readr::read_tsv("/Users/s23jm9/Dropbox/University of Aberdeen/Research-Projects/1.Research Projects/Niche_holes/0181701-210914110416597.txt")

araucaria_clean <- araucaria_data  %>%
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
  distinct(decimallongitude, decimallatitude,specieskey,datasetkey, .keep_all = TRUE) %>%
  glimpse() 

levels(as.factor(araucaria_clean$species))

### creating a SP object to subset env data
araucaria_clean_sp <- araucaria_clean
coordinates(araucaria_clean_sp) <- ~decimallongitude+decimallatitude
crs(araucaria_clean_sp) <- crs(envsdt)

# extracting env data
env_araucaria <- raster::extract(envsdt, araucaria_clean_sp, 
                                 fun = mean, 
                                 na.rm = TRUE, sp = TRUE)
env_araucaria_df <- data.frame(env_araucaria)

## cleaning the data from extra stuff
env_araucaria_df_clean <- env_araucaria_df %>% dplyr::select(species, 
                                                             decimallongitude, 
                                                             decimallatitude,
                                                             bio1:tmin9)
## maintaing in the dataset only those species also present in the ncbi data
ncbipresent_species <- env_araucaria_df_clean$species %in% intersect(env_araucaria_df_clean$species, 
                                                                     names(cleaned_mySequences)) 
env_araucaria_df_clean2 <- env_araucaria_df_clean[ncbipresent_species,] %>%
  na.omit() %>%
  droplevels()

#write.csv(env_araucaria_df_clean2, "Environmental_Data.csv")

env_araucaria_df_clean2 <- read.csv("/Users/s23jm9/Dropbox/University of Aberdeen/Research-Projects/1.Research Projects/Niche_holes/Environmental_Data.csv",
                                    stringsAsFactors = TRUE,
                                    header = TRUE) %>%
  dplyr::select(-X) %>%
  na.omit() %>%
  droplevels()




#### ~~~~ Calculating Hypervolumes and holes therein using persistence homology ~~~~ ####
#### ~~~~~   Hypervolumes 
library(hypervolume)
library(tidyr)
library(purrr)
library(TDAstats)

head(env_araucaria_df_clean_pca2)
colnames(env_araucaria_df_clean_pca2) <- c("species", "lon", "lat", "PC1", "PC2", "PC3")


## Loading the data saved already (to speed up processing time)


env_araucaria_df_clean_hypervolume
## Calculating hypervolumes ##
env_araucaria_df_clean_hypervolume <- env_araucaria_df_clean_pca2 %>%
  nest(-species) %>%
  mutate(hypervol = map(data, function(.x){
    hyper <- hypervolume_gaussian(.x[3:5])
    hyper
  }))

env_araucaria_df_clean_hypervolume2 <-env_araucaria_df_clean_hypervolume %>%
  mutate(randomHyper = map(hypervol, function(.x){
    out <- apply(.x@RandomPoints, 2, base::scale)
  }))









#### ~~ Persistence Homology calculations ~~~ ####
## Functions 
subsample.distance_demo <- function (x, size, d, d.max = NULL, replacement = FALSE, latlong = FALSE, 
                                     echo = FALSE) 
{
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




subsample_distance <- function(df, size, d, ...){
  sub.meuse <- c()
  for(i in 1:ncol(df)){
    for(j in 1:ncol(df)){
      if(i == j | j > i){
        next
      } else {
        testsp <- df
        sp::coordinates(testsp) <- c(i, j)
        sub <- spatialEco::subsample.distance(testsp, size = size/ncol(df), d = d) ##d is for the minimum distance between the points sampled
        sub.meuse <- bind_rows(sub.meuse, as.data.frame(sub))
      }
    }
  }
  return(sub.meuse)
}
####
###### ~~~~~ automating the calculations for PH ~~~~~~~ #####
env_araucaria_df_clean_hypervolume2 <- env_araucaria_df_clean_hypervolume %>%
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


env_araucaria_df_clean_hypervolume3 <- env_araucaria_df_clean_hypervolume2 %>%
  mutate(hull = map(randomHyper, function(.x){
    hullout <- with(as.data.frame(.x), chull(PC1, PC2))
    hullout
  }),
  out = map2(randomHyper, hull, ~ .x[.y,,drop=FALSE]))

### ** ~~~~ to here ~~~~ ** ###




#* You can load the RData to avoid having to calculate homology and hypevol
env_araucaria_df_clean_hypervolume3 <- readRDS("/Users/s23jm9/Dropbox/University of Aberdeen/Research-Projects/1.Research Projects/Niche_holes/env_araucaria_df_clean_hypervolume3.RDS")
env_araucaria_df_clean_hypervolume3


index <- which(env_araucaria_df_clean_hypervolume3$species %in% araucariaTree_rooted$tip.label)
env_araucaria_df_clean_hypervolume3 <- env_araucaria_df_clean_hypervolume3[index,]


## Plots of the point clouds
plot_points <- list()
for(i in 1:nrow(env_araucaria_df_clean_hypervolume3)){
  plot_points[[i]] <- ggplot(as.data.frame(env_araucaria_df_clean_hypervolume3$out[[i]]), 
                             aes(PC1, PC2, col = PC3)) +
    geom_polygon(alpha = 0.05, fill = "royalblue2", col = NA) +
    geom_point(data = as.data.frame(env_araucaria_df_clean_hypervolume3$randomHyper[[i]]),
               aes(PC1, PC2, col = PC3), inherit.aes = FALSE, alpha = 0.07, size = 0.2) +
    geom_point(data = env_araucaria_df_clean_hypervolume3$persistenceHomology[[i]]$hyper_subset,
               aes(PC1, PC2, col = PC3), 
               inherit.aes = FALSE, alpha = 1, size = 1.5) +
    ggtitle(paste0(env_araucaria_df_clean_hypervolume3$species[i])) +
    scale_color_viridis_c() +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 9),
          plot.title = element_text(size = 9, face = "bold.italic", hjust = 0.5))
  
}


## plots of the persistence diagrams ##
plot_persistence <- list()
for(i in 1:nrow(env_araucaria_df_clean_hypervolume3)){
  plot_persistence[[i]] <- ggplot(data = env_araucaria_df_clean_hypervolume3$persistenceHomology[[i]]$persHomol,
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




