######################################################
# Additional Estimations of Functional Trait Diversity
######################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(alphahull)

library(hypervolume)
library(BAT)

library(tidyverse)

rm(list=ls())

# Load native species Hawaiian forest and species by trait matrix

datt_trait = read.csv("C:/_Lectures/Lecture_Jan_Msc/Class_2022/datt_trait.csv")
traits     = read.csv("C:/_Lectures/Lecture_Jan_Msc/Class_2022/Trait_imputed.csv")

#_______________________________________________________________________________
# # QUESTION: Do single traits differ along rainfall gradients?

# TO DO: Single trait analysis: dry and  wet areas across islands

# First explore the data set
colnames(datt_trait)
unique(datt_trait$prec_range)

tmp = dplyr::filter(datt_trait, prec_range == "dry" | prec_range == "wet")

SLA = ggplot(data= tmp, aes(x = SLA, group = prec_range, fill= prec_range)) + geom_density(adjust=1.5, alpha=.4)
WD  = ggplot(data= tmp, aes(x = WD, group = prec_range, fill= prec_range)) + geom_density(adjust=1.5, alpha=.4)
H   = ggplot(data= tmp, aes(x = H, group = prec_range, fill= prec_range)) + geom_density(adjust=1.5, alpha=.4)
N   = ggplot(data= tmp, aes(x = N, group = prec_range, fill= prec_range)) + geom_density(adjust=1.5, alpha=.4)

# TASK: Do single traits differ across islands, i.e., age?

# What is the trait space of wet and dry areas on Hawaiian islands? 
# How much do they overlap? 

#_______________________________________________________________________________
# lets explore the traits

head(traits[2:5])

traits = data.frame(species =  traits$species, 
                    SLA     =  log10(traits$SLA),     # mm^2 / mg^1  (area/mass)
                    WD      =  log10(traits$WD + 10), # g / cm^3 (mass/vol)
                    H       =  log10(traits$H),       # m 
                    N       =  log10(traits$N),       # mg / g  (mass/mass)
                    Scientific_name = traits$Scientific_name)
row.names(traits) = traits$Scientific_name

PCA         <- prcomp(traits[2:5], scale. = TRUE)
biplot(PCA)
summary(PCA)

# Interesting website about PCA 
# https://builtin.com/data-science/step-step-explanation-principal-component-analysis

PCAvalues   <- data.frame(Scientific_name = traits$Scientific_name, PCA$x)    # Extract PC axes for plotting
PCAloadings <- data.frame(Variables = rownames(PCA$rotation), PCA$rotation)   # Extract loadings of the variables

# Plot trait space of all species

PCAvalues %>% ggplot(aes(PC1,  PC2), size = 1)+ # Plot PCA 
  stat_density_2d(geom = "polygon", contour = TRUE, aes(fill = after_stat(level)),  
                  colour = "gray", bins = 10) +
  scale_fill_distiller(palette = "BuGn", direction = 1) +
  geom_jitter(alpha= 0.6,  size = 2  , colour = "turquoise4") +     #  Display the points, i.e., species 
  geom_text(data = PCAloadings, aes(x = PC1*4.7, y = PC2*4.7, label = Variables), size = 2.3) +
  geom_segment(data = PCAloadings, size = 0.2,    # Plots the loadings, i.e., traits 
               aes(x = 0, xend = PC1*4.2, y = 0, yend = PC2*4.2),
               arrow = arrow(length = unit(0.1, "cm")),colour = "black")  +
  ylim(-5, 5.5) + 
  theme_minimal()

# TASK: 
#  1. identify acquisitive and conservative plants in the trait space of Hawaiian trees 
#  2. Which plant strategies according to the species trait combinations can you identify in the trait space?

# Let's plot trait space of dry and wet species

datt_trait = dplyr::left_join(datt_trait, PCAvalues, by =  "Scientific_name") 

dry = dplyr::filter(datt_trait, prec_range == "dry") # Separating the data into dry and wet  
wet = dplyr::filter(datt_trait, prec_range == "wet")

# Check how many species are in each precipitation range, i.e., dry and wet
length(unique(dry$species))  
length(unique(wet$species)) 
length(intersect(dry$species, wet$species))

colnames(dry)

# Summarize / aggregate trait values per species 
dry1 <- dplyr::summarize(group_by(dry, Scientific_name,  prec_range),
                         PC1=max(PC1), PC2=max(PC2), PC3=max(PC3))
wet1 <- dplyr::summarize(group_by(wet, Scientific_name,  prec_range),
                         PC1=max(PC1), PC2=max(PC2), PC3=max(PC3))


dry_p = dry1  %>% ggplot(aes(PC1,  PC2))+ # Plot Tenerife PCA 
  stat_density_2d(geom = "polygon", contour = TRUE, 
                  aes(fill = after_stat(level)), colour = "gray", bins = 3.5) +
  scale_fill_distiller(palette = "Greys", direction = 1) +
  geom_jitter(alpha=0.5,  size = 2  , colour = "gold2") +                      #  Display the points
  ylim(-5, 5.5) + 
  theme_minimal()

dry_p
ggExtra:: ggMarginal(dry_p, type = "density", fill="transparent", size = 15) # add marginal distribution of the PCs

# TASK: Plot the trait space for wet


#_______________________________________________________________________________
# Do traits combinations of Hawaiian forest species overlap across rainfall gradients?
# how similar wet and dry are with respect to their trait diversity

##  https://benjaminblonder.org/hypervolume_faq.html

overall_bandwidth <- estimate_bandwidth(PCAvalues[,2:4])

dry_hv <-hypervolume_gaussian(dry1[,3:5], name = "Dry volume",
                                   kde.bandwidth=overall_bandwidth, quantile.requested=0.95)

wet_hv <-hypervolume_gaussian(wet1[,3:5], name = "Wet volume",
                                 kde.bandwidth=overall_bandwidth, quantile.requested=0.95)

# calculating overlap statistics
HVs <- hypervolume_join (dry_hv, wet_hv)

HV_set <- hypervolume_set(dry_hv, wet_hv, num.points.max = NULL,
                          verbose = TRUE, check.memory = F, distance.factor = 1)
plot(HV_set)

hypervolume_overlap_statistics(HV_set)

#_______________________________________________________________________________
# Functional richness, evenness, dispersion (Blonder et al approach)
# The difference in species number between wet and dry affects FD estimations. 
# To solve this, we rarefied by the common minimum number of species. 

min_rar <- 25

nbperm <- 20 # Number of permutations

rar <- c()   

for (i in 1:nbperm){   # Loop of rarefaction
  # dry 
  dry_i <- dry1[sample(1:nrow(dry1), size = min_rar),
               c("PC1", "PC2", "PC3")]
  
  # Compute hypervolume with a fixed bandwidth
  dry_i <-  hypervolume_gaussian(
    dry_i, kde.bandwidth = overall_bandwidth,
    quantile.requested = 0.95, quantile.requested.type = "probability",
    verbose = FALSE)
  
  dry_i <- data.frame(perm = i,
                      status = "DRY",
                      rich = kernel.alpha(dry_i),
                      eve = kernel.evenness(dry_i),
                      div = kernel.dispersion(dry_i))
    # wet
  wet_i <- wet1[sample(1:nrow(wet1), size = min_rar),
                c("PC1", "PC2", "PC3")]
  
  # Compute hypervolume with a fixed bandwidth
  wet_i <-  hypervolume_gaussian(
    wet_i, kde.bandwidth = overall_bandwidth,
    quantile.requested = 0.95, quantile.requested.type = "probability",
    verbose = FALSE)
  
  wet_i <- data.frame(perm = i,
                      status = "WET",
                      rich = kernel.alpha(wet_i),
                      eve = kernel.evenness(wet_i),
                      div = kernel.dispersion(wet_i))
  # Bind results
  rar <- rbind(rar, dry_i, wet_i)
  
  cat(paste0(round(100*i/nbperm, 0), " %; "))
}


unique(rar$status)
colnames(rar)

# Summarizing rarefied values
# https://www.scribbr.com/statistics/confidence-interval/

sum_fr <- rar %>%       
  group_by(status) %>%
  summarise( n = n(), mean = mean(rich), sd = sd(rich) ) %>%
  mutate( se = sd/sqrt(n))  %>%
  mutate( ic = se * qt((1-0.05)/2 + .5, n-1))  %>%
  mutate( mean, metric = "Functional Richness") 

sum_fe <- rar %>%
  group_by(status) %>%
  summarise( n = n(), mean = mean(eve),sd = sd(eve) ) %>%
  mutate( se = sd/sqrt(n))  %>%
  mutate( ic = se * qt((1-0.05)/2 + .5, n-1))  %>%
  mutate( mean, metric = "Functional Evenness") 

sum_fd <- rar %>%
  group_by(status) %>%
  summarise( n = n(), mean = mean(div),sd = sd(div) ) %>%
  mutate( se = sd/sqrt(n))  %>%
  mutate( ic = se * qt((1-0.05)/2 + .5, n-1))  %>%
  mutate( mean, metric = "Functional Dispersion") 

sum_metrics <-  rbind(sum_fr, sum_fe, sum_fd)

#___________________________________ Plot functional trait diversity metrics

ggplot(sum_metrics) +
  geom_linerange(aes(factor(status), 
                     ymin= mean-ic, ymax= mean+ic, color = status),  size = 0.3,   # Error bar
                 show.legend = FALSE) +
  geom_point    (aes(x = status,  
                     y=mean, color = status),  size = 1,  # error mean point
                 show.legend = FALSE) +
  facet_wrap(~ metric, scales = "free") 
