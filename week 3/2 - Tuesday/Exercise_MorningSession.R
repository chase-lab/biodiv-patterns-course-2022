# use Thompson et al metacommunity model to introduce trait and species diversity 
# trait and species diversity calculated using Chao's numbers equivalent (or Hill number) framework

rm(list=ls())
R.version # check R version

# Load libraries
library(tidyverse)
# devtools::install_github("plthompson/mcomsimr")
library(mcomsimr)
# devtools::install_github("AnneChao/iNEXT.3D")
library(iNEXT.3D)
# devtools::install_github('MoBiodiv/mobr')
library(mobr)

library(cluster)


# number of patches (or replicates in each treatment)
M <- 20

# total number of species
S <- 50

# model assumes exponential distance decay of dispersal, controlled by a rate parameter
# low values mean distance decay is weaker (==high connectivity)
# this is one of the species traits (currently all species have the same value)
d_high <- 0.05

# set up landscape for metacommutiy
meta_landscape <- landscape_generate(patches = M)   # Generates a landscape for metacommunity simulations

# dispersal matrices: It creates a matrix with dispersal probabilities
disp_mat_high <- dispersal_matrix(landscape = meta_landscape,
                                  torus = TRUE,  # We convert these coordinates into a torus to avoid edge effects.To avoid edge effects, the bounding area can be treated as a torus so that, for example, a circle pushed over the left-hand edge will re-enter the bounding area at the right-hand edge. 
                                  kernel_exp = d_high,  # exponential rate at which dispersal decreases 
                                  plot = TRUE)

# generate the time series of the (20) environmental conditions for each patch (same for each treatment)
env_conditions <- env_generate(landscape = meta_landscape,
                               env1Scale = 1000, # scale of environmental autocorrelation between 0 and 1000
                               timesteps = 1000)

# density independent component of model: these are the remaining traits in the model
# all species have the same traits, except for the location of their niche optimum
densInd_niche <- env_traits(species = S, 
                            max_r = 5,                         # max growth rate
                            min_env = 0, max_env = 1,          # range of environment
                            env_niche_breadth = 0.25,          # niche breadth
                            optima_spacing = 'random')         # spacing of z_i (optima)



# species interaction matrix: Generates density dependent matrix of per capita competition
equal_interaction_mat <- species_int_mat(species = S,
                                         intra = 1,        # intraspecific competition coefficient
                                         min_inter = 0,
                                         max_inter = 0.5)


# Simulate metacommunity dynamics
#================================
sim_meta <- simulate_MC(patches=M, species=S, 
                        # env1Scale = 1,
                        landscape = meta_landscape,
                        disp_mat = disp_mat_high, 
                        env.df = env_conditions, 
                        max_r = densInd_niche$max_r,   # intrinsic growth rate in optimal environment
                        env_niche_breadth = densInd_niche$env_niche_breadth, # standard deviation of environmental niche breadth
                        env_optima = densInd_niche$optima,
                        int_mat = equal_interaction_mat,
                        initialization=100, burn_in=300, timesteps=600)


# Extract data
#=============
# Species data
meta_long <- sim_meta[["dynamics.df"]] %>% 
  as_tibble()

# chose a single time point for analysis
sim_dat <- meta_long %>% 
  filter(time==200)

# long to wide
sim_wide <- sim_dat %>% 
  dplyr::select(patch, species, N) %>% 
  arrange(patch, species) %>% 
  spread(key = species, value = N, fill = 0)

# tidy up species names
spp_names <- paste0('sp', colnames(sim_wide)[-1])
colnames(sim_wide)[-1] <- spp_names

# iNEXT.3D also requires the community data to be wide. But this time, species are in rows, and sites are columns!
wide_t <- t(sim_wide[-1])

# remove species that were not observed, and name columns
wide_t <- wide_t[rowSums(wide_t) > 0,] %>% 
  as.matrix()
colnames(wide_t) <- paste0('patch', 1:ncol(wide_t))

# need a trait-distance matrix to calculate trait-group diversity
traits_long <- densInd_niche %>% 
  mutate(dispersal = d_high,
         species = paste0('sp', species))

# only want to use observed species (i.e., N > 0) to calculate distance matrix
# I've used Gower distance, but other distance measures (e.g., Euclidean) are also 
# appropriate for these data as we only have numeric traits. 
# daisy: calculates pairwise dissimilarities (distances) between observations in the data set
gower_distM <- cluster::daisy(x = traits_long %>% 
                                filter(species %in% rownames(wide_t)) %>% 
                                select(-species), 
                              metric = "gower") %>% 
  as.matrix()

# need to fix names of distance matrix (or iNEXT function blows up)
rownames(gower_distM) <- traits_long %>% 
  filter(species %in% rownames(wide_t)) %>% 
  select(species) %>% pull()
colnames(gower_distM) <- rownames(gower_distM)


########################################################
# Anne Chao et al approach

# calcucate taxonomic and trait diversity
TD <- obs3D(data = wide_t, 
                diversity = 'TD', 
                q = c(0, 1, 2), 
                datatype = 'abundance', 
                nboot = 199)

FD <- obs3D(data = wide_t, 
                diversity = 'FD', 
                q = c(0, 1, 2), 
                datatype = 'abundance', 
                nboot = 199, 
                FDdistM = gower_distM, FDtype = 'tau_values', 
                FDtau = NULL # this uses the dmean pairwise distance
                ) 

# If NULL (default), then threshold is set to be the mean distance between any two 
# individuals randomly selected from the pooled assemblage (i.e., quadratic entropy).

# join to plot on same panels
bind_rows(TD %>% 
            select(Order.q, qD, Assemblage) %>% 
            rename(value = qD) %>% 
            mutate(diversity = 'Taxonomic'),
          FD %>% 
            select(Order.q, qFD, Assemblage) %>% 
            rename(value = qFD) %>% 
            mutate(diversity = 'Functional')) %>% 
  ggplot() +
  facet_wrap(~Order.q) +
  geom_boxplot(aes(x = diversity, y = value))


################################################################################
### Task: use different tau values and test how FD metric changes
# E.g. using 
min(gower_distM[gower_distM>0])
max(gower_distM)

