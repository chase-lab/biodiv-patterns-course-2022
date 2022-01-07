library(sp)
library(tidyverse)
library(cowplot)

##---------------------------------------------------------------------------------------------------------------------------------------------------
# Process-based model I: geometric population growth
# simulate population growth and plot for three cases:
# 1) b > d
# 2) b == d
# 3) b < d

# initial population size
N0 = 100
# per-capita birth rate
b = 1.2
# per-capita mortality (e.g., annual plant)
d = 1
# initialise vector to store population size each year, and assign initial population size as first value
N <- NULL
N[1] <- N0

# simulate population growth for 10 year
for(t in 1:9) {
  N[t+1] = N[t] + (b-d) * N[t]
}

# b > d results
N_growing <- N
rm(N)

# b == d results
b = 1; d = 1
N <- NULL
N[1] <- N0

# simulate population growth for 10 year
for(t in 1:9) {
  N[t+1] = N[t] + (b-d) * N[t]
}
N_static <-  N
rm(N)

# b < d results
b = 0.9; d = 1
N <- NULL
N[1] <- N0

# simulate population growth for 10 year
for(t in 1:9) {
  N[t+1] = N[t] + (b-d) * N[t]
}
N_shrinking <-  N
rm(N)

# plot results
plot(1:length(N_growing), N_growing, type = 'l',
     ylim = c(0, 600),
     xlab = 'Time [years]',
     ylab = 'Population size or density')
lines(1:length(N_static), N_static, type = 'l')
lines(1:length(N_shrinking), N_shrinking, type = 'l')

##---------------------------------------------------------------------------------------------------------------------------------------------------
# Stochasticity: extending models to include structured randomness.
# - Dennis et al 1995 Ecological Monographs; Hilborn & Mangel 1997 The Ecological Detective

# Demographic and environmental stochasticity are most commonly discussed forms of stochasticity in Ecology
# other types: catastrophes (e.g., Lande 1993 Am Nat); measurement or observation error (e.g., Denis et al 1991 Ecol Monographs;
# Knape & de Valpine 2011 Ecol Lett);  process error (e.g., Denis et al 1991 Ecol Monographs; Thibaut & Connolly 2020 Ecology)

# Demographic stochasticity describes realised variability in demographic processes such as births, deaths and migration due
# to their probabilistic nature (see also Clark 2009 TREE for argument that stochastic forces exist only in models; e.g., "...processes perceived
# as stochastic at one level of abstraction have explanations at another.").

# Environmental stochasticity describes variablity in extrinsic environmental conditions (e.g., due to variation in temperature,
# precipitation, wind, etc)

# Extend geometric population model to include demographic stochasticity

# We have to decide how to change our model to include stochasticity. For example, we could assume that the parameter values of both the birth
# and death processes in the geometric growth model as random draws from probability distributions.
# Here, we'll instead use the approach taken in the metacommunity model
# that we will study later today (see also Adler & Drake 2008 Am Nat, Shoemaker et al. 2020 Ecology).

# We assume N[t+1] ~ Pois(lambda = N[t] + (b-d)*N[t]),
# N[t+1], the population size next year, is a random draw from a
# Poisson distribution with lambda (i.e., the one and only parameter of the Poisson distribution)
# equal to the expected population size under our model of geometric growth.

# define per-capita birth rate
b = 1.2
# and per-capita mortality
d = 1

# initialise vector to store population size
N <- NULL
N[1] <- N0

# simulate population growth for 10 years
for(t in 1:9) {

  N[t+1] = rpois(n = 1,
                 lambda = N[t] + (b-d) * N[t])
}

plot(1:length(N_growing), N_growing, type = 'l',
     ylim = c(0, 600),
     xlab = 'Time [years]',
     ylab = 'Population size or density')
lines(1:length(N), N,
      col = 'grey')

# do multiple simulations to visualise variability more clearly
nsims = 20
nyears = 9
N_matrix = matrix(data = NA,
                  nrow = nsims,
                  ncol = nyears + 1,
                  byrow = TRUE)
# initial population
N_matrix[,1] = N0

for(sim in 1:nsims){
  # print counter to screen
  print(paste('simulation ', sim, 'of ', nsims))
  for(t in 1:nyears){
    N_matrix[sim, t+1] = rpois(n = 1,
                   lambda = N_matrix[sim, t] + (b-d) * N_matrix[sim, t])
  }
}

# plot
matplot(t(N_matrix), type = 'l')
lines(1:length(N_growing), N_growing, lwd = 3)

## Environmental stochasiticity

# To include environmental stochasticity in our model of geometric population growth, we include an extra
# term, specifically,
# N[t+1] = N[t] + (b - d)*N[t] + sigma[t] * N[t],


# sigma[t] is a time series of environmental variation (e.g., temperature), and follows
# sigma[t+1] = a*sigma[t] + c*phi[t],
# where c scales the magnitude of the environmental variation phi[t] ~ N(0,1).

# set c = (1 - a^2)^0.5, which means var(sigma) is equal for all values of a.
# a controls autocorrelation:
# 0 < a < 1 environmental fluctuations are positively correlated,
# -1 < a < 0 environmental fluctuations are negatively correlated,
# a = 0, environmental fluctuations are uncorrelated.

phi = rnorm(nyears, mean = 0, sd = 1)
a = 0 # can to control autocorrelation of environmental fluctuations
c = (1 - a^2)^0.5
# initial vector to store environmental fluctuations
sigma <- NULL
sigma[1] = 0
# loop to create environmental fluctuations
for(t in 1:nyears){
  sigma[t+1] = a * sigma[t] + c*phi[t]
}

# visual inspection of environmental variation
plot(1:(nyears+1), sigma,
     xlab = 'Time [years]',
     ylab = 'Environmental condition\n(deviation from mean)',
     type = 'b')
abline(c(0,0), lty = 2)

# create plot with positive, negative and no autocorrelation
uncorrelated_sigma = sigma;
rm(sigma)

# positive autocorrelation
a = 1
sigma <- NULL
sigma[1] = 0
# loop to create environmental fluctuations
for(t in 1:nyears){
  sigma[t+1] = a * sigma[t] + c*phi[t]
}
pos_corr_sigma = sigma
rm(sigma)

# negative autocorrelation
a = -0.8
sigma <- NULL
sigma[1] = 0
# loop to create environmental fluctuations
for(t in 1:nyears){
  sigma[t+1] = a * sigma[t] + c*phi[t]
}
neg_corr_sigma = sigma
rm(sigma)

# plot
plot(1:(nyears+1), uncorrelated_sigma,
     xlab = 'Time [years]',
     ylab = 'Environmental condition\n(deviation from mean)',
     type = 'b',
     ylim = c(min(c(uncorrelated_sigma, pos_corr_sigma, neg_corr_sigma)) - 0.1,
              max(c(uncorrelated_sigma, pos_corr_sigma, neg_corr_sigma)) + 0.1))
lines(1:(nyears+1), pos_corr_sigma, col = 2, lty = 2)
lines(1:(nyears+1), neg_corr_sigma, col = 3, lty = 3)
abline(c(0,0), lty = 2)


# Exercises:
# 1) Extend model to include demographic and environmental stochasticity (see Shoemaker et al. 2020 Ecology).

# 2) How would you measure the impact of including stochasticity in a model of geometric growth using simulations?
# Hint: it can be done using linear (statistical) models; see e.g., Blowes & Connolly 2012 Ecological Applications.
# Use simulations and linear models to quantify the impact of including different combinations of demographic and environmental stochasticity.

# 3) Reformulate geometric population growth model to include stochastic variation in births and deaths as separate processes.


##---------------------------------------------------------------------------------------------------------------------------------------------------
## Metacommunity dynamics

# But first, metapopulation dynamics

# Simulate landscape of discrete patches
# Number of sites/patches (P)
P <- 20

# Random distribution of sites
x <- runif(P,0,10)
y <- runif(P,0,10)
plot(x,y,cex=3)

# Other spatial structures
xy <- spsample(SpatialPoints(coords0), P, "regular") # regular, clustered, random
coords <- coordinates(xy)
plot(coords,cex=3)
text(coords[,1],coords[,2],1:P)

# Calculate distances from each patch to every other patch
dists <- as.matrix(dist(coords))

# Population growth rate (density independent)
r <- 0.02

# define the per-capita probabilty of dispersal (i.e., probability of an individual emigrating)
a <- 0.02

# Assume that successfully arriving in another patch (having emigrated) declines exponentially with increasind distance
# Define the rate of exponential distance decay
d <- 0.1
# visual inspection of exponential distance decay in dispersal
plot(x = 1:M,
     y = exp(-d*(1:M)),
     ylim=c(0,1),
     col = 1,
     type = 'b')

lines(x = 1:M,
      y = exp(-(d*2)*(1:M)), # change the d value to see how function shape changes
      ylim=c(0,1),
      col = 2,
      type = 'b')



# calculate dispersal from and to each patch
disp_mat <- exp(-d * dists)
diag(disp_mat) <- 0 # all emigrants are assumed to have left their home patch

# convert to probability of immigrating to each patch from every other patch
disp_mat <- apply(disp_mat, 1, function(x) x / sum(x))

# check that we have a probability (should sum to one)
colSums(disp_mat)

# visualise
landscape <- ggplot(as_tibble(coords) %>%
         rename(x = x1,
                y = x2) %>%
         mutate(patch = 1:nrow(coords))) +
  geom_label(aes(x = x, y = y, label = patch)) +
  coord_fixed() +
  theme_bw() +
  theme()

connectivity <- as_tibble(disp_mat) %>%
  mutate(to.patch = rownames(disp_mat)) %>%
  gather(key = from.patch, value = dispersal, -to.patch) %>%
  mutate(from.patch = as.numeric(as.character(from.patch)),
                to.patch = as.numeric(as.character(to.patch))) %>%
  ggplot(aes(x = from.patch, y = to.patch, fill = dispersal)) +
  geom_tile() +
  labs(x = 'From patch',
       y = 'To patch') +
  scale_fill_viridis_c() +
  coord_fixed() +
  theme_bw()

plot_grid(landscape,
          connectivity)


# Simulate dynamics
#==================

# Distribute it on the landscape: initialise model
# 100 individuals in one site

# Set initial abundance
N <- rep(0,P)
N[sample(1:P,1)] <- 100
plot(coords,cex=3)
points(coords,cex=3,pch=16,col=grey(1-N/max(N)))


# Run model one time step (individuals on the move)

# Calculate number of new individuals (r * N)
rN <- round(r * N)

# update population size
N <- N + rN

# Calculate number of emmigrants ()
E <- round(a * N)

# stochastic alternative
E <- rbinom(length(N), size = N, prob = a)

# create immigration vector (to store number of immigrants)
I <- numeric(length = P)

# Calculate immigration
# loop through patches (for each patch, simulate dispersal)
for(p in 1:P){

  if(E[p]>0) {

    # sample the sites that receive immigrants from site i based on the dispersal kernel
    # this also adds stochasticity to the model
    dest <- sample(1:P, E[p], replace=TRUE, prob=disp_mat[p,])
    # Count the number of immigrants in each site
    Ii <- table(dest)
    # store the number of immigrants in each site
    I[sort(unique(dest))] <- I[sort(unique(dest))] + Ii
  }
}
# Calculate the abundance after dispersal
N <- N - E + I
N[N<0] <- 0

# repeat a few times to see redistribution the landscape
# plot abundance in each site (grey scale)
plot(coords,cex=3)
points(coords,cex=3,pch=16,col=grey(1-N/max(N)))


# Run the model for 100 time steps
tmax = 100
N <- matrix(nrow = tmax, ncol = P)
# initial population 100 in each patch
N[1,] <- 100

#
for(t in 2:tmax){
  print(paste(t, ' of ', tmax))
  rN <- round(r * N[t-1,])
  # what about stochastic population growth?
  # rN <- round(rnorm(1, 0, 0.2) * N[t-1,])
  N[t, ] <- N[t-1] + rN
  E <- rbinom(n = ncol(N), size = N[t, ], prob = d)
  I <- numeric(P)
  for(p in 1:P){
    if(E[p]>0) {
      dest <- sample(1:P,E[p],replace = TRUE,prob=disp_mat[p,])
      Ip <- table(dest)
      I[sort(unique(dest))] <- I[sort(unique(dest))] + Ip
    }
  }
  N[t, ] <- N[t, ] - E + I
  N[N<0] <- 0
}

matplot(N, type = 'l')


#---extension to include the niche in density-independent growth rate

# assume niche is Gaussian
# Gaussian distribution has two parameters: mean (or optima, z), and sd (or breadth)

# maximum density independent growth rate
rmax = 5
# niche optima
z = 0.5
# niche breath
breadth <- 0.01
plot(x = seq(0,1,0.01),
     y = exp(-((z - seq(0,1,0.01))/(2 * breadth))^2) * rmax,
     ylim=c(0,5),
     type = 'b',
     xlab = 'Environment',
     ylab = 'Density independent growth rate (r)')

# with the function from the Thompson model
S = 10 # number of species
env_traits(species = S,
           max_r = 5,                         # max growth rate
           min_env = 0, max_env = 1,          # range of environment
           env_niche_breadth = 0.05,       # niche breadth
           optima_spacing = 'even')           # spacing of z_i (optima)


##--------simulate metacommumity dynamics with the Thompson model
# Number of patches
P = 20
# Number of species
S = 20

# set up landscape for experiment (want the same landscape for both treatments)
meta_landscape <- landscape_generate(patches = M)

# dispersal matrix
# set rate of distance decay
d = 0.1
disp_mat <- dispersal_matrix(landscape = disp_experiment,
                                  torus = TRUE,
                                  kernel_exp = d,
                                  plot = TRUE)

# generate the time series of the environmental conditions for each patch (same for each treatment)
# TODO: create an alternate environment with stronger temporal autocorrelation
env_conditions <- env_generate(landscape = disp_experiment,
                               env1Scale = 1, # temporal autocorrelation in the environment (here the environment is temporally uncorrelated)
                               timesteps = 1000)

# density independent component of model
densInd_niche <- env_traits(species = S,
                            max_r = 5,                         # max growth rate
                            min_env = 0, max_env = 1,          # range of environment
                            env_niche_breadth = 0.25,       # niche breadth
                            optima_spacing = 'even')           # spacing of z_i (optima)

# species interaction matrix:
# TODO: create more matrices for other dynamics in the paper (mixed, competitive dominance, destabilising competition)
equal_interaction_mat <- species_int_mat(species = S,
                                         intra = 1,
                                         min_inter = 1,
                                         max_inter = 1)

stabilising_interaction_mat <- species_int_mat(species = S,
                                         intra = 1,
                                         min_inter = 0,
                                         max_inter = 0.8)

# use simulateMC() function to simulate dynamics
sim_equal_comp <- simulate_MC(patches=P, species=S,
                        landscape = meta_landscape,
                        disp_mat = disp_mat,
                        env.df = env_conditions,
                        max_r = densInd_niche$max_r,
                        env_niche_breadth = densInd_niche$env_niche_breadth,
                        env_optima = densInd_niche$optima,
                        int_mat = equal_interaction_mat,
                        initialization=100, burn_in=300, timesteps=600)

sim_stabil_comp <- simulate_MC(patches=P, species=S,
                              landscape = meta_landscape,
                              disp_mat = disp_mat,
                              env.df = env_conditions,
                              max_r = densInd_niche$max_r,
                              env_niche_breadth = densInd_niche$env_niche_breadth,
                              env_optima = densInd_niche$optima,
                              int_mat = stabilising_interaction_mat,
                              initialization=100, burn_in=300, timesteps=600)

# extract data
sim_equalC_dat <- sim_equal_comp$dynamics.df %>%
  as_tibble() %>%
  # reduce to last 100 time steps
  filter(time > 499 & time < 601)


sim_stabilC_dat <- sim_stabil_comp$dynamics.df %>%
  as_tibble() %>%
  # reduce to last 100 time steps
  filter(time > 499 & time < 601)

# visual inspection

# environmental conditions:
ggplot() +
  geom_line(data = sim_equalC_dat,
            aes(x = time, y = env, colour = env,
                group = patch)) +
  scale_colour_viridis_c()

# patch level population dynamics: equal comp
ggplot() +
  facet_wrap(~patch) +
  geom_line(data = sim_equalC_dat,
            aes(x = time, y = N, colour = optima,
                group = interaction(species, patch))) +
  scale_colour_viridis_c()

ggplot() +
  facet_wrap(~patch) +
  geom_line(data = sim_stabilC_dat,
            aes(x = time, y = N, colour = optima,
                group = interaction(species, patch))) +
  scale_colour_viridis_c()

# Exercises:
# Simulate dynamics according to the classic metacommunity paradigms (Fig 2 in paper) and plot
#
