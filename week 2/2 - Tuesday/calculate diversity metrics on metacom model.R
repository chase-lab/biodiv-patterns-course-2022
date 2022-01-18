
# bolt on diversity metrics
# continuation of Thompson model of tuesday in week 1


sim_equalC_dat
sim_stabilC_dat


# chose a single time point for analysis
high_dat <- high_long %>%
filter(time>0)

low_dat <- low_long %>%
filter(time>0)

ext_dat <- ext_long %>%
filter(time>0)

# label treatments and combine
exp_dat <- bind_rows(sim_equalC_dat %>%
                       mutate(treatment = 'equalizing'),
                     sim_stabilC_dat %>%
                       mutate(treatment = 'stabilizing'))



exp_results_alpha_scale <- exp_dat %>%
  # remove zeroes
  filter(N > 0) %>%
  group_by(treatment, patch, time) %>%
  summarise(N = sum(N),
            S = n_distinct(species),
            S_PIE = mobr::calc_PIE(N, ENS = TRUE),
            S_PIE2 = vegan::diversity(N, index = "invsimpson"),
            shannon = vegan::diversity(N, index = "shannon") # so this is not working
  )

SxS<- pivot_wider(exp_dat, id_cols = c(patch, time, treatment) ,
                  names_from = species,
                  values_from = N,
                  values_fill = 0)

library(vegan)
spec<- SxS[ , -c(1:3) ]
env<- SxS[ , c(1:3) ]

env<- cbind(env,
            shan = exp(diversity(spec, index = "shannon")),
            S = specnumber(spec),
            S_PIE2 = diversity(spec, index = "invsimpson"))

ggplot(data = subset(env, time>500), aes(x = time, y = S, color = treatment, group = patch ))+
  geom_point()

env %>% filter(time == -300) %>%
  ggplot()+
  geom_boxplot(aes(x = treatment, y = shan))



### This is the code that I wrote during the lecture:

highPiv<- sim_high_disp$dynamics.df %>%
  filter(time == 500)%>%
  pivot_wider(sim_equal_compDat,
              id_cols = patch,
              names_from = species,
              values_from = N)

lowPiv<- sim_low_disp$dynamics.df %>%
  filter(time == 500) %>%
  pivot_wider( id_cols = patch,
               names_from = species,
               values_from = N)
#stablePiv2<- reshape2::dcast(subset(sim_stabil_comp$dynamics.df, time == 0),
#                             patch~ species, value.var = "N" )


#alpha
compare<- rbind(
  data.frame(S = specnumber(highPiv[, -1]),
             treatment = "high",
             PIE = diversity(highPiv[, -1], index = "invsimpson"),
             ENSshan = exp(diversity(highPiv[, -1], index = "shannon"))),
  data.frame( S= specnumber(lowPiv[, -1]),
              treatment = "low",
              PIE = diversity(lowPiv[, -1], index = "invsimpson"),
              ENSshan = exp(diversity(lowPiv[, -1], index = "shannon"))))

par(mfrow = c(2,2))
boxplot(compare$S ~ compare$treatment, ylim = c(0,9))
boxplot(compare$PIE ~ compare$treatment, ylim = c(0,7))
boxplot(compare$ENSshan ~ compare$treatment, ylim = c(0,8))

summary(lm(compare$S~ compare$treatment))
summary(lm(compare$PIE~ compare$treatment))


# gamma
specnumber(colSums(highPiv[, -1]))
specnumber(colSums(lowPiv[, -1]))

diversity(colSums(highPiv[, -1]), index = "invsimpson")
diversity(colSums(lowPiv[, -1]), index = "invsimpson")

rarecurve(highPiv[, -1], main = "high")
rarecurve(lowPiv[, -1], main = "low")

mean(vegdist(highPiv[, -1]))
mean(vegdist(lowPiv[, -1]))


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







SxSgam<- pivot_wider(exp_dat, id_cols = c( time, treatment) ,
                     names_from = species,
                     values_from = N,
                     values_fill = 0,
                     values_fn = sum)

specGam<- SxSgam[ , -c(1:2) ]
envGam<- SxSgam[ , c(1:2) ]

envGam<- cbind(envGam,
               shan = exp(diversity(specGam, index = "shannon")),
               S = specnumber(specGam),
               S_PIE = mobr::calc_PIE(specGam, ENS = TRUE),
               S_PIE2 = diversity(specGam, index = "invsimpson"))

ggplot()+
  geom_point(data = env, aes(x = time, y = S, color = treatment, group = patch ), alpha = 0.5, size = 0.1)+
  geom_line(data = envGam, aes(x = time, y = S, color = treatment ))




#########################################################################################
library(mobr)

mob_in<- make_mob_in( subset(spec, env$time == -200), cbind(subset(env, time == -200), meta_landscape) , coord_names = c('x','y') )

mob_in


plot_rarefaction(mob_in, 'treatment', ref_level = 'equalizing', 'sSBR', lwd = 4)

oldpar <- par(no.readonly = TRUE)
par(mfrow = c(1, 2))
plot_rarefaction(mob_in, 'patch', '1', 'IBR',
                 pooled=TRUE, leg_loc='bottomright')
plot_rarefaction(mob_in, 'treatment', 'equalizing', 'IBR',
                 pooled=TRUE, leg_loc='bottomright')



plot_rarefaction(mob_in, 'treatment', 'equalizing', 'IBR', pooled = FALSE, lwd = 2,
                 leg_loc = 'topright')
plot_rarefaction(mob_in, 'treatment', 'equalizing', 'IBR', pooled = TRUE, lwd = 4,
                 leg_loc = NA)


# rank abundance

plot_abu(mob_in, 'treatment', type = 'rad', pooled = FALSE, log = 'x')
plot_abu(mob_in, 'treatment', type = 'rad', pooled = TRUE , log = 'x')
par(oldpar)


# 2-scale analysis
stats <- get_mob_stats(mob_in, group_var = "treatment", ref_level = 'equalizing',
                       n_perm = 199)

plot(stats, 'S')
plot(stats, 'N')
plot(stats, 'S_n')
plot(stats, 'S_PIE')

# continuous scale : Looking at the DIFFERENCE between the

deltaS = get_delta_stats(mob_in, 'treatment', ref_level='equalizing',
                         type='discrete', log_scale=TRUE, n_perm = 199)


plot(deltaS, stat = 'b1', scale_by = 'indiv', display='S ~ effort')
plot(deltaS, stat = 'b1', scale_by = 'indiv', display='stat ~ effort')



# iNEXT
com<- list(com1=X1, com2=X2, gamma = X3)
out <- iNEXT(com, q=c(0, 1, 2), datatype="abundance")
ggiNEXT(out, type=1, facet.var="order")









site1<- c(1, 3, 5, 7,  89, 6, 56, 3, 9)
site2 <-c(5, 5, 3, 23,  0, 12, 0, 6, 0)
site3 <-c(0 ,0, 1, 1,   0, 0,  0, 0, 0 )

allsites<- rbind(site1, site2, site3)

specnumber(site1)
specnumber(allsites)

sum(site1)
sum(site2)

sum(allsites)
rowSums(allsites)


diversity(allsites, index = "shannon")
diversity(allsites, index = "simpson")
diversity(allsites, index = "invsimpson")

rarefy(allsites, 2)

rarecurve(allsites)

vegdist(allsites, method = "jaccard")



site4<- c(100,1,1,1,1,1,1)
site5 <- c(3,3,3,3,3,3,3)
site6 <- c(1,1,1,1,1,1,1)

diversity(site5, index = "invsimpson")
rarecurve(rbind(site4, site5, site6))
