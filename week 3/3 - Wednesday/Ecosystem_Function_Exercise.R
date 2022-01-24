

library(tidyverse)
library(ggplot2)
devtools::install_github("ctkremer/priceTools")
library(priceTools)
library(patchwork)

# Paper: https://doi.org/10.1111/ele.13566
# Code for the paper: https://github.com/emma-ladouceur/SeedAdditionSynthesis
# Data from the paper (Download data here): https://figshare.com/articles/dataset/Reducing_dispersal_limitation_via_seed_addition_increases_species_richness_but_not_aboveground_biomass/12319682


# Download data from above links
# Set your working directory
setwd("~/Downloads/")
# species level data from 12 seed addition experiments
sp <- read.csv("SeedAdd_Sp_level.csv", header= TRUE)


head(sp)
colnames(sp)

# Rename Experiment short names to match the paper
 sp <- sp %>% mutate( Experiment = case_when(Experiment_ == "ASGA_Michigan" ~ "Michigan",
                                             Experiment_ == "California_Invade" ~ "California.1",
                                             Experiment_ == "California_Prop_Limi" ~ "California.2",
                                             Experiment_ == "CCR_04" ~ "Cedar.Creek.4",
                                             Experiment_ ==  "CCR_093" ~ "Cedar.Creek.93",
                                             Experiment_ == "Germany_Montane" ~ "Montane",
                                             Experiment_ == "Halle" ~ "Halle",
                                             Experiment_ ==  "Jena" ~ "Jena",
                                             Experiment_ == "Jena2" ~ "Jena.2",
                                             Experiment_ ==  "Kansas_KUFS_LTER_Hay_Meadow_Exp_2" ~ "Kansas.Hay.Meadow",
                                             Experiment_ == "Kansas_KUFS_LTER_Old_Field_Exp_1" ~ "Kansas.Old.Field",
                                             Experiment_ == "Texas_Temple_Prarie" ~ "Texas.Temple.Prairie")) %>%
   mutate( trt = case_when( trt == "Seed" ~ "Seeds",
                            TRUE ~ as.character(trt)))


# select only the plot-level data
plot <- sp %>% select(Experiment, site, block, plot, yr.trt, trt,
                      seed.rich, biomass.plot, biomass.measure, rich.plot) %>%
  distinct() %>% mutate(Treatment = trt) %>% select(-trt)


head(plot)


# let's have a look at the data
ggplot(data = plot, aes(x = rich.plot, y = biomass.plot)) +
  geom_point(aes(color = Experiment
                 ) , alpha = 0.5, size = 1.2, position = position_jitter(width = 0.95, height = 0.95) ) +
 geom_smooth( aes(linetype= Treatment), color = "black", method = lm, se=FALSE, fullrange=TRUE ) +
  labs(x = 'Realised Richness',
       y = expression(paste('Biomass (g/',m^2, ')')), title= "",
       subtitle="") +
  scale_colour_manual(values = c( "#EE0011FF" , "#EC579AFF", "#15983DFF", "#149BEDFF", "#0C5BB0FF",
                                  "#8F2F8BFF", "#F9B90AFF", "#16A08CFF" ,"#6A7F93FF","#FA6B09FF","#A1C720FF","#9A703EFF" )) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     strip.background = element_rect(colour="black", fill="white"),legend.position="bottom",
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank())



ggplot(data = plot, aes(x = Treatment, y = rich.plot)) +
 geom_point(position = position_jitter(0.2), alpha= 0.5, color = "#C0C0C0") +
  geom_boxplot(alpha = 0.2) +
  labs(y = "Plot Richness") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     strip.background = element_rect(colour="black", fill="white"),legend.position="bottom",
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())




# Now we will go back to the species-level data to examine things a bit more closely

# identify grouping variables
sp$trt<-as.factor(as.character(sp$trt))


treat.vars <- c("trt")
grouped.seed.dat <- sp %>% group_by(.dots=c(treat.vars))

# use the question mark to find out more about a function
#?pairwise.price

res.seed <- pairwise.price(grouped.seed.dat, species = "species", func = "biomass.sp")

# Create a single column keeping track of the paired set of seeding treatments & other grouping variables:
pp.s <- res.seed

pp.s <- group.columns(pp.s , gps=c(treat.vars), drop=T)

head(pp.s)

# prune to only a comparison we want
pp.s <- pp.s %>% filter( trt %in% 'Control Seeds')


# plot!

# BEF style
s1<-leap.zig(pp.s,type='bef',standardize = FALSE,raw.points = T)+
  annotate("text", x = mean(pp.s$x.rich), y = mean(pp.s$x.func),
           label = "*",size=8)+ggtitle('Control Seeds')+theme_classic()


# CAFE style with gains and losses
s2<-leap.zig(pp.s,type='cafe',standardize = FALSE,raw.points = T)+
  annotate("text", x = mean(pp.s$x.rich), y = mean(pp.s$x.func),
           label = "*",size=8)+ggtitle('Control Seeds')+theme_classic()

(s1  + s2)


# use all the raw data


treat.vars <- c("trt", "block",  "plot", "seed.rich")
head(sp)

grouped.seed.dat <- sp %>% filter(Experiment == "Cedar.Creek.93") %>%
  filter(yr.trt == "3") %>%
  group_by(.dots=c(treat.vars))

View(grouped.seed.dat)

grouped.seed.dat$block <- as.factor(grouped.seed.dat$block)
levels(grouped.seed.dat$block)
grouped.seed.dat$yr.trt <- as.factor(grouped.seed.dat$yr.trt)
levels(grouped.seed.dat$yr.trt)

# use the question mark to find out more about a function
#?pairwise.price

res.seed <- pairwise.price(grouped.seed.dat, species = "species", func = "biomass.sp")

# Create a single column keeping track of the paired set of seeding treatments & other grouping variables:
pp.s <- res.seed

pp.s <- group.columns(pp.s , gps=c(treat.vars), drop=T)

head(pp.s)

# prune to only a comparison we want
pp.wrangle <- pp.s %>%
  # unite(trt, trt.x , trt.y) %>%
  # unite(block, block.x , block.y) %>%
  # unite(yr.trt, yr.trt.x , yr.trt.y) %>%
  filter( trt %in% c('Control Seeds'),
  block %in% c('27 27','28 28',
              '29 29','30 30','31 31','32 32','33 33','34 34','35 35','36 36','37 37', '38 38','39 39','40 40','41 41',
              '42 42','43 43','44 44','45 45','46 46','47 47','48 48','49 49','50 50','51 51','52 52','53 53','54 54','55 55',
              '56 56') ) %>%
  mutate(s.loss = (x.rich - c.rich), # calculate species loss
         s.gain = (y.rich - c.rich)) # calculate species gain


head(pp.wrangle)

# plot!
div.5 <- pp.wrangle %>% filter(seed.rich %in% '0 5')

# CAFE style with gains and losses
s5 <- leap.zig(div.5,type='cafe',standardize = FALSE, raw.points = F)+
  annotate("text", x = mean(div.5$x.rich), y = mean(div.5$x.func),
           label = "*",size=8)+ggtitle('Control Seeds')+theme_classic()
s5



div.54 <- pp.wrangle %>% filter(seed.rich %in% '0 54')


s54 <- leap.zig(div.54,type='cafe',standardize = FALSE, raw.points = F)+
  annotate("text", x = mean(div.54$x.rich), y = mean(div.54$x.func),
           label = "*",size=8)+ggtitle('Control Seeds')+theme_classic()
s54


ggplot(data= pp.wrangle, aes(x = seed.rich, y = SG)) +
  geom_boxplot(aes(color = seed.rich))

ggplot(data= pp.wrangle, aes(x = seed.rich, y = CDE)) +
  geom_boxplot(aes(color = seed.rich))

ggplot(data= pp.wrangle, aes(x = seed.rich, y = s.gain)) +
  geom_boxplot(aes(color = seed.rich))

