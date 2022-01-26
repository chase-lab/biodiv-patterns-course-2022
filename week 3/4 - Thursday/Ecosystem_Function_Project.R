####packages#####
library(tidyverse)


# rm(list = ls())
#####data####
# Paper:  https://doi.org/10.1002/ece3.4331
# Data from the paper (Download data here):
#https://datadryad.org/search?utf8=%E2%9C%93&q=Grazing+effect+on+grasslands+escalated+by+abnormal+precipitations+in+Inner+Mongolia

# Download data from above links
# Set your working directory
#setwd("~/Downloads/")
data <- read.csv("Table1.csv", header= TRUE)
####packages#####
library(tidyverse)


# rm(list = ls())
#####data####
# Paper:  https://doi.org/10.1002/ece3.4331
# Data from the paper (Download data here):
#https://datadryad.org/search?utf8=%E2%9C%93&q=Grazing+effect+on+grasslands+escalated+by+abnormal+precipitations+in+Inner+Mongolia

# Download data from above links
# Set your working directory
#setwd("~/Downloads/")
data <- read.csv("Table1.csv", header= TRUE)

data %>% head
data %>% names

# renaming columns to ease the work
data<- data %>% rename(land_use = "Land_use" ,
                       land_use1 = "Land_use1",
                       years = "Years"    ,
                       month = "Month"    ,
                       plots = "Plots"    ,
                       per.bunchgrasses = "PB...."   ,
                       rhizome.grass = "PR...."   ,
                       per.forbs = "PF...."   ,
                       annuals.biennials = "AB...."   ,
                       richness = "Species.richness..no...species...m.2.",
                       stand.density = "Stand.density..no...plants...m.2.",
                       canopy.height = "Canopy.height..cm.",
                       spp.diversity = "Species.diversity..H..",
                       func.diversity = "Functional.diversity..Rao.s.Q.",
                       biomass = "Aboveground.biomass..g...m.2.")

# and now it is up to you to play with the data and the BEF relationships


#  -----------------------------------------------------------------
# Regarding our BEF debate, what can we learn with this study?

