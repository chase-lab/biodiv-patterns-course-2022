#################################################
##  Studying beta-diversity: ecological variation partitioning
######## Wubing Xu; wubing.xu@idiv.de ##############

## install and load packages
install.packages("vegan")
install.packages("ade4")
install.packages("devtools")
library(devtools)
install_github("r-forge/sedar/pkg/PCNM/")

library(vegan)
library(ade4)
library(PCNM)


## as an example, we will use the mite dataset provided by the vegan package
# this dataset provide species composition, environmental conditions and spatial information of 70 sampling sites
# the community composition data table contain abundance of 35 species in 70 sites
data(mite)
# five environmental variables in the sampling sites
# SubsDens (Substrate density), WatrCont (Water content), Substrate (7 unordered classes), Shrub (3 ordered classes), Topo (Microtopography)
data(mite.env)
 # geographic coordinates of sites
data(mite.xy)


################
## quick explore the datasets

dim(mite) # 70 rows (sites) and 35 columns (species)
mite[1:4, 1:5]

# number of species (richness) in each sites
site.rich <- rowSums(mite>0)
hist(site.rich, xlab = "species richness", ylab = "Number of sites", main = "Species richness")

# number of sites occupied by each species
spe.nsites <- colSums(mite>0)
hist(spe.nsites, xlab = "Number of occupied sites", ylab = "Number of species", main = "Species occurrences")

# spatial distributions of sites, with the size of points proportional to richness
plot(mite.xy, asp = 1, cex = 2*site.rich/max(site.rich),
     pch = 21, bg = "skyblue", xlab = "X coordinate", ylab = "Y coordinate", main = "Map of species richness")

# Spatial distribution of four environmental conditions
par(mfrow = c(2, 2))
plot(mite.xy, asp = 1, cex = 2*mite.env$SubsDens/max(mite.env$SubsDens),
     pch = 21, bg = "skyblue", xlab = "X coordinate", ylab = "Y coordinate", main = "Substrate density")

plot(mite.xy, asp = 1, cex = 2*mite.env$WatrCont/max(mite.env$WatrCont),
     pch = 21, bg = "skyblue", xlab = "X coordinate", ylab = "Y coordinate", main = "Water content")

plot(mite.xy, asp = 1, cex = 1,
     pch = 21, bg = rainbow(3)[mite.env$Shrub], xlab = "X coordinate", ylab = "Y coordinate", main = "Shrub")
legend("topright", pch =16, col = rainbow(3), legend = c("None", "Few", "Many"), cex = 0.7)

plot(mite.xy, asp = 1, cex = 1,
     pch = 21, bg = rainbow(2)[mite.env$Topo], xlab = "X coordinate", ylab = "Y coordinate", main = "Microtopography")
legend("topright", pch =16, col = rainbow(2), legend = c("Blanket", "Hummock"), cex = 0.7)
par(mfrow = c(1, 1))



####################
## calculate dissimilarity of species compositions, and spatial and environmental distance

# Bray-Curtis dissimilarity: based on species abundance
mite.bc <- vegdist(mite, method = "bray", binary = FALSE)

# Jaccard dissimilarity: based on species presence/absence
mite.jac <- vegdist(mite, method = "jaccard", binary = TRUE)

# Sorensen dissimilarity (Bray-Curtis for presence/absence data)
mite.sor <- vegdist(mite, method = "bray", binary = TRUE)

# compare different beta indexes
# Bray-Curtis vs. Jaccard; strongly correlated, but a lot of variation
plot(mite.bc, mite.jac, xlab = "Bray-Curtis dissimilarity", ylab = "Jaccard dissimilarity")
abline(0, 1, col = "blue", lwd = 2)

# Sorensen vs. Jaccard; these two index can be derived each other
plot(mite.sor, mite.jac, xlab = "Sorensen dissimilarity", ylab = "Jaccard dissimilarity")

# spatial distance
mite.spadist <- dist(mite.xy)

# environmental distance. For simple illustration, we only used the continuous variables
# standardize different variables: center and scale variables: x_scaled = (x - mean)/sd
mite.envdist <- dist(scale(mite.env[, 1:2]))

# relationship between beta-diversity and spatial and environmental distance
# the relationship seems are non-linear. We don't consider the non-linear relationship here for simplification.
# This is why RDA is a better method for variation partitioning for most cases
plot(mite.bc ~ mite.spadist, xlab = "Spatial distance", ylab = "Bray-Curtis dissimilarity")
abline(lm(mite.bc ~ mite.spadist), lwd = 2, col = "blue")

plot(mite.bc ~ mite.envdist, xlab = "Environmental distance", ylab = "Bray-Curtis dissimilarity")
abline(lm(mite.bc ~ mite.envdist), lwd = 2, col = "blue")


####################
## Variation partitioning of beta-diversity between sites

## calculate different fractions of variation explained :
# regress beta-diversity against spatial and environmental distance. fraction = a+b+c
lm1 <- lm(mite.bc ~ mite.spadist + mite.envdist)

# beta-diversity against spatial distance. fraction =  a+b
lm2 <- lm(mite.bc ~ mite.spadist)

# beta-diversity against environmental distance. fraction =  b+c
lm3 <- lm(mite.bc ~ mite.envdist)

# get the fractions of variance explained. Use the adjusted r.squared. The r.squared is affected by the number of predictors
abc <- summary(lm1)$adj.r.squared
ab <- summary(lm2)$adj.r.squared
bc <- summary(lm3)$adj.r.squared

a <-  abc - bc
b <- ab + bc - abc
c <- abc - ab
d <- 1 - abc
a; b; c; d

# use function varpart to perform variation partitioning based on partial regression
# check whether the fractions are consistent with those that we calculated manually
mite.varp.regression <- varpart(as.vector(mite.bc), as.vector(mite.spadist), as.vector(mite.envdist))
mite.varp.regression

# show the fractions in Venn's diagram
plot(mite.varp.regression, digits = 2, Xnames = c("Spa.", "Env."), bg = c('blue', 'darkred'))


# test the significance of different fractions
# use the function mantel to test the significance of correlation between two dissimilarity matrices
# fraction a + b: beta-diversity vs. spatial distance
mantel(mite.bc, mite.spadist, method="pearson", permutations=999)

# fraction b + c: beta-diversity vs. environmental distance
mantel(mite.bc, mite.envdist, method="pearson", permutations=999)

# use the function mantel.partial to test the significance of correlation between two dissimilarity matrices while controlling for the third matrix
# fraction a
mantel.partial(mite.bc, mite.spadist, mite.envdist, method="pearson", permutations=999)
# fraction c
mantel.partial(mite.bc,  mite.envdist, mite.spadist, method="pearson", permutations=999)



###################
## Variation partitioning of comunity composition data table

# a) generate PCNM spatial eigenfunctions to represnet spatial structur across spatial scales
mite.spadist <- dist(mite.xy)
mite.PCNM.auto <- PCNM(mite.spadist)

# number of positive eigenvalus
length(mite.PCNM.auto$values)

# number of PCNM variables has significant positive spatial auto-correlations
head(mite.PCNM.auto$Moran_I)
table(mite.PCNM.auto$Moran_I$Positive)

# select the PCNM variables with positive spatial auto-correlations
select <- which(mite.PCNM.auto$Moran_I$Positive == TRUE)
mite.pcnm <- as.data.frame(mite.PCNM.auto$vectors)[,select]
dim(mite.pcnm)

# plot some PCNM variables using the function s.value from package "ade4"
op <- par(mfrow=c(2,2))
par(mfrow=c(2,2))
somePCNM <- c(1, 2, 10, 20)
for(i in 1:length(somePCNM)){
  s.value(mite.xy, mite.pcnm[,i], sub=somePCNM[i], csub=2)
}
par(op)


# b) selection of variable in RDA

# species dataset needs to be standardized in RDA
# Hillinger transformation is a common used: divided by the rowsums (total abundance in a site), and then calculate square root
mite.hel <- decostand(mite, "hellinger")

# selection of environmental variable
rda.env1 <- rda(mite.hel ~ 1, data = mite.env) # the model with only intercept
rda.env2 <- rda(mite.hel ~ ., data = mite.env) # the model include all environmental variables
rda.env.step.fore <- ordistep(rda.env1, scope = formula(rda.env2), direction = "forward", perm.max = 200, pstep = 999)
rda.env.step.fore$anova # all variables are selected

# selection of spatial PCNM variable
rda.spa1 <- rda(mite.hel ~ 1, data = mite.pcnm) # the model with only intercept
rda.spa2 <- rda(mite.hel ~ ., data = mite.pcnm) # the model include all PCNM variables
rda.spa.step.fore <- ordistep(rda.spa1, scope = formula(rda.spa2), direction = "forward", perm.max = 200, pstep = 999)
rda.spa.step.fore$anova # some variables are selected

# extract the selected pcnm variables
mite.pcnm.sel_names <- gsub(pattern = "+ ","", rownames(rda.spa.step.fore$anova),fixed = TRUE)
mite.pcnm.sel <- mite.pcnm[, names(mite.pcnm) %in%  mite.pcnm.sel_names]
dim(mite.pcnm.sel)
head(mite.pcnm.sel)


# c) calculate different fractions of variation explained using RDA

# use both spatial and environmental variables: fraction = a+b+c
rda.all <- rda(mite.hel, cbind(mite.pcnm.sel, mite.env))
# use only spatial variables: fraction =  a+b
rda.spa <- rda(mite.hel, mite.pcnm.sel)
# use only environmental variables. fraction =  b+c
rda.env <- rda(mite.hel, mite.env)

# get the fractions of variance explained. Use the adjusted r.squared
abc <- RsquareAdj(rda.all)$adj.r.squared
ab <- RsquareAdj(rda.spa)$adj.r.squared
bc <- RsquareAdj(rda.env)$adj.r.squared

a <-  abc - bc
b <- ab + bc - abc
c <- abc - ab
d <- 1 - abc
a; b; c; d

# use function varpart to perform variation partitioning based on RDA
# check whether the fractions are consistent with those that we calculated manually
mite.varp.rda <- varpart(mite.hel, mite.pcnm.sel, mite.env)
mite.varp.rda

# show the fractions in Venn's diagram
plot(mite.varp.rda, digits = 2, Xnames = c("Spa.", "Env."), bg = c('blue', 'darkred'))

# test the significance of different fractions
# the fractions a+b+c
anova.cca(rda.all)
# the fractions a+b
anova.cca(rda.spa)
# the fractions b+c
anova.cca(rda.env)
# the fraction a
rda.spa.env <- rda(mite.hel, mite.pcnm.sel, mite.env)
anova.cca(rda.spa.env)
# the fraction c
rda.env.spa <- rda(mite.hel, mite.env, mite.pcnm.sel)
anova.cca(rda.env.spa)



###################
# In some empirical studies, there are many environmental variables, which can be strongly correlated.
# PCA are usually used to extract few uncorrelated principle components (new variables) to represent most of variation in raw variables
# we didn't perform PCA for mite dataset due to few environmental variables in this dataset.
# here we use the dataset "varechem" (14 soil variables, 24 sites) as an example to know how perform PCA

data(varechem)
summary(varechem)

# performs a principal components analysis
pca_model <- princomp(scale(varechem))
summary(pca_model) # the first three components totally explain 71.9% of variation

# get the values of principle components for each site
head(pca_model$scores[,1:3])

# the loading of raw environmental variables for each principle component
# they are actually the standardized coefficients from regression of principle component against environmental variables
head(pca_model$loadings[,1:3])

# check whether the loading are the standardized coefficients
y <- pca_model$scores[,1]
x <- scale(varechem)
coef(lm(y~x))[-1]
pca_model$loadings[,1]
