##################################################################
##     R Practical for MLU Block Course: Modelling species      ##
##      distribution and biodiversity patterns - Part II        ##
##                                                              ##
##          ****  BONUS: VARIANCE PARTITIONING *****            ##
##                                                              ##
##                 Wednesday, January 12, 2022                  ##
##               Joint Species Distribution Models              ##
##                    Dr. Kimberly Thompson                     ##
##                                                              ##
##     R code adapted from:                                     ##
##     Ovaskainen, Otso. Joint species distribution modelling:  ##
##        with applications in R. Cambridge : Cambridge         ##
##        University Press, 2020.                               ##
##################################################################


# Call the Hmsc package (needs to be done everytime you start a new R session)
library(Hmsc)

# To disentangle the effects of the environmental and spatial predictors in our model
# of C. monedula abundance, we will build three models - a full model (mFULL: same as in Part IV
# of the main R exercise), a model that includes only environmental predictors but no
# spatial random effect (mENV), and a model that includes only a spatial random effect but
# no environmental predictors (mSPACE).

# The first part of this procedure is identical to Part Iv in the main R exercise (but here
# I have reproduced only the essentail lines of code, without as much explanation.)



# Download the bird data (data.csv) and the grid data (grid_1000.csv) from the class github page if you haven't already.

# Set your working directory to where you stored the folder.
setwd()

# Load the bird data
da_all = read.csv("data.csv", header = TRUE)

# We will examine Corvus monedula during 2014, so first we subset our data frame to show us
# only this year.
da = droplevels(subset(da_all, Year==2014))

# We will examine the relationship between distribution and Habitat and April/May temperature.
# So we create our data frame of predictor variables:
XData = data.frame(hab = da$Habitat, clim = da$AprMay)

# We also need to make sure that R sees our habitat variable as a factor.
XData$hab = as.factor(XData$hab)

# Create the matrix for the response using only the Corvus monedula data. Again, our matrix
# will only have 1 column.
Y = as.matrix(da$Corvus_monedula)

# We additionally need to create a matrix of our x and y spatial coordinates.
xy = as.matrix(cbind(da$x, da$y))

# Define the study design. This time sampling units are grouped into Routes.
studyDesign = data.frame(route = as.factor(da$Route))

# Rename the rows of the xy matrix according to the Route of each sampling unit
rownames(xy) = studyDesign[, 1]

# Define the random spatial effect
rL = HmscRandomLevel(sData = xy)

# Now since we have more than one predictor, we will explicitly define our XFormula before 
# we build the model object. Here we specify a polynomial term for temperature aka climate,
# which then automatically also includes the linear term.
XFormula = ~ hab + poly(clim, degree = 2, raw = TRUE)


#### Now is where the code diverges from Part IV of the main R exercise.


####################################################
##   Creating the model objects and fitting them  ##
####################################################


# Build the model objects
mFULL = Hmsc(Y = Y, XData = XData, XFormula = XFormula, distr = "lognormal poisson",
             studyDesign = studyDesign, ranLevels = list(route = rL))

mENV = Hmsc(Y = Y, XData = XData, XFormula = XFormula, distr = "lognormal poisson")

mSPACE = Hmsc(Y = Y, XData = XData, XFormula = ~ 1, distr = "lognormal poisson",
              studyDesign = studyDesign, ranLevels = list(route = rL))


# Specify the parameters for model fitting
nChains = 4
samples = 250
thin = 10
transient = round(0.5*samples*thin)
verbose = round(0.5*samples*thin)

# (Keep in mind that we are now fitting 3 models, which will increase the computation time)


# We can fit all three models by looping the MCMC sampler over them
# Create a list of the models
models = list(mFULL, mENV, mSPACE)

# Loop the sampler over the 3 models
for (i in 1:3) {
  models[[i]] = sampleMcmc(models[[i]], thin = thin, samples = samples, 
                          transient = transient, nChains = nChains, verbose = verbose,
                          initPar = "fixed effects")
}

# While the code is not reproduced here, you should check each model for convergence 
# following the steps outlined in the main R exercise (lines 961-978 in Part IV).


#####################################################
##   Evaluate the explanatory power of each model  ##
#####################################################

# We again, create a list and a loop to examine the fit of our model

# Create an empty list where the results of our fit test will be stored
MF = list()

# Loop the evaulation of model fit over the 3 models
for (i in 1:3) {
  preds = computePredictedValues(models[[i]], expected = FALSE)
  MF[[i]] = evaluateModelFit(hM = models[[i]], predY = preds)
}

# Examine the model fit
MF.table <- data.frame(Model = c("Model Full", "Model ENV", "Model SPACE"),
                                 O.AUC = c(MF[[1]]$O.AUC, MF[[2]]$O.AUC, MF[[3]]$O.AUC),
                                 C.SR2 = c(MF[[1]]$C.SR2, MF[[2]]$C.SR2, MF[[3]]$C.SR2))

MF.table

# To recap, we look at AUC for discriminating presences and absences (O.AUC, where the O
# stands for 'observed'), and the pseudo-R squared for discriminating abundances 
# conditional on presence (C.SR2, where the C stands for conditional: in other words
# this assesses how well the model did at predicting abundances for the instances
# where C. monedula was present)

# We see that the 3 models have a relatively high AUC and hence a relatively high
# power to discriminate presences from absences. However, this doesn't necessarily mean
# that the models are great, since high explanatory power can also arise from overfitting.
# In contrast, the models have a rather low ability to discriminate among abundance 
# variation, even for the same data that were used to fit the models. We examine this 
# further through variance partitioning.


#####################################################
##             Variance Partitioning               ##
#####################################################

# We will perform variance partitioning among the fixed effects related to habitat and
# climate, and the random effects.

# Examine the design matrix X that Hmsc has constructed from the XData and XFormula.
round(head(models[[1]]$X), 2)

# We next group the explanatory variables. In the design matrix X, columns 2-5 relate
# to habitat type, and the columns 6 and 7 to climate.
# The first column models the intercept, which does not explain any variance so we can
# group this column arbitrarily. Here we will group it with the habitat variable.
groupnames = (c("habitat", "climate"))
group = c(1 , 1, 1, 1, 1, 2, 2)


# Now we can call the function computeVariancePartitioning in a loop over the models
# Create an empty list where the results of our fit test will be stored
VP = list()

# We loop only over the models mFULL AND mENV because model mSPACE includes only the 
# spatial random effect and thus all variance is explaine by the random part. In the 
# model mENV we have included only the environmental predictors with no random effect
# (so the spatial random effect doesn't explain anything in that model).
for (i in 1:2) {
  VP[[i]] = computeVariancePartitioning(models[[i]], group = group,
                                        groupnames = groupnames)
}

# Examine the results
VP.table <- data.frame(Model = c("Model Full", "Model ENV", "Model SPACE"),
                       habitat = c(VP[[1]]$vals[1], VP[[2]]$vals[1], 0.000),
                       climate = c(VP[[1]]$vals[2], VP[[2]]$vals[2], 0.000),
                       random = c(VP[[1]]$vals[3], 0.000, 1.000))

VP.table

# In both FULL and ENV, the climatic variables explain the highest amount of the variation,
# but the habitat variable is important as well. In FULL, the random part explains 
# very little, suggesting that we have included all relevant predictors.


