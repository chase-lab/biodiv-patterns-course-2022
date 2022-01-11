##################################################################
##     R Practical for MLU Block Course: Modelling species      ##
##      distribution and biodiversity patterns - Part II        ##
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

# We will use the Hmsc package (Hierarchical Modelling of Species Communities)
# to explore ways of analyzing multiple species that are affected by not only
# their environment, but also by dispersal, biotic interactions, and stochasticity.

# Documentation for the Hmsc package can be found at
# https://cran.r-project.org/web/packages/Hmsc/Hmsc.pdf


# Install the Hmsc package (installation only needs to be done once)
install.packages("Hmsc")

# Call the Hmsc package (needs to be done everytime you start a new R session)
library(Hmsc)

# We will also need a few other packages that you may also need to install.
install.packages("MASS")
install.packages("sp")
install.packages("gridExtra")
install.packages("abind")
install.packages("ggplot2")

library(MASS)
library(sp)
library(gridExtra)
library(abind)
library(ggplot2)



##################################################################
##    Part I: A Univariate Example with Simulated Data          ##
##################################################################

# We will simulate data for a single covariate (i.e. predictor) x on n=50 sampling
# units (i.e. sites). We then will analyze this date with three different models: 
# a linear (normal) model, a probit model (for presence absence data), and a logpoisson 
# normal model (for abundance data). 


#################################
##   Step 1: Simulating Data   ##
#################################

# Set the random number seed so that our results (in which we will generate
# some random numbers) are reproducible
set.seed(1)

# Specify the number of sampling units
n = 50

# Generate random values of the predictor variable x, one value for each of our
# sampling units. 
# Predictor variable x represents some aspect of the environment that we are measuring
# (or in this case simulating). An example would be the mean temperature at a site. We assume
# that this environmental condition varies by site, but that all species at a given site 
# experience the same environmental condition, though importantly species will not respond to 
# the condition in the same way (--> distributional differences).
# Note that the default mean and standard deviation of the function
# 'rnorm' are 0 and 1, respectively.
x = rnorm(n)

# We know the default mean and standard deviation of the rnorm function, but we can confirm this
# to illustrate how the rnorm function works.
mean(x)     # This should be approximately 0.
sd(x)       # This should be approximately 1.
plot(x)

# Define the intercept
beta1 = 0

# Define the slope
beta2 = 1

# Define the general form of our linear predictor (you should recognize this as
# an equation following the basic y=mx+b format where b(and beta 1) = intercept and
# m (and beta2) = slope)
L = beta1 + beta2*x

# Complete the linear (normal) model by adding a term that corresponds to the residuals
# Since we have n=50 sampling units, we generate a residual value for each unit.
# Recall from the above comment that the default mean and standard deviation of the function
# 'rnorm' are 0 and 1, respectively.
y1 = L + rnorm(n)

# Generate the responses for the probit model (presence absence model) 
y2 = 1*((L + rnorm(n)) > 0)

# Note what we just did in the above line of code. We really just generated a new set of 
# responses for our linear normal model (L + rnorm(x)), asked R to tell us which of these values
# were greater than 0, and recoded them as 0 and 1 accordingly. We can break up this line of 
# code into smaller sections to illustrate this point.
random.response <- L + rnorm(n)
random.response
responses.greater.than.zero <- random.response > 0
responses.greater.than.zero
new.responses <- 1*responses.greater.than.zero
new.responses

# Remove these vectors since they were only used to illustrate the point and we won't need them
# for our analysis.
rm(random.response, responses.greater.than.zero, new.responses)

# Generate the responses for the lognormal Poisson model (abundance model)
y3 = rpois(n = n, lambda = exp(L + rnorm(n)))

# Note that the function rpois, which generates random values from the Poisson distribution, 
# takes 2 arguments. The first n, corresponds to the number of observations. In our case,
# we have defined n=50 sampling units. The second argument, lambda, refers to the
# estimated rate of events for the distribution; this is expressed as average events per period.
# In our case, what we are saying then is for each of our 50 sampling units (e.g. we observed 50 
# different forest sites for a specified amount of time) lambda is the rate at which we
# observed individuals. So if n=1 and lambda=2, in our first sampling unit we observed 2 
# individuals.
# Examine values of y3 you just generated to illustrate this point.
y3

# Now that we have generated our data, we can plot the predictor (x), and the responses 
# (y1, y2, and y3) to visualize our simulated data.
par(mfrow = c(1,3)) # sets up a graphing panel with 1 row and 3 columns
plot(x, y1, main = "Normal", ylab = "Response", xlab = "Explanatory Variable x") # plot 1  
plot(x, y2, main = "Probit", ylab = "Presence Absence", xlab = "Explanatory Variable x") # plot 2
plot(x, y3, main = "Lognormal Poisson", ylab = "Abundance", xlab = "Explanatory Variable x") # plot 3



#####################################################
##   Step 2: Building the Linear (Normal) Model    ##
#####################################################

# For a point of comparison, let's explore our normal model without the package Hmsc.

# Make a dataframe out of our predictor variable and our normally distributed response.
df = data.frame(x, y1)

# Define the linear model
m.lm = lm(y1 ~ x, data = df) 

# Look at the parameter estimates
summary(m.lm)

# You can see that the estimate for the intercept is approximately 0, while the estimate for the
# slope is approximately 1, as we specified when we generated the data. 

# How well does this linear model fit our simulated data? The R squared value in this summary is
# a good indicator of this since it represents the proportion of the variation that's explained
# by x.

par(mfrow = c(1,1)) # Returns plot space to default
plot(x, y1)
abline(lm(y1~x))




# Let's now conduct an analogous analysis with Hmsc.
# The inputs for Hmsc are a little different since the aim of the package is to model multiple
# species and multiple predictors. Consequently, we will convert our y1 and x vectors 
# accordingly.

# Convert y1 to a matrix. Of course, since we are examining only 1 hypothetical species, the 
# dimensions of the matrix will be 50x1 (50 sampling units (i.e. sites), 1 species). 
# With more species we will have more columns in our matrix.
Y = as.matrix(y1)

# Convert x to a dataframe. Similarly, with more predictors, the number of columns in our
# dataframe would increase.
XData = data.frame(x = x)

# Build the Hmsc model object. We are specifying here the response, the predictor, 
# how we want to relate them, and the distribution of the data.
m.normal = Hmsc(Y = Y, XData = XData, XFormula = ~x, distr = "normal")

# Note that while the lm function both constructs and fits the model with one line of code, 
# the Hmsc function only constructs the model object. 

# Hmsc uses Bayesian inference, and while there are many concpetual and statistical details
# we could delve into regarding the Bayesian framework, we will focus on only general 
# concepts and the details necessary to build the model today. Overall, it is important to 
# remember that our results for this type of analysis will not be a specific parameter 
# estimate as in the m.lm linear model we just ran, but rather a distributon (called the 
# posterior distribution), from which our parameter likely comes from.

# To create the posterior distribution, analyses that use a Bayesian framework make use of a 
# sampling algorithm which through repeated runs generates the values that in turn form the
# distribution. A common sampling algorithm is Markov Chain Monte Carlo (MCMC) sampling. This
# is also the algorithm used by the Hmsc package.

# We have built our Hmsc model object, but now we have to do the MCMC sampling.
# Before we can do this though, we have to specify some additional inputs 
# to tell R how we would like to analyze this model. These include:
# how many chains to sample (nChains)
# how many samples to obtain per chain (samples)
# how much thinning to apply (thin)
# what length of a transient (also called burn-in) to include (transient)
# how frequently we want to see the progress of this process (verbose)

# What do these inputs mean? We will discuss this a little further down in the code.

nChains = 2
thin = 5
samples = 1000
transient = 500*thin
verbose = 500*thin

# Now that we have specified these inputs, we can do the MCMC sampling.
m.normal = sampleMcmc(m.normal, thin = thin, samples = samples, transient = transient,
                      nChains = nChains, verbose = verbose)

# Once the sampling is complete, our model object m.normal includes the estimated parameters,
# just like our linear model object m.lm did when we used the lm function.

# To view our Hmsc model estimates, we extract the posterior distribution from the model object.
# This is done by using the convertToCodaObject function.
mpost = convertToCodaObject(m.normal)

# We can now view the parameter estimates using the summary function, and we are particularly 
# interested in the beta estimates (i.e. slope and intercept).
summary(mpost$Beta)

# Again, we see that the estimate of the intercept is close to 0 and the estimate of the slope
# is close to 1, as we expect given what we simulated.


# How well does the Hmsc model fit our simulated data? We can calculate R squared in the Hmsc 
# package with the following steps:

# Compute the predicted values of the posterior distribution
preds = computePredictedValues(m.normal)

# Calculate the model fit
MF = evaluateModelFit(hM = m.normal, predY = preds)

# View the R squared value
MF$R2

# Similar to the slope and the intercept, the R squared value given by Hmsc is also highly 
# consistent with that given by the lm function. BUT, the two approaches are not identical as 
# discussed above. The lm function applies a maximum-likelihood (i.e. frequentist) framework 
# (and thus yields confidence intervals), while the Hmsc package applies a Bayesian framework
# (and thus yields credible intervals).


#####################################################
##     Step 3: Checking the Hmsc linear model      ##
#####################################################

# How do we know whether the modeling process we just employed with Hmsc was applied correctly?
# With Hmsc (and more generally with any MCMC sampling method) we have to confirm that the MCMC
# chains converged. If they have converged, the samples provided by the MCMC chain will be 
# representative of the true posterior distribution (and thus the inference from the model can
# be trusted). 
# Lack of convergence can lead to biased parameter estimates as well as a biased view on the 
# amount of uncertainty in the estimates. In other words, without convergence we can not use any
# of the model results.

# We check for convergence both visually and quantitatively.

# For visual inspection, we use trace plots.
plot(mpost$Beta)

# The trace plots (i.e. the top left and bottom left plots) will help us get a better 
# understanding of the additional inputs we had to specify for the MCMC sampling.
# Recall that we specified nChains = 2. The black and red lines show our two independent MCMC 
#   chains.One reason for running multiple chains is that any individual chain might converge
#   toward one target, while another chain might converge elsewhere, which would be a problem. 
# Recall that we set transient = 500*thin where thin = 5. On the x-axis, our iterations do not 
#   start at 1 because we chose not to store the values before 2500 iterations. So transient 
#   refers to how many of the initial iterations are discarded.
# Recall that we specified samples = 1000 and thin = 5. This means that we are selecting 1000 
#   samples, but we are only going to record every 5th iteration. Consequently, we need 5000 
#   iterations total and this is why the maximum value on the x-axis is 7500.

# What are we looking for in these trace plots? 
#   In general, we look for a plot that shows random scatter around a mean value and that the 
#   chains mix well together (there is strong overlap between them throughout all iterations).
#   Here, our trace plots look good.


# Next look at the top right and bottom right plots. In our lm model, we visualized the modeled
# results by plotting the best fit line. Our results from Hmsc though are instead two 
# distributions, one for the intercept and one for the slope.


# Checking for convergence quantitatively in two ways. First with effective sample sizes and 
# second with potential scale reduction factors.
effectiveSize(mpost$Beta)

# For all parameters the effective sample sizes are very close to the actual sample size
# which is 2000 (1000 per chain). This indicates that there is very little autocorrelation among
# consecutive samples.

gelman.diag(mpost$Beta, multivariate = FALSE)$psrf

# psrf stands for potential scale reduction factors. Since for both parameters they are very 
# close to 1, this means that the two chains gave consistent results, which is what we want to 
# see.


#####################################################
##      Step 4: The Probit Model with Hmsc         ##
#####################################################

# Convert the vector of simulated response values for the probit model to a matrix
Y = as.matrix(y2)

# From the work we did in Step 2, our environmental predictor x is already formatted in a 
# dataframe called XData so we do not need to make any other adjustments to it.

# Build the Hmsc model object
m.probit = Hmsc(Y = Y, XData = XData, XFormula = ~x, distr = "probit")

# Conduct the MCMC sampling
# We will leave our inputs for the MCMC sampling (nchains, samples, etc.) the same as in Step 2.
m.probit = sampleMcmc(m.probit, thin = thin, samples = samples, transient = transient,
                      nChains = nChains, verbose = verbose)

# Extract the posterior distribution from the model object with the convertToCodaObject function.
mpost = convertToCodaObject(m.probit)

# Evaluate convergence visually with trace plots.
plot(mpost$Beta)

# Evaluate convergence quantitatively by examining effective samples sizes.
effectiveSize(mpost$Beta)

# Compared to the linear model, we observe that the effective sample size is somewhat smaller. 
# You may find this to be the case for one or both of your parameters. This is because achieving
# MCMC convergence is generally more challenging for non-normal models than for normal models.
# So this result for a non-normal model is ok.

# Evaluate convergence quantitatively by examining potential scale reduction factors.
gelman.diag(mpost$Beta, multivariate = FALSE)$psrf

# We again see that for both parameters the PSRF estimates are very close to 1, so we can feel
# confident that the two chains gave consistent results.


#############################################################
##  Step 5: Checking the Probit Model's explanatory power  ##
#############################################################

# For presence absence models, Hmsc provides 3 different metrics to evaluate how well the model
# explains the data.
# Root Mean Square Error (RMSE) : Measures the error of a model in predicting data by answering 
#       the question, "How far off should we expect our model to be on its next prediction?" 
#       RMSE measures the standard deviation of a typical observed value from our model's 
#       prediction. Units of RMSE are the same as the response variable in the model.
#       The output of presence absence models is typically in the form of probabilities of 
#       occurrence. So the RMSE value here is telling us by how much the probability of 
#       occurrence is likely to be over- or underestimated by our model.
# Area Under the ROC Curve (AUC) : AUC tells how much the model is capable of distinguishing 
#       presences and absences correctly. AUC ranges from 0 to 1. The higher the AUC, the better
#       the model is at predicting 0 classes as 0 and 1 classes as 1. An AUC of 0.5 means the 
#       model performs as well as by chance.
# Tjur R squared : Similar to R squared used in linear models, but adapted for generalized 
#       linear mixed models with binary outcomes (like the probit model). A Tjur R2 of 0 means 
#       the model performs as well as by chance, while a value of 1 indicated perfect 
#       discrimination between presences and absences.

# Compute the predicted values of the posterior distribution
preds = computePredictedValues(m.probit)

# Evaluate the model's explanatory power
evaluateModelFit(hM = m.probit, predY = preds)

# The output of presence absence models is typically in the form of probabilities of occurrence.
# So the RMSE value here is telling us by how much the probability of occurrence is likely to 
# be over- or underestimated by our model.


#####################################################
##  Step 6: The Lognormal Poisson Model with Hmsc  ##
#####################################################

# Convert the vector of simulated response values for the lognormal poisson model to a matrix
Y = as.matrix(y3)

# From the work we did in Step 2, our environmental predictor x is already formatted in a 
# dataframe called XData so we do not need to make any other adjustments to it.

# Build the Hmsc model object
m.lognormal.poisson = Hmsc(Y = Y, XData = XData, XFormula = ~x, distr = "lognormal poisson")

# Conduct the MCMC sampling
# We will leave our inputs for the MCMC sampling (nchains, samples, etc.) the same as in Step 2.
m.lognormal.poisson = sampleMcmc(m.lognormal.poisson, thin = thin, samples = samples,
                                 transient = transient, nChains = nChains, verbose = verbose)

# Extract the posterior distribution from the model object with the convertToCodaObject function.
mpost = convertToCodaObject(m.lognormal.poisson)

# Evaluate convergence visually with trace plots.
plot(mpost$Beta)

# Our chains are still mixing fairly well, but you should now be able to notice some differences
# between the lognormal poisson model trace plots in which the chains are not fully overlapping
# one another throughout all the iterations and the trace plots of the normal and probit models
# where there was more mixing between the two chains.

# Evaluate convergence quantitatively by examining effective samples sizes.
effectiveSize(mpost$Beta)

# Compared to the linear model (and now the probit model), we observe that the effective sample 
# sizes are much smaller. This is still ok though and is to be expected because the normal model
# is the easiest case for MCMC sampling.

# Evaluate convergence quantitatively by examining potential scale reduction factors.
gelman.diag(mpost$Beta, multivariate = FALSE)$psrf

# We again see that for both parameters the PSRF estimates are very close to 1, so we can feel 
# confident that the two chains gave consistent results.


########################################################################
##  Step 7: Checking the Lognormal Poisson Model's explanatory power  ##
########################################################################

# There are many more ways to evaluate model fit for count data as you will see when running the
# below. We will not go into these in detail today, but for reference SR2 can be viewed as a
# pseudo-R squared, metrics with "O." in front of them refer to how well the model predicts 
# occurrences, and metrics with a "C." in front of them signify that the metric is conditional
# on presence. In other words, "C." metrics examine how well the model is able to predict
# abundance variation in sampling units where the species is present, ignoring absences.

# Compute the predicted values of the posterior distribution
preds = computePredictedValues(m.lognormal.poisson)

# Evaluate the model's explanatory power
evaluateModelFit(hM = m.lognormal.poisson, predY = preds)


####################################################
##       Step 8: Visualizing Our Results          ##
####################################################

# While we can plug in any value for x and use any of our 3 models to predict a value for y, we 
# often just want to visualize how the response value changes within the range of the 
# explanatory variable.

par(mfrow = c(1, 3)) # sets up a graphing panel with 1 row and 3 columns

# we can create 3 graphs (1 for each of our models) at once by using a for loop
for (i in 1:3) {     
  m = switch(i, m.normal, m.probit, m.lognormal.poisson)
  # constructGradient constructs a range of x values from which to predict
  Gradient = constructGradient(m, focalVariable = "x")
  # predict function generates predictions based on the constructed x values
  predY = predict(m, Gradient = Gradient, expected = TRUE) 
  plotGradient(m, Gradient, pred = predY, measure = "Y", index=1, showData = TRUE,
               main = c("Normal", "Probit", "Lognormal Poisson")[i] )
  }

# In these plots, the lines show the mean of the posterior distribution and the shaded area
# represents the 95% credible interval of the model prediction. Note that many points fall
# outside of the 95% credible interval because the shaded areas illustrate the credible 
# intervals for the model predictions rather than for the data.




##################################################
##    Part II: Hierarchical Random Effects      ##
##################################################


# From our simulated examples above, what if the 50 sites we sampled were grouped in a way that
# might introduce additional variation to our response not captured by our predictor variable?
# For example, imagine that half of our sites are in one country while half are in a different
# country. The country in which the site is found might affect
# our observed response (e.g. whether a species is present or absent, how many individuals
# we observe) if there are different land management strategies or other differences we did 
# not include as predictors.
# Even if our sites were all in the same geopolitical area, there may be attributes of 
# individual sites that introduce variation in our response. This is where random effects come 
# in.
# While there are much more statistically rigorous definitions, you can think of random effects 
# as parameters that could introduce variation but that you did not explicitly measure. 
# In contract, a fixed effect is a parameter you explicitly measure (i.e. we have one measure
# of x, our predictor variable, for each of our 50 sampling units. Through measuring x in 
# all of our sampling units, we can quantify the variation that x produces in the response y.)
# We have not however conducted any measurements that account for the site's effect on the 
# response y. So we could include site in the model as a random effect.


#################################
##   Step 1: Simulating Data   ##
#################################

# To illustrate this, let us take a new example in which n=120 sampling units have been
# sampled in 12 plots (~ 10 sampling units per plot). We are going to assume that plots 
# have an additive effect to the response variable and simulate the data accordingly.

# Specify the number of sampling units
n = 120 

# Generate random values of the predictor variable x, one value for each of our
# sampling units.  
x = rnorm(n)

# Define the intercept
beta1 = 0

# Define the slope
beta2 = 1

# Define the standard deviation of our observations
sigma = 1

# Define the number of plots
np = 12

# Define the standard deviation of our plot-level effect
sigma.plot = 1

# Define the general form of our linear predictor
L = beta1 + beta2*x

# Generate a vector of plot ID's (length of 120, which we will use to assign
# each of our observations to a particular plot) 
plot.id = sample(1:np, n, replace = TRUE)

# We will simulate the random effect of each plot, by generating a 12 random values
ap = rnorm(np, sd = sigma.plot)

# Create a vector of random effects (length 120) by matching the value of the effect (ap) 
# to the plot.id
a = ap[plot.id]

# Generate the response values by adding the linear predictor value, the random effect, and a
# term for the residual error
y = L + a + rnorm(n, sd = sigma)


# To see the plot-level effects, we can graph the data with a best-fit line for each plot.
par(mfrow = c(1,1))
cols = rainbow(np)
plot(x, y, col = cols[plot.id], las = 1)
for (p in 1:np){
  abline(beta1+ap[p], beta2, col = cols[p])
}


####################################################
##      Step 2: Fitting the model with Hmsc       ##
####################################################

# Designate plot ID as a factor 
plot.id = as.factor(plot.id)

# Create a sample ID for each sampling unit
sample.id = as.factor(1:n)

# Convert vector of predictor values to a dataframe (again it will have only 1 column)
XData = data.frame(x = x)

# Convert response values to a matrix (again it will only have 1 column)
Y = as.matrix(y)

# Create a dataframe that shows how sampling units and plot IDs correspond to one another
studyDesign = data.frame(sample = sample.id, plot = plot.id)

# Specify the random effect
rL = HmscRandomLevel(units = studyDesign$plot)

# Build the Hmsc model object. Now, we include studyDesign, which defines the nature of the
# study design, and ranLevels, which includes a list of random effects to be included in
# the model. (In this case, the random effect is the plot, which we created above with the
# function HmscRandomLevel.)
m = Hmsc(Y = Y, XData = XData, XFormula = ~x, studyDesign = studyDesign,
         ranLevels = list("plot" = rL))

# Fit the model by conducting the MCMC sampling
# We will again leave our inputs for the MCMC sampling (nchains, samples, etc.) the same 
# as our earlier models.
# Note that even though we have simulated data for a linear (normal) model (i.e. the
# simplest model to fit, including a random effect increases the computation time.
m = sampleMcmc(m, thin = thin, samples = samples, transient = transient, nChains = nChains,
               verbose = verbose)


# Did our model converge? Can we trust our model results?

# Extract the posterior distribution from the model object with the convertToCodaObject function.
mpost = convertToCodaObject(m)

# Evaluate convergence visually with trace plots.
plot(mpost$Beta)

# Evaluate convergence quantitatively by examining effective samples sizes.
effectiveSize(mpost$Beta)

# Now that we are back to a linear (normal) model, the actual sample size
# is again ~2000 (1000 per chain).

# Evaluate convergence quantitatively by examining potential scale reduction factors.
gelman.diag(mpost$Beta, multivariate = FALSE)$psrf

# We again see that for both parameters the PSRF estimates are very close to 1, so we can feel 
# confident that the two chains gave consistent results.


# How well does the Hmsc model fit our simulated data? 
# We can compute the R squared value to find out.

# Compute the predicted values of the posterior distribution
preds = computePredictedValues(m)

# Calculate the model fit
MF = evaluateModelFit(hM = m, predY = preds)

# View the R squared value
MF$R2


####################################################
##   Step 3: Evaluating Predictive Performance    ##
####################################################

# We know that our model converged, but how good is it at prediction?

# To determine predictive accuracy, we use a method called k-fold cross-validation.
# The idea is that we will first divide our data into a certain number of groups, so imagine
# 2 groups, #1 and #2. We will first make predictions for group #1, so we will only use group
# #2 to build the model, and then assess how accurate the predictions are. Then we repeat the
# process, but this time only building the model with the data from group #1. After this, 
# the predictions are combined to provide one matrix of predictions for all sampling units.

# The more groups (folds) you have, the more times you will have to fit the model (i.e. 
# k-folds = k times the model has to be fit). So to minimize computation time, we will stick
# with creating 2 folds only.

# Partition the data into 2 groups.
# The createPartition function randomly does this.
partition = createPartition(m, nfolds = 2, column = "sample")

# Generate model predictions according to this partition
preds = computePredictedValues(m, partition = partition)

# We have reduced computation time by limiting our evaluation to only 2 folds; however, as you
# increase the number of folds, you increase the predictive power of the model. This is because
# more data is then available to fit the model. 
# With two folds, and 120 sampling units, we fit the model using 60 data points, then predict
#   for the remaining 60.
# If we had 10 folds, we would use 108 data points to fit the model, and then predict for the
#   remaining 12.
# Consequently, when we evaluate the predictive R squared with only 2 folds, it is likely
# to be less than the explanatory R squared value we observed in Step 2.

# Calculate the model fit
MF = evaluateModelFit(hM = m, predY = preds)

# View the R squared value
MF$R2



##################################################
##      Part III: Spatial Random Effects        ##
##################################################

# Often when we are investigating ecological communities and/or metacommunities, our data 
# varies in space. In our modeling, we therefore need to account for spatial autocorrelation,
# which is the tendency for areas or sites that are close together to have similar
# values. For example, temperature is a parameter that is highly autocorrelated. 



#################################
##   Step 1: Simulating Data   ##
#################################

# To illustrate this, let us take a new example in which we have n=100 sampling units.
# We are going to assume that the amount of spatial autocorrelation decreases exponentially
# with increasing distance between points. We are also going to assume that there is a strong
# spatial effect and simulate the data such that the spatial effect is 4 times higher than the
# residual variation.

# Specify the number of sampling units
n = 100 

# Generate random values of the predictor variable x, one value for each of our
# sampling units.  
x = rnorm(n)

# Define the intercept
beta1 = 0

# Define the slope
beta2 = 1

# Define the standard deviation of our observations
sigma = 1

# Define the standard deviation of spatial effect
sigma.spatial = 2

# Define the spatial scale over which the autocorrelation decreases
alpha.spatial = 0.5

# Define the general form of our linear predictor
L = beta1 + beta2*x

# Define the x and y coordinates
# The runif function generates random values of the uniform distribution (i.e. a 
# rectangular distribution in which all values are equally likely) with a minimum value of
# 0 and a maximum value of 1. Here we are telling R to generate 200 such values and divide
# them into a matrix with two columns.
xycoords = matrix(runif(2*n), ncol = 2)

# Visually examine our simulated spatial coordinates
plot(xycoords)

# Define the exponentially decreasing spatial autocorrelation 
Sigma = sigma.spatial^2 * exp(-as.matrix(dist(xycoords)) / alpha.spatial)

# What we have just done with the above line is create a 100x100 matrix, which specifies
# the degree of autocorrelation between each spatial coordinate and all the other spatial
# coordinates. 

Sigma[1:5, 1:5]

# Examine the first five rows and columns. Note that the diagonal (i.e. value for 1,1; 2,2;
# 3,3; 4,4; 5,5; etc.) is always 4. This is because we set the variance of the spatial random
# effect as 4 above when we set the standard deviation (sigma.spatial) = 2 (Recall that
# standard deviation squared = variance. 

# Generate the random spatial effect for each sampling unit
a = mvrnorm(mu = rep(0, n), Sigma = Sigma)

# Generate the response values by adding the linear predictor value, the random spatial effect,
# and a term for the residual error
y = L + a + rnorm(n, sd = sigma)



####################################################
##      Step 2: Fitting the model with Hmsc       ##
####################################################

# Create a sample ID for each sampling unit
sample.id = as.factor(1:n)

# Create a dataframe of the study design which this time will only contain sampling ID
studyDesign = data.frame(sample = sample.id)

# Rename the rows of the xycoords matrix according to sample ID
rownames(xycoords) = sample.id

# Specify the random spatial effect
rL = HmscRandomLevel(sData = xycoords)

# Convert vector of predictor values to a dataframe (again it will have only 1 column)
XData = data.frame(x = x)

# Convert response values to a matrix (again it will only have 1 column)
Y = as.matrix(y)

# Build the Hmsc model object. 
m = Hmsc(Y = Y, XData = XData, XFormula = ~x, studyDesign = studyDesign,
         ranLevels = list("sample" = rL))

# Fit the model by conducting the MCMC sampling
# We will again leave our inputs for the MCMC sampling (nchains, samples, etc.) the same 
# as our earlier models.
# Note again that even though we have simulated data for a linear (normal) model (i.e. the
# simplest model to fit, including a random effect increases the computation time.
m = sampleMcmc(m, thin = thin, samples = samples, transient = transient, nChains = nChains,
               verbose = verbose)


# Did our model converge? Can we trust our model results?

# Extract the posterior distribution from the model object with the convertToCodaObject function.
mpost = convertToCodaObject(m)

# Evaluate convergence visually with trace plots.
plot(mpost$Beta)

# Evaluate convergence quantitatively by examining effective samples sizes.
effectiveSize(mpost$Beta)

# Now that we are back to a linear (normal) model, the actual sample size
# is again ~2000 (1000 per chain).

# Evaluate convergence quantitatively by examining potential scale reduction factors.
gelman.diag(mpost$Beta, multivariate = FALSE)$psrf

# We again see that for both parameters the PSRF estimates are very close to 1, so we can feel
# confident that the two chains gave consistent results.


# How well does the Hmsc model fit our simulated data? 
# We can compute the R squared value to find out.

# Compute the predicted values of the posterior distribution
preds = computePredictedValues(m)

# Calculate the model fit
MF = evaluateModelFit(hM = m, predY = preds)

# View the R squared value
MF$R2


####################################################
##   Step 3: Evaluating Predictive Performance    ##
####################################################

# We will now evaluate the predictive power of the model through two-fold cross validation.

# Partition the data into 2 groups.
partition = createPartition(m, nfolds = 2, column = "sample")

# Generate model predictions according to this partition
preds = computePredictedValues(m, partition = partition)

# Calculate the model fit
MF = evaluateModelFit(hM = m, predY = preds)

# View the R squared value
MF$R2

# Given what you know about spatial autocorrelation, why might using a spatial model improve
# predictive accuracy?



##################################################
##      Part IV: A Read World Example           ##
##################################################

# We will now move beyond simulated data to a real world example using data on Finnish birds.
# Again though, we will conduct a univariate analysis (so investigating the distribution of
# only 1 species).

# Before we begin we clear our workspace and memory to improve efficiency.
rm(list = ls() ) 
gc() #releases memory

# Download the bird data (data.csv) and the grid data (grid_1000.csv) from the class github page if you haven't already.

# Set your working directory to where you stored the folder.
setwd()

# Load the bird data
da_all = read.csv("data.csv", header = TRUE)

# We are going to analyze data for Corvus monedula in 2014, but you can look at what other
# species, years, and predictors are present in the data frame.
names(da_all)
unique(da_all$Year)

# You can see that we have some identifying attributes of the sampling units (Route, Year,
# x coord, and y coord), survey effort (helpful for quantifying detection probability),
# and some environmental predictors (Habitat, and different measures of temperature). Then
# we have data on many different bird species.

# Let's look at the values for some of these species to see if we are dealing with presence
# absence or count data.
da_all[1:10, 10:13]

# We will examine Corvus monedula during 2014, so first we subset our data frame to show us
# only this year.
da = droplevels(subset(da_all, Year==2014))

# We will examine the relationship between distribution and Habitat and April/May temperature.
# So we create our data frame of predictor variables:
XData = data.frame(hab = da$Habitat, clim = da$AprMay)

# We also need to make sure that R sees our habitat variable as a factor.
XData$hab = as.factor(XData$hab)

# Examining our XData shows us that the habitat variable is categorical and is comprised of
# 5 different codes.
unique(XData$hab)

# These correspond to broadleaved forest, coniferous forest, open habitat, urban habitat, 
# and wetlands.

# April/May clim refers to mean spring temperature in April and May.
hist(XData$clim)

# Create the matrix for the response using only the Corvus monedula data. Again, our matrix
# will only have 1 column.
Y = as.matrix(da$Corvus_monedula)

# We additionally need to create a matrix of our x and y spatial coordinates.
xy = as.matrix(cbind(da$x, da$y))

# We can plot our environmental predictors and our bird data to get a sense of how each vary
# throughout space. To do this we first have to create a spatial points dataframe.
spdf = SpatialPointsDataFrame(coords = xy, data = da[, c(6, 9, 59)])

par(mfrow = c(1,3))
plot(spdf, pch = spdf$Habitat, main = "Habitat")
plot(spdf, pch = 21, col = spdf$AprMay, main = "Spring Temp")
plot(spdf, pch = 21, bg = spdf$Corvus_monedula, main = "C. monedula Abundance")



#########################################
##   Setting up and fitting the model  ##
#########################################

# Since we have count data we will use a lognormal Poisson model. As fixed effects, we will
# include habitat and temperature. Since sometimes temperature has a non-linear effect, we will
# include both the linear and squared (i.e. quadratic) values of temperature. Including the
# squared effect of temperature allows species' abundance to peak at intermediate climate
# conditions.

# We will also include a spatial random effect to account for the spatial nature of the 
# study design.

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

# Build the model object
mFULL = Hmsc(Y = Y, XData = XData, XFormula = XFormula, distr = "lognormal poisson",
             studyDesign = studyDesign, ranLevels = list(route = rL))

# Specify the parameters for model fitting
nChains = 4
samples = 250
thin = 10
transient = round(0.5*samples*thin)
verbose = round(0.5*samples*thin)

# Model fitting with MCMC
# We are adding an additional argument here with initPar. While there is a more techincal
# explanation of what this term means, overall we are using it here to help ensure
# convergence.
mFULL = sampleMcmc(mFULL, thin = thin, samples = samples, transient = transient, 
                   nChains = nChains, verbose = verbose, initPar = "fixed effects")


# Did our model converge? Can we trust our model results?
# Evaluate MCMC convergence for the beta parameters.

# Extract the posterior distribution from the model object with the convertToCodaObject function.
# Here we use spNamesNumbers to refer to the species by their names instead of by their numbers.
# Similarly, covNamesNumbers is used to refer to the explanatory variables by their names
# instead of their numbers.
mpost = convertToCodaObject(mFULL, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))

# Evaluate convergence quantitatively by examining effective samples sizes.
effectiveSize(mpost$Beta)

# We had 137 data points of abundances but some of these abundance values were zero, especially
# in th northern part of the country where Wetlands predominate. We can see the lower effective
# sample size reflected in this for Wetlands. 

# Evaluate convergence quantitatively by examining potential scale reduction factors.
gelman.diag(mpost$Beta, multivariate = FALSE)$psrf



# We will also evaluate the predictive power of our model through two-fold cross validation.
# Partition the data into 2 groups.
partition = createPartition(mFULL, nfolds = 2, column = "route")

# Generate model predictions according to this partition
preds = computePredictedValues(mFULL, partition = partition)

# Calculate the model fit
MF = evaluateModelFit(hM = mFULL, predY = preds)

# View the fit metrics
MF

# What is your AUC value? What does this tell you about how well the model will predict
# C. monedula distribution?



#########################################################
##  Visualizing the results: Environmental Gradients   ##
#########################################################

# We can look at parameter estimates in a table; however, these results are often difficult to 
# interpret. To understand what these parameter estimates mean, we will construct prediction
# plots.
round(summary(mpost$Beta, quantiles = c(0.025, 0.5, 0.975))
      [[2]],2)


#### C. monedula and Temperature

par(mfrow = c(1,2))

# Gradient constructs a range of x values from which to predict for our first plot
# Setting the non-focal variable to hab=1, will make these predictions for the most 
# common habitat type in the data, which is coniferous forest.
Gradient = constructGradient(mFULL, focalVariable = "clim",
                             non.focalVariables = list(hab = 1))

# Generate predictions based on the constructed x values for our first plot
predY = predict(mFULL, Gradient = Gradient, expected = TRUE) 

# Create the first plot
plotGradient(mFULL, Gradient, pred = predY, measure = "Y", index=1, showData = TRUE)

# Gradient constructs a range of x values from which to predict for our second plot
# Setting the non-focal variable to hab=2, means that the values of habitat are modelled
# as functions of climate, using generalized linear models.
# In other words, when moving along the climatic gradient (from coldest to warmest), 
# we assume that under coldest conditions that the species survey took place in a wetland,
# whereas under the other climatic conditions it has been conducted in a coniferous forest.
Gradient = constructGradient(mFULL, focalVariable = "clim",
                             non.focalVariables = list(hab = 2))

# Generate predictions based on the constructed x values for our second plot
predY = predict(mFULL, Gradient = Gradient, expected = TRUE) 

# Create the second plot
plotGradient(mFULL, Gradient, pred = predY, measure = "Y", index=1, showData = TRUE)


#### C. monedula and Habitat

# Note the plotting procedure is a little different for this type of plot

# Construct the range of x values from which to predict 
Gradient = constructGradient(mFULL, focalVariable = "hab",
                             non.focalVariables = list(clim = 1))

# Generate predictions based on the constructed x values for our first plot
predY = predict(mFULL, Gradient = Gradient, expected = TRUE) 

# Create the first plot
a <- plotGradient(mFULL, Gradient, pred = predY, measure = "Y", index=1, showData = TRUE)

# Construct the range of x values from which to predict 
Gradient = constructGradient(mFULL, focalVariable = "hab",
                             non.focalVariables = list(clim = 2))

# Generate predictions based on the constructed x values for our second plot
predY = predict(mFULL, Gradient = Gradient, expected = TRUE) 

# Create the second plot
b <- plotGradient(mFULL, Gradient, pred = predY, measure = "Y", index=1, showData = TRUE)

# Arrange the two plots together
gridExtra :: grid.arrange(a, b, ncol = 2)


###########################################################################
##  Visualizing the results: Distribution of C. monedula across Finland  ##
###########################################################################

# To predict the distribution of C. monedula across Finland, we will use a grid
# of 10000 Finnish locations, with information about the habitat type, climatic
# conditions, and spatial coordinates for each location.

# Read in the grid dataframe
grid = read.csv("grid_1000.csv")

# Subset the dataframe to exclude the marine habitat category. While there are locations
# in the dataframe for this habitat type, the data we trained our model on did not 
# include it, and so our model will not be able to make predictions for it.
grid = droplevels(subset(grid, !(Habitat=="Ma")))

# Create a matrix of x and y coordinates
xy.grid = as.matrix(cbind(grid$x, grid$y))

# Create a dataframe of our predictor variables
XData.grid = data.frame(hab = grid$Habitat, clim = grid$AprMay)

# Construct the range of x values from which to predict 
Gradient = prepareGradient(mFULL, XDataNew = XData.grid, sDataNew = list(route = xy.grid))

# Generate predictions based on the constructed x values
predY = predict(mFULL, Gradient = Gradient)

# Summarize the posterior predicted distribution by its mean value (i.e. expected count
# of C. monedula at each location).
EpredY = apply(abind(predY, along = 3), c(1,2), mean)


# Mapping the results
# Create a dataframe for the data we will need to plot and build 3 separate graphs
# to show the variation in habitat type, April/May temperature, and the predicted
# abundance of C. monedula

# Create the dataframe with the data we'll need for plotting
mapData=data.frame(x=xy.grid[,1],y=xy.grid[,2], EpredY,H=XData.grid$hab, C=XData.grid$clim)

# Habitat Map
a <- ggplot(data = mapData, aes(x=x, y=y, color=H))+geom_point(size=1.5) +
  ggtitle("Habitat") + scale_color_discrete(name="") +
  theme(legend.position="bottom")

# Climate Map
b <- ggplot(data = mapData, aes(x=x, y=y, color=C))+geom_point(size=1.5) +
  ggtitle("Climate") + scale_color_gradient(low="blue", high="red", name = "") +
  theme(legend.position="bottom") +labs(y="")

# Abundance Map
c <- ggplot(data = mapData, aes(x=x, y=y, color=sp1))+geom_point(size=1.5) + 
  ggtitle("Corvus monedula")+ scale_color_gradient(low="blue", high="red", name = "") +
  theme(legend.position="bottom")+labs(y="")

# Arrange the three plots together
gridExtra :: grid.arrange(a, b, c, ncol = 3)


