# BayesRecipe (Bayesian Reciprocal Regularization)

## Introduction

This R package implements Bayesian reciprocal regularization and variable selection for regression and classification. Currently, it includes a set of computationally efficient MCMC algorithms (Gibbs and slice samplers) for solving the Bayesian reciprocal LASSO in linear models. It also includes a modified S5 algorithm to solve the (frequentist) reduced reciprocal LASSO problem in linear regression. For technical details, see our papers referenced below.


## Dependencies

`BayesRecipe` requires the following `R` package: `devtools` (for installation only). Please install it before installing `BayesRecipe`, which can be done as follows (execute from within a fresh R session):
```r
install.packages("devtools")
library(devtools)
```

## Installation

Once the dependencies are installed, `BayesRecipe` can be loaded using the following command:
```r
devtools::install_github("himelmallick/BayesRecipe")
library(BayesRecipe)
```

In addition, to run the slice sampler, we need to install and load the following:
```r
devtools::install_github("JingyuHe/bayeslm")
library(bayeslm)
```


## Basic Usage
```r
BayesRLasso(x, y)
```


## Input

- **x**: A numeric matrix with standardized predictors in columns and samples in rows.
- **y**: A mean-centered continuous response variable with matching rows with `x`.
- **method**: If Gibbs sampling is desired, 'SMU' or 'SMN' must be selected as the data augmentation scheme or mixture representation, although not recommended when the data size is as big as hundreds of covariates. By default it uses an ellipitical slice sampler ('slice') which is scalable for large-scale problems.
- **lambda.estimate**: Estimating lambda by empirical bayes ('EB'), MCMC ('MCMC'), or apriori ('AP'). Default is 'AP' (currently only available option when method = 'slice').
- **update.sigma2**: Whether sigma2 should be updated. Default is TRUE (currently only available option when method = 'slice').
- **max.steps**: Number of MCMC iterations. Default is 11000.
- **n.burn**: Number of burn-in iterations. Default is 1000.
- **n.thin**: Lag at which thinning should be done. Default is 1 (no thinning).
- **ridge.CV**: When method = 'SMU' and X is rank deficient, a ridge parameter is added to the diagonal of the crossproduct (XtX) to allow for proper calculation of the inverse. If TRUE, the ridge parameter is estimated by cross-validation using glmnet. Otherwise, it falls back to adding a small number (1e-05) without the cross-validation. Default is TRUE. This parameter is ignored when method = 'SMN' or 'slice'.
- **a**: If lambda.estimate = 'MCMC', shape hyperparameter for the Gamma prior on lambda. Default is 0.001. This parameter is ignored when method = 'slice'.
- **b**: IIf lambda.estimate = 'MCMC', rate hyperparameter for the Gamma prior on lambda. Default is 0.001. This parameter is ignored when method = 'slice'.
- **posterior.summary.beta**: Posterior summary measure for beta (mean, median, or mode). Default is 'mean'. 
- **posterior.summary.lambda**: Posterior summary measure for lambda (mean, median, or mode). Default is 'median'. 
- **beta.ci.level**: Credible interval level for beta. Default is 0.95 (95%).
- **lambda.ci.level**: Credible interval level for lambda. Default is 0.95 (95%).
- **seed**: Seed value for reproducibility. Default is 1234.



## Output

A list containing the following components is returned:

- **time**: Computational time in minutes.
- **beta**: Posterior estimates of beta.
- **lowerbeta**: Lower limit of the credible interval of beta.
- **upperbeta**: Upper limit of the credible interval of beta.
- **lambda**: Posterior estimate of lambda.
- **lambdaci**: Posterior credible interval of lambda.
- **beta.post**: Post-burn-in posterior samples of beta.
- **sigma2.post**: Post-burn-in posterior samples of sigma2.
- **lambda.post**: Post-burn-in posterior samples of lambda.


## Illustration of Reciprocal Bayesian LASSO

Let's consider the `prostate` dataset from the `ElemStatLearn` package.

```r

#########################
# Load Prostate dataset #
#########################

library(ElemStatLearn)
prost<-ElemStatLearn::prostate

###########################################
# Scale data and prepare train/test split #
###########################################

prost.std <- data.frame(cbind(scale(prost[,1:8]),prost$lpsa))
names(prost.std)[9] <- 'lpsa'
data.train <- prost.std[prost$train,]
data.test <- prost.std[!prost$train,]

```

First, we standardize the variables and center the response variable. We also separate out training and test samples using a train/test split included in the dataset.

```r

##################################
# Extract standardized variables #
##################################

y.train   = data.train$lpsa - mean(data.train$lpsa)
y.test <- data.test$lpsa - mean(data.test$lpsa)
x.train = scale(as.matrix(data.train[,1:8], ncol=8))
x.test = scale(as.matrix(data.test[,1:8], ncol=8))


```

Let's apply reciprocal Bayesian LASSO and calculate predictions in the test set.

```r

#############################
# Reciprocal Bayesian Lasso #
#############################

library(BayesRecipe) 
fit_BayesRLasso<-BayesRLasso(x.train, y.train)
y.pred.BayesRLasso<-x.test%*%fit_BayesRLasso$beta

```

Let's compare with various frequentist methods such as OLS, LASSO, adaptive LASSO, classical bridge, elastic net, SCAD, MCP, and the frequentist reciprocal LASSO (rLASSO).

```r 

####################
# Compare with OLS #
####################

ols.lm<-lm(y.train ~ x.train)
y.pred.ols<-x.test%*%ols.lm$coefficients[-1]

##################################
# Compare with Frequentist LASSO #
##################################

set.seed(1234)
library(glmnet)
lasso.cv=cv.glmnet(x.train,y.train,alpha = 1)  
lambda.cv.lasso=lasso.cv$lambda.min          
lasso.sol=glmnet(x.train,y.train,alpha = 1)    
lasso.coeff=predict(lasso.sol,type="coef",s=lambda.cv.lasso)
y.pred.lasso=x.test%*%lasso.coeff[-1]

###########################################
# Compare with Frequentist Adaptive LASSO #
###########################################

alasso.cv=cv.glmnet(x.train,y.train,alpha = 1,penalty.factor=1/abs(ols.lm$coefficients))  
lambda.cv.alasso=alasso.cv$lambda.min          
alasso.sol=glmnet(x.train,y.train,alpha = 1,penalty.factor=1/abs(ols.lm$coefficients))    
alasso.coeff=predict(alasso.sol,type="coef",s=lambda.cv.alasso)
y.pred.alasso=x.test%*%alasso.coeff[-1]

###################################
# Compare with Frequentist Bridge #
###################################

library(grpreg)
group<-c(rep(1,8))
fit.bridge <- gBridge(x.train, y.train,group=group)
bridge.coeff<-grpreg::select(fit.bridge,'AIC')$beta
y.pred.bridge<-x.test%*%bridge.coeff[-1]

########################################
# Compare with Frequentist Elastic Net #
########################################

library(glmnet)
enet.cv=cv.glmnet(x.train,y.train,alpha = 0.5)  
lambda.cv.enet=enet.cv$lambda.min          
enet.sol=glmnet(x.train,y.train,alpha = 0.5)    
enet.coeff=predict(enet.sol,type="coef", s=lambda.cv.enet)
y.pred.enet=x.test%*%enet.coeff[-1]

#################################
# Compare with Frequentist SCAD #
#################################

library(ncvreg)
scad.cv <- cv.ncvreg(x.train,y.train,penalty="SCAD", dfmax=1000,max.iter=10^4)
lambda.cv.scad<- scad.cv$lambda.min
scad.coeff<- ncvreg(X=x.train,y=y.train, penalty='SCAD',dfmax=1000,lambda=lambda.cv.scad)$beta[-1,1]
y.pred.scad = x.test%*%scad.coeff

#################################
# Compare with Frequentist MCP #
#################################

mcp.cv <- cv.ncvreg(x.train,y.train,penalty="MCP", dfmax=1000,max.iter=10^4)
lambda.cv.mcp<- mcp.cv$lambda.min
mcp.coeff<- ncvreg(X=x.train,y=y.train, penalty='MCP',dfmax=1000,lambda=lambda.cv.mcp)$beta[-1,1]
y.pred.mcp = x.test%*%mcp.coeff

###################################
# Compare with Frequentist rLASSO #
###################################

fit_rLASSO = rrLASSO.S5(x.train,y.train)
rlasso.coeff<-fit_rLASSO$beta
y.pred.rlasso<-x.test%*%rlasso.coeff[-1]

```

Let's now compare the performance in test data, which reveals that the reciprocal Bayesian LASSO attains the minimum test error in this dataset.

```r 

###########################################
# Prediction Accuracy (MSE on Test Data ) #
###########################################

mean((y.pred.ols-y.test)^2) # 0.5421042 OLS
mean((y.pred.lasso-y.test)^2) # 0.5150523 LASSO
mean((y.pred.alasso-y.test)^2) # 0.7114202 ALASSO
mean((y.pred.bridge-y.test)^2) # 0.5168394 BRIDGE
mean((y.pred.enet-y.test)^2) # 0.506828 ENET
mean((y.pred.scad-y.test)^2) # 0.5388646 SCAD
mean((y.pred.mcp-y.test)^2) # 0.5388779 MCP
mean((y.pred.rlasso - y.test)^2) # 0.5413048 rLASSO
mean((y.pred.BayesRLasso - y.test)^2) # 0.4868719 BayesRLasso
```

Let's check MCMC convergence of the reciprocal Bayesian LASSO estimator through two visualizations: trace plots and histograms.


```r 

######################################
# Visualization of Posterior Samples #
######################################

##############
# Trace Plot #
##############

library(coda)
par(mar=c(2,2,2,2))
plot(mcmc(fit_BayesRLasso$beta.post),density=FALSE,smooth=TRUE)

```

![plot of chunk traceplot](https://github.com/himelmallick/BayesRecipe/blob/master/misc/traceplot_prostate.png)

```r 

#############
# Histogram #
#############

library(psych)
multi.hist(fit_BayesRLasso$beta.post,density=TRUE,main="")

```

![plot of chunk histogram](https://github.com/himelmallick/BayesRecipe/blob/master/misc/histogram_prostate.png)

It is satisfactory to observe that for this benchmark dataset, the samples traverse the posterior space very fast confirming that the sampler has good mixing properties.

## Issues and Next Steps

As above, note that this tutorial is still under development. Meanwhile, we are happy to troubleshoot any issue with the package. Please reach out directly via email or create an issue in the repository.

## Citation

If you use `BayesRecipe` in your work, please cite the following as appropriate:

1. Mallick H, Alhamzawi R, Paul E, Svetnik V (2021). [The Reciprocal Bayesian LASSO](https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.9098). *Statistics in Medicine*, 40(22):4830-4849.

2. Alhamzawi R, Mallick H (2022). [Bayesian Reciprocal LASSO Quantile Regression](https://www.tandfonline.com/doi/abs/10.1080/03610918.2020.1804585). *Communications in Statistics - Simulation and Computation*, 51(11): 6479-6494.
  
3. Alhamzawi R, Alhamzawi A, Mallick H (2023). [The Reciprocal Elastic Net](https://www.tandfonline.com/doi/abs/10.1080/23737484.2023.2278106). *Communications in Statistics – Case Studies and Data Analysis*, 9(4), 422–436.

4. Paul E, Mallick H (2024+). [Unified Reciprocal LASSO Estimation Via Least Squares Approximation](https://www.tandfonline.com/doi/abs/10.1080/03610918.2022.2146723). *Communications in Statistics - Simulation and Computation*, In Press.

5. Paul E, He J, Mallick H (2024+). [Accelerated Bayesian Reciprocal LASSO](https://www.tandfonline.com/doi/abs/10.1080/03610918.2023.2276050). *Communications in Statistics - Simulation and Computation*, In Press.




