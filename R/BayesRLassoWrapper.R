#' Bayesian Reciprocal Regularization
#'
#' Wrapper for Bayesian Reciprocal LASSO using MCMC
#' 
#' @param x A numeric matrix with standardized predictors in columns and samples in rows.
#' @param y A mean-centered continuous response variable with matching rows with x.
#' @param method If Gibbs sampling is desired, 'SMU' or 'SMN' must be selected as the data augmentation scheme 
#' or mixture representation, although not recommended when the data size is as big as hundreds of covariates.
#' By default it uses an ellipitical slice sampler ('slice') which is scalable for large-scale problems.
#' @param update.sigma2 Whether sigma2 should be updated. Default is TRUE 
#' (currently only available option when method = 'slice').
#' @param lambda.estimate Estimating lambda by empirical bayes ('EB'), MCMC ('MCMC'), or apriori ('AP'). 
#' Default is 'AP' (currently only available option when method = 'slice').
#' @param max.steps Number of MCMC iterations. Default is 11000.
#' @param n.burn Number of burn-in iterations. Default is 1000.
#' @param n.thin Lag at which thinning should be done. Default is 1 (no thinning).
#' @param ridge.CV When method = 'SMU' and X is rank deficient, a ridge parameter is added to the diagonal of the 
#' crossproduct (XtX) to allow for proper calculation of the inverse. If TRUE, the ridge parameter
#' is estimated by cross-validation using glmnet. Otherwise, it falls back to adding a small number (1e-05) 
#' without the cross-validation. Default is TRUE. This parameter is ignored when method = 'SMN' or 'slice'. 
#' @param a If lambda.estimate = 'MCMC', shape hyperparameter for the Gamma prior on lambda. Default is 0.001.
#' This parameter is ignored when method = 'slice'. 
#' @param b If lambda.estimate = 'MCMC', rate hyperparameter for the Gamma prior on lambda. Default is 0.001.
#' This parameter is ignored when method = 'slice'. 
#' @param posterior.summary.beta Posterior summary measure for beta (mean, median, or mode). Default is 'mean'. 
#' @param posterior.summary.lambda Posterior summary measure for lambda (mean, median, or mode). Default is 'median'. 
#' @param beta.ci.level Credible interval level for beta. Default is 0.95 (95\%).
#' @param lambda.ci.level Credible interval level for lambda. Default is 0.95 (95\%).
#' @param seed Seed value for reproducibility. Default is 1234.
#' @param na.rm Logical. Should missing values (including NaN) be omitted from the calculations?
#' 
#' @return A list containing the following components is returned:
#' \item{time}{Computational time in minutes.}
#' \item{beta}{Posterior estimates of beta.}
#' \item{lowerbeta}{Lower limit of the credible interval of beta.}
#' \item{upperbeta}{Upper limit of the credible interval of beta.}
#' \item{lambda}{Posterior estimate of lambda.}
#' \item{lambdaci}{Posterior credible interval of lambda.}
#' \item{beta.post}{Post-burn-in posterior samples of beta.}
#' \item{sigma2.post}{Post-burn-in posterior samples of sigma2.}
#' \item{lambda.post}{Post-burn-in posterior samples of lambda.}
#' 
#' @examples
#' \dontrun{
#' 
#' #########################
#' # Load Prostate dataset #
#' #########################
#' 
#' library(ElemStatLearn)
#' prost<-prostate
#' 
#' ###########################################
#' # Scale data and prepare train/test split #
#' ###########################################
#' 
#' prost.std <- data.frame(cbind(scale(prost[,1:8]),prost$lpsa))
#' names(prost.std)[9] <- 'lpsa'
#' data.train <- prost.std[prost$train,]
#' data.test <- prost.std[!prost$train,]
#'
#' ##################################
#' # Extract standardized variables #
#' ##################################
#' 
#' y.train   = data.train$lpsa - mean(data.train$lpsa)
#' y.test <- data.test$lpsa - mean(data.test$lpsa)
#' x.train = scale(as.matrix(data.train[,1:8], ncol=8))
#' x.test = scale(as.matrix(data.test[,1:8], ncol=8))
#' 
#' #############################
#' # Reciprocal Bayesian LASSO #
#' #############################
#' 
#' fit_BayesRLasso<- BayesRLasso(x.train, y.train)
#' y.pred.BayesRLasso<-x.test%*%fit_BayesRLasso$beta
#' mean((y.pred.BayesRLasso - y.test)^2) # Performance on test data
#'
#' ######################################
#' # Visualization of Posterior Samples #
#' ######################################
#' 
#' ##############
#' # Trace Plot #
#' ##############
#' 
#' library(coda)
#' plot(mcmc(fit_BayesRLasso$beta.post),density=FALSE,smooth=TRUE)
#' 
#' #############
#' # Histogram #
#' #############
#' 
#' library(psych)
#' multi.hist(fit_BayesRLasso$beta.post,density=TRUE,main="")
#' 
#' }
#' 
#' @export
#' @keywords SMU, SMN, Ellipitical slice sampler, MCMC, Bayesian regularization, Reciprocal LASSO

BayesRLasso<-function(x,
                      y, 
                      method = 'slice',
                      lambda.estimate = 'AP',
                      update.sigma2 = TRUE, 
                      max.steps = 11000,
                      n.burn = 1000,
                      n.thin = 1, 
                      ridge.CV = TRUE,
                      a = 0.001,
                      b = 0.001,
                      posterior.summary.beta = 'mean',
                      posterior.summary.lambda = 'median',
                      beta.ci.level = 0.95,
                      lambda.ci.level = 0.95,
                      seed = 1234,
                      na.rm = TRUE) {
  
  if (method == "slice") {
    
    # Fit using elliptical slice sampler
    fit <- tryCatch({
      fit1 <- BayesRLassoSlice(x, 
                             y, 
                             lambda.estimate = lambda.estimate,
                             update.sigma2 = update.sigma2, 
                             max.steps = max.steps,
                             n.burn = n.burn,
                             n.thin = n.thin, 
                             ridge.CV = ridge.CV,
                             a = a,
                             b = b,
                             posterior.summary.beta = posterior.summary.beta,
                             posterior.summary.lambda = posterior.summary.lambda,
                             beta.ci.level = beta.ci.level,
                             lambda.ci.level = lambda.ci.level,
                             seed = seed,
                             na.rm = na.rm)
    }, error=function(err){
      message(err) # If an error occurs, try running again 
      message('\n An error occurred during the MCMC sampling. Retrying ...\n')
      fit1 <- try({ fit1 <- BayesRLassoSlice(x, 
                                           y, 
                                           lambda.estimate = lambda.estimate,
                                           update.sigma2 = update.sigma2, 
                                           max.steps = max.steps,
                                           n.burn = n.burn,
                                           n.thin = n.thin, 
                                           ridge.CV = ridge.CV,
                                           a = a,
                                           b = b,
                                           posterior.summary.beta = posterior.summary.beta,
                                           posterior.summary.lambda = posterior.summary.lambda,
                                           beta.ci.level = beta.ci.level,
                                           lambda.ci.level = lambda.ci.level,
                                           seed = seed, 
                                           na.rm = na.rm)}) 
      return(fit1)
    })
  } else if (method == "SMU"){
  
    # Fit using SMU
    fit <- tryCatch({
      fit1 <- BayesRLassoSMU(x, 
                             y, 
                             lambda.estimate = lambda.estimate,
                             update.sigma2 = update.sigma2, 
                             max.steps = max.steps,
                             n.burn = n.burn,
                             n.thin = n.thin, 
                             ridge.CV = ridge.CV,
                             a = a,
                             b = b,
                             posterior.summary.beta = posterior.summary.beta,
                             posterior.summary.lambda = posterior.summary.lambda,
                             beta.ci.level = beta.ci.level,
                             lambda.ci.level = lambda.ci.level,
                             seed = seed,
                             na.rm = na.rm)
    }, error=function(err){
      message(err) # If an error occurs, try running again with fixed sigma2.
      message('\n An error occurred during the MCMC sampling. Retrying with a fixed sigma2...\n')
      fit1 <- try({ fit1 <- BayesRLassoSMU(x, 
                                           y, 
                                           lambda.estimate = lambda.estimate,
                                           update.sigma2 = FALSE, 
                                           max.steps = max.steps,
                                           n.burn = n.burn,
                                           n.thin = n.thin, 
                                           ridge.CV = ridge.CV,
                                           a = a,
                                           b = b,
                                           posterior.summary.beta = posterior.summary.beta,
                                           posterior.summary.lambda = posterior.summary.lambda,
                                           beta.ci.level = beta.ci.level,
                                           lambda.ci.level = lambda.ci.level,
                                           seed = seed, 
                                           na.rm = na.rm)}) 
      return(fit1)
    })
    } else if (method == "SMN"){
    
    # Fit using SMN
    fit <- tryCatch({
      fit1 <- BayesRLassoSMN(x, 
                             y, 
                             lambda.estimate = lambda.estimate,
                             update.sigma2 = update.sigma2, 
                             max.steps = max.steps,
                             n.burn = n.burn,
                             n.thin = n.thin, 
                             ridge.CV = ridge.CV,
                             a = a,
                             b = b,
                             posterior.summary.beta = posterior.summary.beta,
                             posterior.summary.lambda = posterior.summary.lambda,
                             beta.ci.level = beta.ci.level,
                             lambda.ci.level = lambda.ci.level,
                             seed = seed,
                             na.rm = na.rm)
    }, error=function(err){
      message(err) # If an error occurs, try running again with fixed sigma2.
      message('\n An error occurred during the MCMC sampling. Retrying with a fixed sigma2...\n')
      fit1 <- try({ fit1 <- BayesRLassoSMN(x, 
                                           y, 
                                           lambda.estimate = lambda.estimate,
                                           update.sigma2 = FALSE, 
                                           max.steps = max.steps,
                                           n.burn = n.burn,
                                           n.thin = n.thin, 
                                           ridge.CV = ridge.CV,
                                           a = a,
                                           b = b,
                                           posterior.summary.beta = posterior.summary.beta,
                                           posterior.summary.lambda = posterior.summary.lambda,
                                           beta.ci.level = beta.ci.level,
                                           lambda.ci.level = lambda.ci.level,
                                           seed = seed,
                                           na.rm = na.rm)}) 
      return(fit1)
    })
    
  } else {
      print("Unrecognized method.  Use \"slice\", \"SMU\" or \"SMN\".")
        }
  fit
}



