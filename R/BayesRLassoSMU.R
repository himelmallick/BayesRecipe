#' Bayesian Reciprocal Regularization
#'
#' MCMC for Bayesian Reciprocal LASSO using inverse SMU mixture 
#' 
#' @param x A numeric matrix with standardized predictors in columns and samples in rows.
#' @param y A mean-centered continuous response variable with matching rows with x.
#' @param lambda.estimate Estimating lambda by empirical bayes ('EB'), MCMC ('MCMC'), or apriori ('AP'). Default is 'AP'.
#' @param update.sigma2 Whether sigma2 should be updated. Default is TRUE.
#' @param max.steps Number of MCMC iterations. Default is 11000.
#' @param n.burn Number of burn-in iterations. Default is 1000.
#' @param n.thin Lag at which thinning should be done. Default is 1 (no thinning).
#' @param ridge.CV If X is rank deficient, a ridge parameter is added to the diagonal of the 
#' crossproduct (XtX) to allow for proper calculation of the inverse. If TRUE, the ridge parameter
#' is estimated by cross-validation using glmnet. Otherwise, it falls back to adding a small number (1e-05) 
#' without the cross-validation. Default is TRUE.
#' @param a If lambda.estimate = 'MCMC', shape hyperparameter for the Gamma prior on lambda. Default is 0.001.
#' @param b If lambda.estimate = 'MCMC', rate hyperparameter for the Gamma prior on lambda. Default is 0.001.
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
#' ###################################
#' # Reciprocal Bayesian LASSO (SMU) #
#' ###################################
#' 
#' fit_BayesRLasso_SMU<- BayesRLasso(x.train, y.train, method = 'SMU')
#' y.pred.BayesRLasso_SMU<-x.test%*%fit_BayesRLasso_SMU$beta
#' mean((y.pred.BayesRLasso_SMU - y.test)^2) # Performance on test data
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
#' plot(mcmc(fit_BayesRLasso_SMU$beta.post),density=FALSE,smooth=TRUE)
#' 
#' #############
#' # Histogram #
#' #############
#' 
#' library(psych)
#' multi.hist(fit_BayesRLasso_SMU$beta.post,density=TRUE,main="")
#' 
#' }
#' 
#' @export
#' @keywords SMU, MCMC, Bayesian regularization, Reciprocal LASSO

BayesRLassoSMU<-function(x,
                         y, 
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
  # Set random seed
  set.seed(seed)

  # Set parameter values and extract relevant quantities
	x <- as.matrix(x)
 	n <- nrow(x)
	p <- ncol(x)

	# Calculate cross-product quantities for OLS
	XtX <- t(x) %*% x	
	xy <- t(x) %*% y

	# Initialize the coefficient vectors and matrices
	betaSamples <- matrix(0, max.steps, p)
	sigma2Samples <- rep(0, max.steps)
	uSamples <- matrix(0, max.steps, p)
	lambdaSamples <- rep(0, max.steps)

	# Initial parameter values 
	# Add ridge penalty if X is rank-deficient
	v <- try(solve(XtX), silent = TRUE)
	if(inherits(v, "try-error")){
	  warning('A ridge penalty had to be used to calculate the inverse crossproduct of the predictor matrix')
	  eps <- 1e-05
	  if (ridge.CV == TRUE) eps<-cv.glmnet(x,y, alpha = 0)$lambda.min
	  beta <- drop(backsolve(XtX + diag(eps,p), xy))
	  invA <- solve(XtX + diag(eps,p))
	} else{
	  beta <- drop(backsolve(XtX + diag(nrow=p), xy))
	  invA <- v
	}
	residue <-drop(y-x%*%beta)
	sigma2 <- drop((t(residue) %*% residue) / n)
	u<-1/abs(beta) # Double Pareto Mode
	if (lambda.estimate == 'AP') {
	  cat("Finding the best lambda using the apriori method. The Gibbs sampler will start soon...\n")
	  lambda<-hyper_par_BayesRLasso(x, y)
	} else{
	    lambda <- rgamma(1, shape = 2*p, rate = sum(1/abs(beta)))
	}
	
	
	# Empirical Bayes update of lambda
	  if(lambda.estimate=='EB'){
	    a <-0;
	    b <-0;
	    tol <-10^-3
	    L0 <- Q(p, u, lambda)
	  }
	
	
	#################
	# MCMC SAMPLING #
	#################
	
  # Counter for the MCMC sampling iterations
	k <- 0
	
	# Track time
	cat(c("Job started at:",date()),fill=TRUE)
	start.time <- Sys.time()
	
	# Gibbs sampler
	while (k < max.steps) {
	  k <- k + 1
	  
	  if (k %% 1000 == 0) {
	    cat('Iteration:', k, "\r")
	}
	  
	  # Sample beta
		mean <- drop(invA%*%xy)
		varcov <- sigma2*invA
		truncation.point = 1/u
		beta<- as.vector(t(rtmvnorm_midtruncated(1,		                    
		                                           Mean = mean,
		                                           Sigma = varcov,
		                                           lower = - truncation.point,
		                                           upper = truncation.point)))
		betaSamples[k,] <- beta

		# Sample sigma2
		if (update.sigma2 == TRUE){
		  shape.sig <- (n-1)/2
		  residue <- drop(y -x%*%beta)
		  rate.sig <- (t(residue) %*% residue)/2
		  sigma2 <- 1/rgamma(1, shape.sig, 1/rate.sig)
		}
		sigma2Samples[k] <- sigma2

		# Sample u
		T<- abs(1/beta)
		lambdaPrime <- lambda
		for (i in 1:p){
		  u[i] <- rexp(1,lambdaPrime) + T[i]
		  }
		uSamples[k, ] <- u

		# Update lambda
		  if(lambda.estimate=='EB'){
		    lambda1 <- (2*p)/sum(u)
		    L1 <-Q(p, u, lambda1)
		    if(abs(L0-L1)>tol){
		      lambda <- lambda1
		      L0 <- L1
		    } else {
		      lambda <- lambda
		      L0 <- L0
		    }
		  } 
		if(lambda.estimate=='MCMC'){
		    shape.lamb =  a + 2*p
		    rate.lamb = b + sum(u)
		    lambda <- rgamma(1, shape=shape.lamb, rate=rate.lamb)
		  }
		lambdaSamples[k] <- lambda
}
	
	# Collect all quantities and prepare output
	beta.post=betaSamples[seq(n.burn+1, max.steps, n.thin),]
	sigma2.post=sigma2Samples[seq(n.burn+1, max.steps, n.thin)]
	lambda.post=lambdaSamples[seq(n.burn+1, max.steps, n.thin)]

	# Posterior estimates of beta
	if (posterior.summary.beta == 'mean') {
	  beta <-apply(beta.post, 2, mean, na.rm = na.rm)
	} else if (posterior.summary.beta == 'median') {
	  beta <-apply(beta.post, 2, median, na.rm = na.rm)
	} else if (posterior.summary.beta == 'mode') {
	  beta <-posterior.mode(mcmc(beta.post), na.rm = na.rm)
	} else{
	  stop('posterior.summary.beta must be either mean, median, or mode.')
	}
	
	# Posterior estimates of lambda
	if (posterior.summary.lambda == 'mean') {
	  lambda<-mean(lambda.post, na.rm = na.rm)
	} else if (posterior.summary.lambda == 'median') {
	  lambda<-median(lambda.post, na.rm = na.rm)
	} else if (posterior.summary.lambda == 'mode') {
	  lambda <-posterior.mode(mcmc(lambda.post), na.rm = na.rm)
	} else{
	  stop('posterior.summary.lambda must be either mean, median, or mode.')
	}
	
	# Credible intervals for beta
	lowerbeta<-colQuantiles(beta.post,prob=(1 - beta.ci.level)/2, na.rm = na.rm)
	upperbeta<-colQuantiles(beta.post,prob=(1 + beta.ci.level)/2, na.rm = na.rm)
	
	# Assign names
	names(beta)<-names(lowerbeta)<-names(upperbeta)<-colnames(beta.post)<-colnames(x)
	
	# Credible intervals for lambda
	lambdaci<-credInt(lambda.post, cdf = NULL,conf=lambda.ci.level, type="twosided")
	
	# Track stopping time
	stop.time <- Sys.time()
	time<-round(difftime(stop.time, start.time, units="min"),3)
	cat(c("Job finished at:",date()),fill=TRUE)
	cat("Computational time:", time, "minutes \n")
	
	# Return output
	return(list(time = format(time),
	            beta = beta,
	            lowerbeta = lowerbeta,
	            upperbeta = upperbeta,
	            lambda = lambda,
	            lambdaci = lambdaci,
	            beta.post = beta.post,
	            sigma2.post = sigma2.post,
	            lambda.post = lambda.post))
}


