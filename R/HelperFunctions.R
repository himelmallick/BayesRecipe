############################################################
# Multivariate normal samples under rectangular constraint #
# Extracted and modified from the mombf package            #
############################################################

rtmvnorm_midtruncated <- function(n, 
                                  Mean, 
                                  Sigma, 
                                  lower, 
                                  upper, 
                                  within = FALSE, 
                                  method='Gibbs', 
                                  burnin=0) {
  # Input
  # - n: number of draws
  # - Mean: multivariate normal mean
  # - Sigma: multivariate normal covariance
  # - lower: vector with lower truncation points
  # - upper: vector with upper truncation points
  # - within: if TRUE, each variable is truncated to be >=lower and <= upper. If FALSE, it's truncated to be <lower or >upper
  # - method: set method=='Gibbs' for Gibbs sampling, and method=='MH' for independent proposal MH
  # - burnin: number of burn-in iterations
  # Output: n draws obtained via Gibbs sampling after orthogonalization
  if (length(lower)==1) lower <- rep(1,length(Mean))
  if (length(upper)==1) upper <- rep(1,length(Mean))
  if (length(lower)!=length(Mean)) stop('Length of lower and Mean do not match')
  if (length(upper)!=length(Mean)) stop('Length of upper and Mean do not match')
  if (nrow(Sigma)!=length(Mean) | ncol(Sigma)!=length(Mean)) stop('Dimensions of Mean and Sigma do no match')
  if (!(method %in% c('Gibbs','MH'))) stop('Method should be Gibbs or MH')
  method <- as.integer(ifelse(method=='Gibbs',1,2))
  ans <- .Call("rtmvnormCI",
               as.integer(n), 
               as.double(Mean), 
               as.double(Sigma), 
               as.double(lower), 
               as.double(upper), 
               as.integer(within), 
               method, 
               PACKAGE = 'mombf')
  matrix(ans,ncol=length(Mean))
}

##############################################
# Helper function for data-adaptive lambda   #
# selection adapted from the BayesS5 Package #
##############################################

hyper_par_BayesRLasso <-function(x,
                                 y,
                                 threshold = NULL){
  
  # Extract data dimensions
  n =nrow(x)
  p =ncol(x)
  
  # Set threshold
  threshold<-ifelse(is.null(threshold), p^-0.5, threshold)
  
  # Generate null distrbution of beta
  betas = matrix(0,3,50000)
  for(k in 1:50000){
    sam = sample(1:p,3)
    ind = sample(1:n,n)
    betas[,k] = as.vector(solve(crossprod(x[ind,sam]))%*%crossprod(x[ind,sam],y))
  }
  
  # Calculate the optimal value of lambda
  res=y
  corr = as.vector(cor(res,x))
  ind.ix = sort.int(abs(corr),decreasing=TRUE,index.return=TRUE)$ix
  s = ind.ix[1:3]
  beta.hat =solve(crossprod(x[,s]))%*%crossprod(x[,s],y)
  sig.hat = crossprod(y - x[,s]%*%beta.hat)/n
  betas=as.vector(betas)
  lambda.cand = seq(0.1,(sd(y)+0.1),length.out=300)^2
  pro = rep(0,300)
  for(k in 1:300){
    lambda = lambda.cand[k]
    den = function(x){lambda*x^-2*exp(-1*lambda/abs(x))/2} # rLASSO Density
    den.null1 = density(betas)
    data = list(x=den.null1$x,y=den.null1$y)
    den.null = approxfun(data[[1]], data[[2]], rule=1,method = "linear")
    f = function(x){den(x) - den.null(x)}
    tryCatch({
      th=1
      a = uniroot(f,interval = c(0.001,max(betas))) 
      th = a$root
      loc = integrate(den.null,lower = th, upper =max(betas)-0.001)$value
      nonloc =  integrate(den,lower = 0, upper = th)$value
      pro[k] = loc + nonloc}, error=function(e){})
  }
  lambda=1
  B = lambda.cand[which.min((pro-threshold)^2)]
  return(B)
}

########################################################
# Objective Function for Bayesian rLASSO EB Estimation #
########################################################

Q<-function(p, 
            u, 
            lambda){
  return(2*p*log(lambda)-lambda*sum(u))
}


###############################################################
# Prediction in new samples by Bayesian Model Averaging (BMA) #
###############################################################

predict.BMA<-function(beta.post, 
                      x.test, 
                      na.rm = TRUE){
  
  # Sanity check
  if (ncol(beta.post) != ncol(x.test)) stop('x.test and beta.post should have matching columns.')
  
  # Initialize
  y.pred.mat<-matrix(nrow = nrow(beta.post), ncol = nrow(x.test))
  
  # Calculate
  for (i in 1:nrow(beta.post)){
    y.pred.mat[i,]<-x.test%*%beta.post[i,]
  }
  
  # Summarize and return
  y.pred<-as.vector(colMeans(y.pred.mat, na.rm = na.rm))
  return(y.pred)
}

############################################################
# Prediction in new samples by a vanilla posterior summary #
############################################################

predict.vanilla<-function(beta.post, 
                          x.test, 
                          posterior.summary.beta = 'mean',
                          na.rm = TRUE){
  
  # Sanity check
  if (ncol(beta.post) != ncol(x.test)) stop('x.test and beta.post should have matching columns')
  
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
  
  # Summarize and return
  y.pred<-x.test%*%beta
  return(y.pred)
}

##########################################
# Variable Selection by Back Propagation #
##########################################

variable.selection.BP<-function(x, 
                                y, 
                                penalty ='rLASSO',
                                beta.post,
                                lambda.post,
                                posterior.summary.beta = 'mean',
                                posterior.summary.lambda = 'median',
                                na.rm = TRUE){
  
  # Sanity check
  if (ncol(beta.post) != ncol(x)) stop('x and beta.post should have matching columns')
  
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
  
  # Summarize and return 
  if (penalty == 'rLASSO'){
    rLASSO.sol= rrLASSO.S5(x, y, lam = lambda)  
    newbeta = rLASSO.sol$beta[-1]
  } else if (penalty == 'LASSO'){
    LASSO.sol = glmnet(x, y, alpha = 1)    
    newbeta = predict(LASSO.sol, type="coef", s=lambda)[-1]
  } else {
    stop('penalty must be either LASSO or rLASSO.')
  }
    
  # Save and Return
  beta[which(newbeta==0)]<-0
  beta.bayes<-beta
  beta.freq<-newbeta
  return(list(beta.bayes = beta.bayes, beta.freq = beta.freq))
}

#################################################
# Variable Selection by Median Probablity Model #
#################################################

variable.selection.MPM<-function(x, 
                                 y, 
                                 beta.post, 
                                 lambda.post, 
                                 posterior.summary.beta = 'mean', 
                                 cutoff = 0.5, 
                                 na.rm = TRUE){
  
  # Sanity check
  if (ncol(beta.post) != ncol(x)) stop('x and beta.post should have matching columns')
  
  # Initialize
  beta.mat<-matrix(0, nrow = nrow(beta.post), ncol = ncol(beta.post))
  
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
  
  # Calculate sparse models for varying lambda
  for (i in 1:length(lambda.post)){
    fit_rLASSO = rrLASSO.S5(x, y, lam = lambda.post[i])
    beta.mat[i,]<-ifelse(fit_rLASSO$beta[-1]!=0, 1, 0)
  }
  
  # Summarize and return
  beta.mat.summary<-colMeans(beta.mat)/nrow(beta.mat)
  beta.index<-which(beta.mat.summary>cutoff)
  beta[!beta.index]<-0
  return(list(beta = beta, beta.mat = beta.mat))
}

###################################
# Variable Selection based on DSS #
###################################

variable.selection.DSS<-function(x, 
                                 beta.post, 
                                 posterior.summary.beta = 'mean',
                                 lambda.select = 'min', 
                                 na.rm = TRUE){

    
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
    
    ##############################
    # Only run if there is no NA #
    #############################
    
  if (any(!is.na(beta))){
    # Perform variable selection using DSS
    y.pred<-x%*%beta
    alasso.cv = glmnet::cv.glmnet(x, y.pred, alpha = 1, penalty.factor=1/abs(beta)) 
    if (lambda.select == 'min'){
      lambda.cv.alasso.dss = alasso.cv$lambda.min
    } else if (lambda.select == '1se'){
      lambda.cv.alasso.dss = alasso.cv$lambda.1se          
    } else {
      stop('posterior.summary.beta must be either mean, median, or mode.')
    }
    
    # Summarize and return
    alasso.sol.dss = glmnet::glmnet(x, y.pred, alpha = 1, penalty.factor = 1/abs(beta))    
    alasso.coeff.dss = glmnet::predict.glmnet(alasso.sol.dss, type="coef", s=lambda.cv.alasso.dss)
    nonactive.set.dss <- which(alasso.coeff.dss[-1] == 0)
    beta[nonactive.set.dss] <- 0
  }
  return(beta)
}


##################################
# Variable Selection based on CI #
##################################
variable.selection.CI<-function(beta.post, 
                                posterior.summary.beta = 'mean',
                                beta.ci.level = 0.95,
                                na.rm = TRUE){
  
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
  
  # Credible intervals for beta
  lowerbeta<-colQuantiles(beta.post,prob=(1 - beta.ci.level)/2, na.rm = na.rm)
  upperbeta<-colQuantiles(beta.post,prob=(1 + beta.ci.level)/2, na.rm = na.rm)
  
  # Summarize and return
  beta[which(sign(lowerbeta) != sign(upperbeta))]<-0
  return(beta)
}


######################################
# Bayesian Lasso by Park and Casella #
######################################

##################################################################
# Code modified from https://cs.gmu.edu/~pwang7/gibbsBLasso.html #
##################################################################

BayesLassoSMN <-function(x, 
                         y, 
                         a = 1, 
                         b = 1,
                         max.steps = 11000,
                         n.burn = 1000,
                         n.thin = 1,
                         posterior.summary.beta = 'mean',
                         posterior.summary.lambda = 'median',
                         beta.ci.level = 0.95,
                         lambda.ci.level = 0.95,
                         seed = 1234,
                         na.rm = TRUE) {
  # Set random seed
  set.seed(seed)
  
  # FOR PRIOR
  X <- as.matrix(x)
  p<-ncol(X);
  n<-nrow(X);
  xx <- t(X)%*%X
  xy <- t(X)%*%y
  lambda.sq  <- rgamma(1,shape=a, rate=b) 
  sig.sq     <- runif(1,0.1,10)
  tau.sq     <- rexp(p, rate=lambda.sq/2)
  mean.be    <- rep(0,p)
  cov.be     <- sig.sq*diag(tau.sq)
  beta.l     <- rmvnorm(1,mean=mean.be,sigma=cov.be)
  
  
  # FOR POSTERIOR
  sigsq.post <- lambda.post <- NULL 
  tausq.post <- rbind( tau.sq,matrix(rep(NA,max.steps*p),ncol=p) )
  beta.p     <- rbind( beta.l,matrix(rep(NA,max.steps*p),ncol=p) )
  
  # POSTERIOR
  start.time <- Sys.time()
  cat(c("Job started at:",date()),fill=TRUE)
  # MCMC LOOPING
  for (M in 1:max.steps)  {
    # if (M %% 100 == 0) print(M)
    
    # Full Conditional Posteriors  
    # beta
    dtau.inv   <- diag(1/tau.sq)
    cov.be     <- sig.sq * solve(xx+dtau.inv) 
    mean.be    <- 1/sig.sq * cov.be%*%t(X)%*%(y)
    beta.p[M+1,] <- rmvnorm(1,mean=mean.be,sigma=cov.be)
    
    # tau.sq
    gam <- c()
    for (j in 1:p){
      repeat{
        gam[j]  <- rinvGauss(1, nu=sqrt(lambda.sq * sig.sq/beta.p[M+1,j]^2), lambda=lambda.sq)
        if (gam[j] > 0) break    	
      }
      tau.sq[j] <- 1/gam[j]
    }
    tausq.post[M+1,] <- tau.sq 	
    
    # sig.sq
    sh.sig     <- (n-1+p)/2
    sc.sig     <- 1/2*t(y-X%*%beta.p[M+1,])%*%(y-X%*%beta.p[M+1,])+ 1/2*t(beta.p[M+1,])%*%diag(1/tau.sq)%*%beta.p[M+1,]
    sig.sq     <- rinvgamma(1, shape=sh.sig, scale=sc.sig)
    sigsq.post <- c(sigsq.post, sig.sq)
    
    # lambda
    sh.lam      <- p + a 
    sc.lam      <- 1/2*sum(tau.sq) + b 
    lambda.sq   <- rgamma(1, shape=sh.lam, rate=sc.lam)
    lambda.post <- c(lambda.post, lambda.sq)
  }
  
  # Collect all quantities and prepare output
  beta.post <- beta.p[seq(n.burn+1, max.steps,n.thin ),]
  lamb.post<-sqrt(lambda.post[seq(n.burn+1,max.steps,1)])
  sigsq.post<-sigsq.post[seq(n.burn+1,max.steps,1)]
  
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
  return(list(time = time,
              beta = beta,
              lowerbeta = lowerbeta,
              beta.post = beta.post,
              upperbeta = upperbeta,
              lambda = lambda,
              lambdaci = lambdaci,
              sigma2.post = sigsq.post,
              lambda.post = lamb.post))
  
}

#################################################################
# Generate simulated datasets according to specified parameters #
#################################################################

generate_simdata<-function(ntrain, 
                           nfeature, 
                           ntest = 200,
                           sigma = 1.5,
                           rho = 0.5, 
                           design = 'AR',
                           model = 'shin',
                           nrep = 100, 
                           seed = 1234){
  
  # Set random seed
  set.seed(seed)
  
  # Initialize list of datasets
  dataset<-vector("list", nrep)
  
  # Initialize coefficients
  beta0 = rep(0, nfeature)
  
  # Set Model-specific ground truth
  if (model == 'shin'){
    if (nfeature<5) stop ('For Shin 2018 simulation model, p must be greater than 5!')
    beta0[1:5] = c(0.5, 0.75, 1.00, 1.25, 1.50)*sample(c(1,-1), 5, replace=TRUE)
  }
  
  if (model == 'tibsA'){
    beta0[1]<-3
    beta0[2]<-1.5
    beta0[5]<-2
  }
  
  if (model == 'tibsB'){
    beta0<-rep(0.85, nfeature)
  }
  
  if (model == 'tibsC'){
    beta0[1]<-5
  }

  # Residual variance
  sigma2 = sigma**2
  
  # Set Mu and Sigma of MVN
  mu<-rep(0, nfeature)
  cov<-diag(1, nfeature, nfeature)
  
  # Loop over nrep
  for (k in 1:nrep){
    
    # Design-specific covariance
    if (design == 'AR'){
      for (i in 1:nfeature){for (j in 1:nfeature){if(i!=j) cov[i,j]=rho**(abs(i-j))}}
    }
    if (design == 'CS'){
      for (i in 1:nfeature){for (j in 1:nfeature){if(i!=j) cov[i,j]=rho}}
    }
    
    # Generate X (train)
    X.train<-as.matrix(MASS::mvrnorm(n = ntrain, mu, cov))
    colnames(X.train)<-paste("V", 1:nfeature, sep ='')
    
    # Generate X (test)
    X.test<-as.matrix(MASS::mvrnorm(n = ntest, mu, cov))
    colnames(X.test)<-paste("V", 1:nfeature, sep ='')
    
    # Generate y (train)
    Y.train = X.train%*%beta0 + rnorm(ntrain)*sigma2 
    Y.train<-as.vector(Y.train)
    
    # Generate y (test)
    Y.test = X.test%*%beta0 + rnorm(ntest)*sigma2 
    Y.test<-as.vector(Y.test)
    
    # Store y and X
    dataset[[k]]<-list(Y.train = Y.train, 
                       Y.test = Y.test,
                       X.train = X.train,
                       X.test = X.test)
    
    # Reset rho if orthogonal design
    if (design == 'OR') rho = 0
  }
  return(list(dataset = dataset, 
              beta0 = beta0, 
              ntrain = ntrain, 
              ntest = ntest, 
              nfeature = nfeature, 
              sigma = sigma,
              rho = rho, 
              design = design,
              model = model,
              nrep = nrep, 
              seed = seed))
}


