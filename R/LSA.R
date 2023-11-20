##################
# Load libraries #
##################

library(expm)
library(tidyverse)

###############################################
# Given model object get LSA-transformed data #
###############################################

LSA<-function(obj){
  
  #############################
  # Bare minimum sanity check #
  #############################
  
  supported_models<-c('lm', 'glm', 'negbin', 'zeroinfl', 'coxph')
  if (!class(obj)[1] %in% supported_models){
    stop(paste('Object not recognized - supported model classes are:', paste(supported_models,collapse=" ")))
  }
  
  ##############################################
  # Extract MLE and VCOV from the model object #
  ##############################################
  
  b0<-coefficients(obj) 
  Sigma0 <- solve(vcov(obj))
  
  ####################
  # Handle intercept #
  ####################
  
  intercept_index<-which(grepl("Intercept", names(b0)))
  
  if (rlang::is_empty(intercept_index)){
    b<-b0 # Working Beta
    Sigma<-Sigma0 # Working Sigma
    } else {
    a11 <- Sigma0[intercept_index, intercept_index] # a11 <- Sigma0[1,1]; n1 <- dim(Sigma0)[1];
    a12 <- Sigma0[-intercept_index, intercept_index] # a12 <- Sigma0[2:n1,1]
    a22 <- Sigma0[-intercept_index, -intercept_index] # a22 <- Sigma0[2:n1,2:n1]
    b <- b0[-intercept_index] # Working Beta 
    if (length(intercept_index)==1){
      Sigma<-a22-outer(a12,a12)/a11  # Working Sigma
    } else{ 
      Sigma<-a22-a12%*%solve(a11)%*%t(a12) # Working Sigma
      }
    beta0 <- b0[intercept_index] # Working Intercept
  }
  
  #########################
  # Get x_star and y_star #
  #########################
  
  cov.star<-sqrtm(Sigma) 
  y.star<-as.vector(cov.star%*%b) # Transformed y
  x.star<-as.matrix(cov.star) # Transformed x matrix
  colnames(x.star)<-colnames(Sigma)
  
  ##########################
  # Return relevant output #
  ##########################
  
  if (rlang::is_empty(intercept_index)){
    return(list(y.star = y.star, 
                x.star = x.star,
                Sigma = Sigma))
    } else{
      return(list(y.star = y.star, 
                  x.star = x.star, 
                  b = b,
                  beta0 = beta0,
                  a12 = a12,
                  a22 = a22,
                  Sigma = Sigma))
      }
}

##############################################
# Append intercept to LSA-based coefficients #
##############################################

append_intercept<-function(post_LSA_beta, LSA_obj){
  
  # Reset the name
  Sigma<-LSA_obj$Sigma
  names(post_LSA_beta)<-colnames(Sigma)
  
  if (length(LSA_obj)!=3){
    beta<-post_LSA_beta
    b<-LSA_obj$b
    beta0<-LSA_obj$beta0
    a12<-LSA_obj$a12
    a22<-LSA_obj$a22
    beta0 <- beta0 + t(a12)%*%solve(a22)%*%(beta - b) # Conditional mean of Normal
    post_LSA_beta<-c(beta0, post_LSA_beta)
    names(post_LSA_beta)[1:length(beta0)]<-'(Intercept)'
  }
  return(post_LSA_beta)
}