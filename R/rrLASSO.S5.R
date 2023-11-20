#' Bayesian Reciprocal Regularization
#'
#' Simplified Shotgun Stochastic Search with Screening (S5) algorithm to
#' solve the reduced reciprocal LASSO, based on a source code  
#' adapted from `BayesS5` by Minsuk Shin. 
#' 
#' @param X A numeric matrix with standardized predictors in columns and samples in rows.
#' @param y A mean-centered continuous response variable with matching rows with X.
#' @param S A screening size of variables. Default is 30.
#' @param lam A tuning parameter for the rrLASSO objective function. Default is 1.
#' @param intercept If TRUE, intercept is included in the final OLS fit. The default is TRUE.
#' @param verbose If TRUE, the function prints the current status of the S5 in each temperature. The default is FALSE.
#' @param seed Seed value for reproducibility. Default is 1234.
#' 
#' @return A list containing the following components is returned:

#' \item{hppm}{Index of the highest posterior probablity model.}
#' \item{marg.prob}{Marginal posterior probablities of the coefficients.}
#' \item{beta}{OLS coefficients from the final rLASSO model.}
#' \item{time}{Computation time in seconds.}

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
#' #################################
#' # Reduced Reciprocal LASSO (S5) #
#' #################################
#' 
#' rrLasso<- rrLASSO.S5(x.train, y.train)
#' y.pred.rrLasso<-x.test%*%rrLasso$beta[-1]
#' mean((y.pred.rrLasso - y.test)^2) # Performance on test data
#'
#' }
#' 
#' @export
#' @keywords Bayesian regularization, Reciprocal LASSO, S5

rrLASSO.S5 = function(X,
                      y,
                      S = 30, # Screening size. Default is 30.
                      lam = 1, # Tuning parameter for tau in rLASSO. 
                      intercept = TRUE, 
                      verbose = FALSE, 
                      seed = 1234){
  
  ##########################################
  # Set some intial parameters and options #
  ##########################################
  
  set.seed(seed)
  r0=1;a0=0.01;b0=0.01
  ind_fun = ind_fun_rLASSO
  tau = 1
  n = dim(X)[1]
  p = dim(X)[2]
  ITER=20
  IT=20
  tem =  seq(0.4,1,length.out=IT)^2
  IT.seq = rep(ITER,IT)
  A3 = S
  gam = rep(0,p)
  p.g=sum(gam)
  ind2= which(gam==1)
  
  
  ###############################
  # Set up screening parameters #
  ###############################
  
  curr = -1000000
  GAM = rep(1,p)
  OBJ = -1000000
  OBJ.id = -1000000
  OBJ.m0 = rep(-1000000,p)
  OBJ.p0 = rep(-1000000,p)
  ID = -100
  ID.obj = -100
  it=1
  if(p.g>0){
    fit = solve(crossprod(X[,ind2])+diag(p.g))%*%crossprod(X[,ind2],y)
    res = y-X[,ind2]%*%fit
    corr = as.vector(crossprod(res,X))
    ind.ix = sort.int(abs(corr),decreasing=TRUE,index.return=TRUE)$ix
    s = c(ind2,ind.ix)
  }else{res=y
  corr = as.vector(crossprod(res,X))
  ind.ix = sort.int(abs(corr),decreasing=TRUE,index.return=TRUE)$ix
  s = ind.ix
  }
  
  if(p<50){p00=10}else{p00=round(p/10)}
  size = A3
  IND = s[1:size]
  p.ind = length(IND)
  C.p = rep(-1000000,p)
  C.m = rep(-1000000,p)
  obj=-10000    
  
  ############################
  # Set up parameters for S5 #
  ############################
  
  tau = lam
  B = lam
  GAM = rep(1,p)
  OBJ = -1000000
  OBJ.id = -1000000
  OBJ.m0 = rep(-1000000,p)
  OBJ.p0 = rep(-1000000,p)
  ID = -100
  ID.obj = -100
  it=1
  ID0 = NULL
  INT = NULL
  K=0
  pmt = proc.time()
  for(it in 1:IT){   
    IT0 = IT.seq[it] 
    pq=0
    for(iter in 1:IT0){
      id = sum(2^(2*log(ind2)))
      id.ind = which(id==ID)
      leng = length(id.ind)
      
      if(leng==0){
        ID = c(ID,id)
        C.p = rep(-100000,p)
        for(i in (p.g+1):p.ind){
          j=IND[i]
          gam.p = gam;gam.p[j]=1;ind.p=which(gam.p==1)
          int  = ind_fun(ind.p,B,y,X,tau,a0,b0,p)     
          obj.p =  c(int)
          if(is.na(obj.p)==TRUE){obj.p = -100000}
          C.p[j] = obj.p
        }
        p.g = sum(gam)
        C.m = rep(-100000,p)
        IND.m = ind2
        p.ind.m = length(IND.m)
        for(i in 1:p.g){
          j=ind2[i]
          gam.m = gam;gam.m[j]=0;ind.m=which(gam.m==1)     
          int  = ind_fun(ind.m,B,y,X,tau,a0,b0,p)
          obj.m =  c(int)
          if(is.na(obj.m)==TRUE){obj.m = -100000}
          C.m[j] = obj.m  
        }    
        OBJ.p0 = cbind(OBJ.p0,C.p)
        OBJ.m0 = cbind(OBJ.m0,C.m)  
      }else{
        pq= pq+1
        C.p = OBJ.p0[,(id.ind[1])];C.m = OBJ.m0[,(id.ind[1])]
      }
      
      prop = exp(tem[it]*(C.p-max(C.p)))
      sample.p = sample(1:length(prop),1,prob=prop)
      obj.p = C.p[sample.p]
      prop = exp(tem[it]*(C.m-max(C.m)))
      sample.m = sample(1:length(prop),1,prob=prop)
      obj.m = C.m[sample.m]
      
      l = 1/(1+exp(tem[it]*obj.m-tem[it]*obj.p))
      if(l>runif(1)){ gam[sample.p]=1;obj = obj.p;curr=obj.p
      }else{
        gam[sample.m]=0;obj = obj.m;curr=obj.m
      } 
      ind2 = which(gam==1)
      p.g = sum(gam)
      
      
      
      if(p.g>0){
        fit = solve(crossprod(X[,ind2])+diag(p.g))%*%crossprod(X[,ind2],y)
        res = y-X[,ind2]%*%fit
        corr = as.vector(crossprod(res,X))
        ind.ix = sort.int(abs(corr),decreasing=TRUE,index.return=TRUE)$ix
        s = c(ind2,ind.ix)
      }else{res=y
      corr = as.vector(crossprod(res,X))
      ind.ix = sort.int(abs(corr),decreasing=TRUE,index.return=TRUE)$ix
      s = ind.ix
      }   
      size = A3
      IND = s[1:size];
      p.ind = length(IND)
      
      id = sum(2^(2*log(ind2)))
      id.ind = which(id==ID.obj)
      
      leng = length(id.ind)
      if(leng==0){
        ID.obj = c(ID.obj,id)
        OBJ = c(OBJ,obj)
        GAM= cbind(GAM,gam)
      }
      gam = rep(0,p);gam[ind2]=1
      #if(crossprod(gam-xd0)==0){K=1}
    }
    
    if(length(OBJ)>100){
      prop = exp(tem[it]*(OBJ-max(OBJ)))
      sample.p = sample(1:length(prop),5,replace=FALSE,prob=prop)
      gam = GAM[,sample.p[1]]
      gam[which(GAM[,sample.p[2]]==1)]=1 
      gam[which(GAM[,sample.p[3]]==1)]=1
      gam[which(GAM[,sample.p[4]]==1)]=1
      gam[which(GAM[,sample.p[5]]==1)]=1
      ind2 = which(gam==1);p.g = sum(gam)
      id = sum(2^(2*log(ind2)))
      id.ind = which(id==ID.obj)
      leng = length(id.ind)
      if(leng==0){ID.obj = c(ID.obj,id)
      curr = ind_fun(ind2,B,y,X,tau,a0,b0,p)
      if(is.na(curr)==TRUE){curr=-100000}
      OBJ= c(OBJ,curr)
      GAM = cbind(GAM,gam)
      }else{curr=OBJ[id.ind[1]]}
    }
    ind2 = which(gam==1)
    p.g = sum(gam)
    if(p.g>0){
      fit = solve(crossprod(X[,ind2])+diag(p.g))%*%crossprod(X[,ind2],y)
      res = y-X[,ind2]%*%fit
      corr = as.vector(crossprod(res,X))
      ind.ix = sort.int(abs(corr),decreasing=TRUE,index.return=TRUE)$ix
      s = c(ind2,ind.ix)
    }else{res=y
    corr = as.vector(crossprod(res,X))
    ind.ix = sort.int(abs(corr),decreasing=TRUE,index.return=TRUE)$ix
    s = ind.ix
    }   
    
    size = A3
    IND = s[1:size];
    p.ind = length(IND)
    
    gam.pr = GAM[,which.max(OBJ)]
    obj.pr = max(OBJ)
    ind2.pr = which(gam.pr==1)
    if(verbose==TRUE){
      curr = ind_fun(ind2,B,y,X,tau,a0,b0,p)  
      mse=crossprod(lm(y~X[,ind2.pr])$residuals)/n
      print("Inverse Temperature");print(tem[it]);print("The Selected Variables in the Searched MAP Model");print(ind2.pr);print("The Evaluated Object Value at the Searched MAP Model");print(obj.pr);
      print("Current Model");print(ind2);  print("The Evaluated Object Value at the Current Model");print(curr);
      print("SSE of the Searched MAP model");print(mse);print("The Number of Total Searched Models");print(length(OBJ));print(pq)
    }
    
  }
  
  time  = proc.time() - pmt
  marg.gam = rep(0,p)
  for(u in 2:ncol(GAM)){
    marg.gam = marg.gam + GAM[,u]*exp(OBJ[u]-max(OBJ))
  }
  marg.gam = marg.gam / sum(exp(OBJ-max(OBJ)))
  gam0 = GAM[,which.max(OBJ)]
  ind2 = which(gam0==1)
  post = exp(OBJ-max(OBJ))/sum(exp(OBJ-max(OBJ)))
  hppm = 1/sum(exp(OBJ-max(OBJ)))
  
  ########################################################################
  # Calculate OLS coefficients for final set of coefficients (debiasing) #
  ########################################################################
  
  hppm = which(gam0==1)
  marg.prob = marg.gam
  beta = rep(0, ncol(X))
  if (length(hppm)>=1){
    X_hppm<-as.matrix(X[,hppm])
    if (intercept){
      ols.lm = lm(as.vector(y) ~ X_hppm) 
      beta[hppm]<-ols.lm$coefficients[-1]
      beta<-c(ols.lm$coefficients[1], beta)
    } else{
      ols.lm = lm(as.vector(y) ~ X_hppm - 1) 
      beta[hppm]<-ols.lm$coefficients
    }
  } else{
    if (intercept){
      ols.lm = lm(as.vector(y) ~ as.matrix(X)) 
    } else{
      ols.lm = lm(as.vector(y) ~ as.matrix(X)-1) 
    }
    beta<-ols.lm$coefficients
  }
  if (intercept){
    names(beta)<-c('(Intercept)', colnames(X))
  } else {
    names(beta)<-colnames(X)
  }
  
  #################
  # Return Output #
  #################
  
  return(list(hppm = hppm, 
              marg.prob = marg.prob,
              beta = beta,
              time = time))
  
}


ind_rLASSO= function(x,ind2,B,y,X){
  p.g = length(x)
  if(p.g>1){
    ress = crossprod(y-X[,ind2]%*%x)
    int = (ress) + sum(B/(abs(x)))
  }else{
    if(p.g==1){
      ress = crossprod(y-X[,ind2]*x)
      int = (ress) + (B/abs(x))
    }
  }
  return(int)
} 


ind_fun_rLASSO=function(ind2,B,y,X,tau,a0,b0,p){ 
  a0=0.01;b0=0.01
  tau = 1
  p.g=length(ind2)
  sb = lbeta(1+p.g,1+p-p.g)
  if(p.g >1){
    fit = solve(crossprod(X[,ind2]))%*%crossprod(X[,ind2],y)
    ress = crossprod(y-X[,ind2]%*%fit)
    initial_x = fit
    wrapper <- function(theta) ind_rLASSO(theta,ind2,B,y,X)
    o <- optim(initial_x, wrapper,method="L-BFGS-B" ) 
    f_x = o$value
    x0 = o$par
    
    int = -1*f_x
    
  }else{
    
    if(p.g==1){
      fit = ((crossprod(X[,ind2])+1)^-1*crossprod(X[,ind2],y))[1]
      ress = crossprod(y-X[,ind2]*fit)
      initial_x = fit
      wrapper <- function(theta) ind_rLASSO(theta,ind2,B,y,X)
      o <- optim(initial_x, wrapper,method="L-BFGS-B") 
      f_x = o$value
      int = -1*f_x
      
    }else{if(p.g==0){
      int = -1*crossprod(y) 
    }
    }
  }
  return(int)
}

