

gamAR <- function (formula, data, p.ar=1, starts = starts, w = rep(1,NN),
                       family = "poisson",  cc = 0.5,de = 0.01,
                       control=list(...),...)
{
  
  
  if(family!="poisson")
    stop("sorry, only poisson family is currently implemented")
  family=poisson()
  control <- do.call("glm.control", control)
  times <- control$maxit
  epsilon <- control$epsilon
  if (missing(data)) 
    data <- environment(formula)
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  term.labels <- attr(attributes(mf)$terms,"term.labels")
  mt <- attr(mf, "terms")
  y <- model.response(mf, "any")
  
  X <- model.matrix(mt, mf, contrasts)
  xnames <- colnames(X)
  ynames <- names(y)
  
  b <-ncol(X)-1
  p <- p.ar
  NN <- NROW(y)
  pp <-rep(0,p)
  
  r<-matrix(0,times+1,1+b+p)
  
  
  r[1,]<-c(starts,pp)
  
  
  
  #aa is part of design matrix for AR resar pearson residuals 
  aa<-matrix(0,NN-p,p)
  #Q is partial derivative of eta (linear predictor) on every parameter 
  Q<-matrix(0,NN-p,b+p+1);
  QQ<-matrix(0,NN-p,b+1)
  #transpose of X
  tX<-t(X)
  #all are variables used in the calculation                        
  a<-mm<-n<-u<-numeric(NN);
  
  
  
  #yy=max(cc,y) avoid that y<0
  yy<-y                    
  yy[y<cc]<-cc    
  pan <- 0
  
  
  
  
  k=1        
  while(k<(times+1)&pan==0)
  {
    
    if(k==1)
      #while in the first iteration
      #give the initial parameters.
    { 
      bb<-r[k,c(1:(b+1))];
      pp<-r[k,c((b+2):(b+p+1))];           
      #calculate a, the intermediate variabel a
      a<-log(yy)-X%*%bb;
      for(i in 1:p)
        aa[,p-i+1]<-a[i:(NN-p+i-1)]     
      n[(p+1):NN]<-(X%*%bb)[(p+1):NN]+aa%*%pp+log(w)[(p+1):NN]
      u<-exp(n)          
      #log partial likelihood lpl=ln(L)    
      lple<-y[(1+p):NN]*n[(1+p):NN]-u[(1+p):NN]
      lpl<-sum(lple)
      lpl0<-lpl;
    }
    
    
    #Q are partial derivatives of eta on each coefficient
    #Since the first p time points don't have the full AR terms 
    #Q matrix begin at time (p+1), Q[t,i] is the partial derivative of eta on the ith coefficient at time point (t+p).
    rp<-rev(pp)    
    for(t in 1:(NN-p))
      QQ[t,]<-as.matrix(tX[,t:(t+p-1)])%*%rp  
    Q[,1:(b+1)]<-X[(1+p):NN,]-QQ                
    Q[,(2+b):(1+b+p)]<-aa
    
    
    #dpl is the partial derivative of partial likelihood on each coefficient
    edpl<-(y[(1+p):NN]-u[(1+p):NN])*Q;
    dpl<-apply(edpl,2,sum)                
    
    
    #tt is fisher information matrix
    tt<-matrix(0,b+p+1,b+p+1); 
    for(t in 1:(NN-p))
    {tt[1:(1+b),1:(1+b+p)]<-tt[1:(1+b),1:(1+b+p)]+u[t+p]*Q[t,1:(1+b)]%*%t(Q[t,1:(1+b+p)])
     tt[(b+2):(b+p+1),(b+2):(b+p+1)]<-tt[(b+2):(b+p+1),(b+2):(b+p+1)]+u[t+p]*Q[t,(b+2):(b+p+1)]%*%t(Q[t,(b+2):(b+p+1)])
     for(j in (b+2):(b+p+1))
       tt[1:(1+b),j]<-tt[1:(1+b),j]+(y[t+p]-u[t+p])*X[t+p+b+1-j,1:(1+b)]
    } 
    tt[(2+b):(1+b+p),1:(1+b)]<-t(tt[1:(1+b),(2+b):(1+b+p)])
    
    
    #eigen decomposition of tt  
    ev<-eigen(tt); 
    #val is the eigen values of tt                     
    val<-ev$values;                     
    #vec is the eigen vector matrix of tt
    vec<-ev$vectors;                                                   
    
    # if all eigen values are larger than de, tt1 is the inverse matrix of tt
    # else 
    if(all(val>de)) 
      tt1<-solve(tt)    else  
      {  tt[c(1:(b+1)),c((b+2):(b+p+1))]<-matrix(0,b+1,p);
         tt[c((b+2):(b+p+1)),c(1:(b+1))]<-matrix(0,p,b+1);
         tt1<-solve(tt); 
      }  
    
    #l indicates whether it is the first time to change the coefficients in this iteration
    l<-1;                    
    #to get the new coefficients, add adr to the former coefficients
    adr<-tt1%*%dpl;                     
    
    
    #when the mode of adr is less than epsilon, make r[k+1]=r[k]
    if(t(adr)%*%adr>epsilon)      
    {        #when the new lpl is larger than the old lpl, then r[k+1]= new coefficients
      while(lpl<=lpl0)                    
      { if(l==1)                         
      {
        #new r=old r+adr
        r[k+1,]<-r[k,]+adr
        #change the indicator l     
        l=l+1;
      } else      r[k+1,]<-r[k,]+runif(1,0,1.5)*adr;
        #new r=old r+adr   
        #get lpl for comparison 
        
        bb<-r[k+1,c(1:(b+1))];
        pp<-r[k+1,c((b+2):(b+p+1))];           
        a<-log(yy)-X%*%bb;
        for(i in 1:p)
          aa[,p-i+1]<-a[i:(NN-p+i-1)]     
        n[(p+1):NN]<-(X%*%bb)[(p+1):NN]+aa%*%pp+log(w)[(p+1):NN]
        u<-exp(n)             
        lple<-y[(1+p):NN]*n[(1+p):NN]-u[(1+p):NN]
        lpl<-sum(lple)
      }
    }  else
    {r[k+1,]<-r[k,];
     pan<-1;
    }
    lpl0<-lpl;k=k+1
  }
  
  
  fit<-list()
  Estimate <-r[k,]
  Std.Error <-sqrt(diag(tt1))
  zvalue <- Estimate/Std.Error
  Pr <-  2 * pnorm(-abs(zvalue))
  fit$coefficients <- cbind(Estimate, Std.Error, zvalue, Pr)
  dimnames(fit$coefficients) <-list(c(xnames,paste("AR", 1:p.ar, sep = "")),c("Estimate","Std.Error","z value","Pr(>|z|)"))
  
  names(n) <- 1:NN
  names(u) <- 1:NN
  
  n <- n[-(1:p)]
  u <- u[-(1:p)]
  
  fit$linear.predictor <- n
  fit$fitted.values <- u
  y <- y[-(1:p)]
  w <- w[-(1:p)]
  dev.resids <- family$dev.resids
  aic <- family$aic
  Pearson.res <- (u-y)/sqrt(u)
  phi <-sum(Pearson.res^2)/(NN-1-b-2*p)
  
  
  
  fit$aic <- aic(y, mu=u, wt=w) + 2 * (1+b+p)
  fit$dev <- sum(dev.resids(y, mu=u, wt=w))
  fit$Pearson.residuals <- Pearson.res
  fit$phi <- phi
  fit$X <- X
  fit$loglikelihood <- -(1/2)*aic(y, mu=u, wt=w)
  ll<- c("Intercept",term.labels,"AR")
  fc <- c(attr(X, "assign"),rep(max(attr(X,"assign"))+1,p))
  aaa <- factor(fc, labels = ll)
  asgn <- split(order(fc), aaa)
  nterms <- length(asgn)
  
  predictor <- matrix(ncol = nterms, nrow = NN)
  dimnames(predictor) <- list(rownames(X), names(asgn))
  
  aab<-rbind(matrix(0,p,p),aa)
  X1 <-cbind(X,aab)
  
  for (i in seq(1,nterms)) {
    iipiv <- asgn[[i]]   
    predictor[, i] <-  
      X1[, iipiv, drop = FALSE] %*% Estimate[iipiv]
  }
  
  
  fit$term.predictors <- predictor[-(1:p),] 

  fit$pan <- pan
  
  fit
  
}
