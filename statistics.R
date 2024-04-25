# if(!require(glmnet)){
#   install.packages("glmnet")
# }
# if(!require(Matrix)){
#   install.packages("Matrix")
# }
# if(!require(SIS)){
#   install.packages("SIS")
# }
# if(!require(parallel)){
#   install.packages("parallel")
# }
# if(!require(hdi)){
#   install.packages("hdi")
# }
# 
# library(SIS)
# library(scalreg)
# library(parallel)
# library(hdi)

score_test <- function(v, y, test.set = NULL, model = 'linear', method = 'asymptotic',
                       multiplier = 'gaussian', M = 500, orthogonalization = FALSE, scale = FALSE, lambda = 0, threshold = TRUE) {
  if(scale) {
    v <- scale(v); y <- scale(y)
  }
  if(is.null(test.set)) {
    test.set <- 1:ncol(v)
  }
  x <- v[,test.set]; z <- v[,-test.set]
  n <- nrow(x); p.beta <- ncol(x); p.gamma <- ncol(z)
  if(lambda <= 0){
    lambda <- 1.25/sqrt(n)
  }
  
  ## group testing with L2 norm
  # derive epsilon
  if(p.gamma == 0){
    if(model == 'linear'){
      err <- as.numeric(y)
    }
    if(model == 'logit'){
      err <- as.numeric(y - 0.5)
    }
    if(model == 'poisson'){
      err <- as.numeric(y - 1)
    }
  }else{
    if(model == 'linear'){
      if(p.gamma <= n/10){
        gamma.hat <- lm.fit(as.matrix(z),y)$coefficients
      }else{
        cv.model <- cv.glmnet(z, y, intercept=FALSE)
        gamma.hat <- as.numeric(coef(cv.model, s="lambda.min"))[-1]
      }
      err <- as.numeric(y-z%*%gamma.hat)
    }
    if(model == 'logit'){
      if(p.gamma <= n/10){
        gamma.hat <- glm.fit(z, y, family = binomial(), intercept = FALSE)$coefficients
      }else{
        cv.model <- cv.glmnet(z, y, intercept=FALSE, family = "binomial")
        gamma.hat <- as.numeric(coef(cv.model, s="lambda.min"))[-1] 
      }
      err <- as.numeric(y-sigmoid(z%*%gamma.hat))
    }
    if(model == 'poisson'){
      if(p.gamma <= n/10){
        gamma.hat <- glm.fit(z, y, family = poisson(), intercept = FALSE)$coefficients
      }else{
        cv.model <- cv.glmnet(z, y, intercept=FALSE, family = "poisson")
        gamma.hat <- as.numeric(coef(cv.model, s="lambda.min"))[-1] 
      }
      err <- as.numeric(y-exp(z%*%gamma.hat))
    }
  }
  
  # derive x.mat
  if((p.gamma >= 1) & orthogonalization){
    if(p.gamma <= n/10){
      W.hat <- solve(t(z) %*% z) %*% t(z) %*% x
      # adding threshold lambda
      if(threshold){
        W.hat[abs(W.hat)<=lambda] <- 0
      }
    }else{
      W.hat <- matrix(0, p.gamma, p.beta)
      for(i in 1:p.beta){
        cv.model.xz <- cv.glmnet(z, x[,i], intercept=FALSE)
        W.hat[,i] <- as.numeric(coef(cv.model.xz, s="lambda.min"))[-1]
      }
    }
    eta <- x - z %*% W.hat
    x.mat <- eta %*% t(eta)
  }else{
    x.mat <- x %*% t(x) 
  }
  
  # construct test statistic
  err.mat <- outer(err, err, "*")
  Tn <- (sum(err.mat*x.mat) - sum(diag(err.mat*x.mat)))/n
  
  # construct p-value
  if (method == 'asymptotic'){
    if (model == 'linear'){
      tr.Sigx2.hat <- (sum(x.mat^2) - sum(diag(x.mat^2)))/(n*(n-1))
      sigma.hat <- mean(err^2)
      Tn.std <- Tn/(sigma.hat*sqrt(2*tr.Sigx2.hat))
    }else{
      err.square.mat <- err.mat^2; x.square.mat <- x.mat^2
      Rn <- (sum(err.square.mat*x.square.mat) - sum(diag(err.square.mat*x.square.mat)))/(n*(n-1))
      Tn.std <- Tn/sqrt(2*Rn)
    }
    pval <- 1-pnorm(Tn.std)
  }
  if(method == 'bootstrap'){
    if(multiplier == 'gaussian'){
      Tn.bootstrap <- sapply(1:M, function(i){
      e <- rnorm(n); e.mat <- outer(e, e, "*")
      return((sum(err.mat*x.mat*e.mat) - sum(diag(err.mat*x.mat*e.mat)))/n)
      })
    }
    if(multiplier == 'rademacher'){
      Tn.bootstrap <- sapply(1:M, function(i){
        e <- 2*(rbinom(n,1,0.5) - 0.5); e.mat <- outer(e, e, "*")
        return((sum(err.mat*x.mat*e.mat) - sum(diag(err.mat*x.mat*e.mat)))/n)
      })
    }
    pval <- (sum(Tn.bootstrap >= Tn) + 1)/(M + 1)
  }
  return(pval)
}

sigmoid <- function(t){exp(t) / (1 + exp(t))}

orthogonal_score_test <- function(v,y,test.set,screen.num = length(test.set)) {
  x <- v[,test.set]; z <- v[,-test.set]
  n <- nrow(x); p.beta <- ncol(x); p.gamma <- ncol(z)
  
  ## group testing with orthogonal score
  cv.model <- cv.glmnet(z, y, intercept=FALSE)
  gamma.hat <- as.numeric(coef(cv.model, s="lambda.min"))[-1]
  err.mat <- outer(as.numeric(y-z%*%gamma.hat), as.numeric(y-z%*%gamma.hat), "*")
  sigma.hat <- mean((y-z%*%gamma.hat)^2)
  
  # feature screening
  cor.xz <- cor(x, z)
  screen.index <- order(abs(cor.xz)[,1], decreasing=TRUE)[1:screen.num]
  W.hat <- matrix(0, p.gamma, p.beta)
  for (i in screen.index) {
    cv.model.xz <- cv.glmnet(z, x[,i], intercept=FALSE)
    W.hat[,i] <- as.numeric(coef(cv.model.xz, s="lambda.min"))[-1]
  }
  
  x.decor.mat <- (x-z%*%W.hat) %*% t(x-z%*%W.hat)
  Tn.descore <- (sum(err.mat*x.decor.mat) - sum(diag(err.mat*x.decor.mat)))/n
  tr.Sigxdecor2.hat <- (sum(x.decor.mat^2) - sum(diag(x.decor.mat^2)))/(n*(n-1))
  Tn.descore.std <- Tn.descore/(sigma.hat*sqrt(2*tr.Sigxdecor2.hat))
  pval.descore <- 1-pnorm(Tn.descore.std)
  return(pval.descore)
}

oracle_orthogonal_score_test <- function(v,eta,y,test.set) {
  x <- v[,test.set]; z <- v[,-test.set]
  n <- nrow(x); p.beta <- ncol(x); p.gamma <- ncol(z)
  
  ## group testing with L2 norm
  cv.model <- cv.glmnet(z, y, intercept=FALSE)
  gamma.hat <- as.numeric(coef(cv.model, s="lambda.min"))[-1]
  eta.mat <- eta %*% t(eta)
  err.mat <- outer(as.numeric(y-z%*%gamma.hat), as.numeric(y-z%*%gamma.hat), "*")
  Tn.oracle <- (sum(err.mat*eta.mat) - sum(diag(err.mat*eta.mat)))/n
  tr.Sigeta2.hat <- (sum(eta.mat^2) - sum(diag(eta.mat^2)))/(n*(n-1))
  sigma.hat <- mean((y-z%*%gamma.hat)^2)
  Tn.oracle.std <- Tn.oracle/(sigma.hat*sqrt(2*tr.Sigeta2.hat))
  pval.oracle <- 1-pnorm(Tn.oracle.std)
  
  return(pval.oracle)
}



ST <- function(X.f, Y.f, sub.size=nrow(X.f)*0.3,ncores = 1) {
  n <- dim(X.f)[1]
  p <- dim(X.f)[2]
  
  n1 <- sub.size
  n0 <- n-floor(n1)
  S1 <- sample(1:n, floor(n1), replace=FALSE)
  X.sub <- X.f[S1,]
  Y.sub <- Y.f[S1]
  cvfit <- cv.glmnet(X.sub, Y.sub, intercept=FALSE)
  cf <- as.numeric(coef(cvfit, s="lambda.min"))[-1]
  # cvfit <- cv.ncvreg(X.sub, Y.sub, penalty=penalty, nfolds=nfolds)
  # model.cvfit <- ncvfit(X.sub, Y.sub, penalty=penalty, lambda=cvfit$lambda.min)
  # cf <- as.numeric(model.cvfit$beta)
  set1 <- (1:p)[abs(cf)>0]
  resi <- Y.sub-X.sub%*%cf
  beta.m <- t(standardize(X.sub[,-set1]))%*%resi
  screen.set <- sort(order(abs(beta.m),decreasing=TRUE)[1:(n0-1-length(set1))])
  a <- (1:p)[-set1]
  screen.set <- union(a[screen.set],set1)
  X <- X.f[-S1,screen.set]
  Y <- Y.f[-S1]
  
  score.nodewiselasso = getFromNamespace("score.nodewiselasso", "hdi")
  node <- score.nodewiselasso(X, wantTheta=TRUE, verbose=FALSE, lambdaseq="quantile",
                              parallel=TRUE, ncores=ncores, oldschool = FALSE, lambdatuningfactor = 1)
  Theta <- node$out
  return(list(X=X, Y=Y, n0=n0, screen.set=screen.set, Theta=Theta))
}

ST_res <- function(test.set, X, Y, n0, screen.set, Theta, M=500) {
  Gram <- t(X)%*%X/n0
  sreg <- scalreg(X,Y)
  beta.hat <- sreg$coefficients
  sigma.sq <- sum((Y-X%*%beta.hat)^2)/(n0-sum(abs(beta.hat)>0))
  test.set.i <- intersect(screen.set,test.set)
  index <- screen.set%in%test.set.i
  
  Omega <- diag(Theta%*%Gram%*%t(Theta))*sigma.sq
  beta.db <- beta.hat+Theta%*%t(X)%*%(Y-X%*%beta.hat)/n0
  margin.st <- sqrt(n0)*abs(beta.db[index])/sqrt(Omega[index])
  margin.nst <- sqrt(n0)*abs(beta.db[index])
  stat.st <- max(margin.st)
  stat.nst <- max(margin.nst)
  
  stat.boot.st <- stat.boot.nst <- rep(NA,M)
  for (i in 1:M) {
    e <- rnorm(n0)
    xi.boot <- Theta[index,]%*%t(X)%*%e*sqrt(sigma.sq)/sqrt(n0)
    stat.boot.nst[i] <- max(abs(xi.boot))
    stat.boot.st[i] <- max(abs(xi.boot/sqrt(Omega[index])))
  }
  
  # if (stat.nst>quantile(stat.boot.nst,1-alpha)) rej.nst <- "reject" else rej.nst <- "fail to reject"
  # if (stat.st>quantile(stat.boot.st,1-alpha)) rej.st <- "rejct" else rej.st <- "fail to reject"
  # result <- list(stat.nst, rej.nst, stat.st, rej.st)
  # names(result) <- c("non-studentized test","non-studentized test","studentized test","studentized test")
  # rej.nst <- stat.nst > quantile(stat.boot.nst, 1-alpha)
  pval.nst <- sum(stat.boot.nst >= stat.nst)/M
  # rej.st <- stat.st > quantile(stat.boot.st, 1-alpha)
  pval.st <- sum(stat.boot.st >= stat.st)/M
  # result <- c(rej.nst, rej.st)
  # return(result)
  pvalue <- matrix(c(pval.nst, pval.st),nrow = 1,ncol = 2)
  colnames(pvalue) <- c("NST","ST")
  return(pvalue)
}

max_norm_test <- function(v,y,test.set,M=500,sub.size = nrow(v)*0.3,ncores = 1){
  ST.para <- ST(v, y, sub.size=sub.size, ncores = 1)
  return(ST_res(test.set, ST.para$X, ST.para$Y, ST.para$n0, 
                ST.para$screen.set, ST.para$Theta, M=500))
}

CGZ_test <- function(v, y, test.set = NULL, model = 'linear', orthogonalization = FALSE, scale = FALSE){
  if(scale) {
    v <- scale(v); y <- scale(y)
  }
  if(is.null(test.set)) {
    test.set <- 1:ncol(v)
  }
  x <- v[,test.set]; z <- v[,-test.set]
  n <- nrow(x); p.beta <- ncol(x); p.gamma <- ncol(z)
  
  ## group testing with L2 norm
  # derive epsilon
  if(ncol(z) == 0){
    err <- as.numeric(y)
  }else{
    if(model == 'linear'){
      if(p.gamma <= n/10){
        gamma.hat <- lm.fit(as.matrix(z),y)$coefficients
      }else{
        cv.model <- cv.glmnet(z, y, intercept=FALSE)
        gamma.hat <- as.numeric(coef(cv.model, s="lambda.min"))[-1]
      }
      err <- as.numeric(y-z%*%gamma.hat)
    }else{
      print("CGZ test is only applicable to linear models!")
      break
    }
  }
  
  # derive x.mat
  if((p.gamma >= 1) & orthogonalization){
    if(p.gamma <= n/10){
      W.hat <- solve(t(z) %*% z) %*% t(z) %*% x
    }else{
      W.hat <- matrix(0, p.gamma, p.beta)
      for(i in 1:p.beta){
        cv.model.xz <- cv.glmnet(z, x[,i], intercept=FALSE)
        W.hat[,i] <- as.numeric(coef(cv.model.xz, s="lambda.min"))[-1]
      }
    }
    eta <- x - z %*% W.hat
    x.mat <- eta %*% t(eta)
  }else{
    x.mat <- x %*% t(x) 
  }
  
  # construct test statistic
  x.centered <- scale(x,scale = FALSE); x.square.mat <- diag(diag(x.mat)) %*% matrix(1,n,n)
  delta_x.mat <- x.centered %*% t(x.centered) + (x.square.mat + t(x.square.mat) - 2 * x.mat)/(2*n)
  err.centered <- scale(err,scale = FALSE)
  err.mat <- outer(err, err, "*"); err.square.mat <- diag(err^2) %*% matrix(1,n,n)
  delta_y.mat <- (err.centered %*% t(err.centered)) + (err.square.mat + t(err.square.mat) - 2 * err.mat)/(2*n)
  
  
  Tn <- ((1 - 2/n)^{-2})*(choose(n,2)^{-1}) * (sum(delta_x.mat * delta_y.mat) - sum(diag(delta_x.mat * delta_y.mat)))/2
  
  # construct p-value
  
  tr.Sigx2.hat <- (sum(x.mat^2) - sum(diag(x.mat^2)))/(n*(n-1))
  sigma.hat <- mean(err^2)
  Tn.std <- n*Tn/(sigma.hat*sqrt(2*tr.Sigx2.hat))
  
  pval <- 1-pnorm(Tn.std)
  return(pval)
}