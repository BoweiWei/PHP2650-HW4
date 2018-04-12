X <- as.matrix(read.csv("X.csv", header = F, sep = ","))
y <- as.matrix(read.csv("y.csv", header = F, sep = ","))

library(glmnet)


soft_thresholding <- function(x,a){
  ## This could be done more efficiently using vector multiplication
  ## See the forumula in slides
  ##  sign(x)*pmax(abs(x) - a, 0)
  result <- numeric(length(x))
  result[which(x > a)] <- x[which(x > a)] - a
  result[which(x < -a)] <- x[which(x < -a)] + a
  return(result)
}



lasso_kkt_check <- function(X,y,beta,lambda, tol=1e-3){
  ## check convergence 
  beta <- as.matrix(beta); X <- as.matrix(X)
  ## Assuming no intercepts 
  G <- t(X)%*%(y-X%*%beta)/length(y)
  ix <- which(beta == 0 )
  iy <- which(beta != 0)
  if (any(abs(G[ix]) > (lambda + tol) )) { return(pass=0) }
  if (any(abs( G[iy] - lambda*sign(beta[iy] ) ) > tol)) { return(pass=0) }  
  return(pass=1)
}



lasso.cd <- function(X,y,beta,lambda,tol=1e-8,maxiter=100000,quiet=FALSE){
  # note that the LS part  in this function is the one in slides divided by length(y) = n 
  ## Or equivalently  lambda here = n * lambda in class
  beta <- as.matrix(beta); X <- as.matrix(X)
  obj <- numeric(length=(maxiter+1))
  betalist <- list(length=(maxiter+1))
  betalist[[1]] <- beta
  N = length(y)
  
  for (j in 1:maxiter){
    for (k in 1:length(beta)){
      r <- y - X[,-k]%*%beta[-k]
      beta[k] <- (1/norm(as.matrix(X[,k]),"F")^2)*soft_thresholding(t(r)%*%X[,k],length(y)*lambda)
      #beta[k] <- (1/N)*soft_thresholding(t(r)%*%X[,k],length(y)*lambda)
    }
    betalist[[(j+1)]] <- beta
    obj[j] <- (1/2)*(1/length(y))*norm(y - X%*%beta,"F")^2 + lambda*sum(abs(beta))
    if (norm(betalist[[j]] - beta,"F") < tol) { break }
  } 
  check <- lasso_kkt_check(X,y,beta,lambda,tol=tol) 
  
  if (quiet==FALSE){
    if (check==1) {
      cat(noquote("Minimum obtained.\n"))
    }
    else { cat(noquote("Minimum not obtained.\n")) } 
  }
  return(list(obj=obj[1:j],beta=beta)) 
}


elastic.cd <- function(X,y,beta,lambda,tol=1e-10,maxiter=1000000,quiet=FALSE, alpha = 0.95){
  # note that the LS part  in this function is the one in slides divided by length(y) = n 
  ## Or equivalently  lambda here = n * lambda in class
  beta <- as.matrix(beta); X <- as.matrix(X)
  obj <- numeric(length=(maxiter+1))
  betalist <- list(length=(maxiter+1))
  betalist[[1]] <- beta
  
  for (j in 1:maxiter){
    for (k in 1:length(beta)){
      r <- y - X[,-k]%*%beta[-k]
      beta[k] <- (1/norm(as.matrix(X[,k]),"F")^2)*soft_thresholding(t(r)%*%X[,k],length(y)*lambda*alpha)/ (1 + lambda*(1-alpha))
    }
    betalist[[(j+1)]] <- beta
    obj[j] <- (1/2)*(1/length(y))*norm(y - X%*%beta,"F")^2 + lambda*sum(abs(beta)) 
    if (norm(betalist[[j]] - beta,"F") < tol) { break }
  } 
  check <- lasso_kkt_check(X,y,beta,lambda, tol=tol) 
  
  if (quiet==FALSE){
    if (check==1) {
      cat(noquote("Minimum obtained.\n"))
    }
    else { cat(noquote("Minimum not obtained.\n")) } 
  }
  return(list(obj=obj[1:j],beta=beta)) 
}


lasso_cd <- lasso.cd(X, y, rep(0, 400), 1)
idx_lasso_cd = which(lasso_cd$beta !=0)

glm_lasso = glmnet(X,y,lambda = 1, alpha = 1)
idx_glm_lasso = which(glm_lasso$beta !=0)

ela_cd_lasso <- elastic.cd(X, y, rep(0, 400), 1, alpha =1)
idx_ela_cd_lasso = which(ela_cd_lasso$beta !=0)


idx_glm_lasso
idx_lasso_cd
idx_ela_cd_lasso

glm_lasso$beta[idx_glm_lasso]
lasso_cd$beta[idx_lasso_cd]
ela_cd_lasso$beta[idx_ela_cd_lasso]



glmet_ela = glmnet(X,y,lambda = 1, alpha = 0.95, intercept = F, standardize = F)
idx_glmnet_ela = which(glmet_ela$beta !=0)

ela_cd <- elastic.cd(X, y, rep(0, 400), 1)
idx_ela_cd = which(ela_cd$beta !=0)



idx_glmnet_ela
idx_ela_cd
glmet_ela$beta[idx_glmnet_ela]
ela_cd$beta[idx_ela_cd]