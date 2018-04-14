X <- as.matrix(read.csv("X.csv", header = F, sep = ","))
y <- as.matrix(read.csv("y.csv", header = F, sep = ","))

scale_X <- scale(X, center = T)
scale_y <- scale(y, center = T)

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


lasso.cd <- function(X,y,beta = rep(0, 400),lambda,tol=1e-6,maxiter=10000,quiet=FALSE){
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


elastic.cd <- function(X,y,beta,lambda,tol=1e-6,maxiter=10000,quiet=FALSE, alpha = 0.95){
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
    obj[j] <- (1/2)*(1/length(y))*norm(y - X%*%beta,"F")^2 + lambda*(sum(abs(beta))*alpha + 1/2*(1-alpha)*norm(beta,"F")^2) 
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







test_lasso.cd <- function(X,y,beta= rep(0, 400),lambdas = c(0.01,0.1,1,10),
                          tol=1e-8,maxiter=10000,quiet=FALSE, alpha = 1, 
                          intercept = F, standardize = F, thresh = 1e-8){
  beta_output <- data.frame(matrix( rep(NA, 2*400*length(lambdas)), nrow = 400, ncol = 8))
  for (i in 1:length(lambdas)){
    beta_output[, 2*i-1] <- lasso.cd(X,y,beta,lambda = lambdas[i],tol,maxiter,quiet)$beta
    beta_output[, 2*i] <- as.matrix(glmnet(X,y,lambda = lambdas[i], alpha = alpha, intercept = intercept, 
                                           standardize = standardize, thresh = thresh )$beta)
  }
  colnames(beta_output) <- paste(rep(c("cd","glmnet"), 4),rep( c(0.01,0.1,1,10),each=2), sep = " ")
  return(beta_output)
}


###genereate a table to compare results from lasso.cd and glmnet
lasso_cd_glmnet <- test_lasso.cd(scale_X, scale_y,lambdas = c(0.01,0.1,1,10))
write.csv(lasso_cd_glmnet, "lasso_vs_glmnet.csv")



test_elastic.cd <- function(X,y,beta= rep(0, 400),lambdas = c(0.01,0.1,1,10),
                            tol=1e-8,maxiter=10000,quiet=FALSE, alpha = 0.95, 
                            intercept = F, standardize = F, thresh = 1e-8){
  beta_output <- data.frame(matrix( rep(NA, 2*400*length(lambdas)), nrow = 400, ncol = 8))
  for (i in 1:length(lambdas)){
    beta_output[, 2*i-1] <- elastic.cd(X,y,beta,lambda = lambdas[i],tol,maxiter,quiet, alpha)$beta
    beta_output[, 2*i] <- as.matrix(glmnet(X,y,lambda = lambdas[i], alpha = alpha, intercept = intercept, 
                                           standardize = standardize, thresh = thresh )$beta)
  }
  colnames(beta_output) <- paste(rep(c("elastic","glmnet"), 4),rep( c(0.01,0.1,1,10),each=2), sep = " ")
  return(beta_output)
}


###genereate a table to compare results from elastic.cd and glmnet
elastic_cd_glmnet <- test_elastic.cd(scale_X, scale_y,lambdas = c(0.01,0.1,1,10))
write.csv(elastic_cd_glmnet, "elastic_vs_glmnet.csv")





###another way to check by comparing the non-zero coefficients
glmet_ela = glmnet(scale_X,scale_y,lambda = 0.1, alpha = 0.95, intercept = F, standardize = F, thresh = 1e-6)
idx_glmnet_ela = which(glmet_ela$beta !=0)

ela_cd <- elastic.cd(scale_X,scale_y, rep(0, 400), 0.1, alpha = 0.95)
idx_ela_cd = which(ela_cd$beta !=0)



idx_glmnet_ela
idx_ela_cd
glmet_ela$beta[idx_glmnet_ela]
ela_cd$beta[idx_ela_cd]


lasso_cd <- lasso.cd(scale_X, scale_y,lambda = 0.1, rep(0, 400))
idx_lasso_cd = which(lasso_cd$beta !=0)

glm_lasso = glmnet(scale_X, scale_y,lambda = 0.1, alpha = 1, intercept = F, standardize = F, thresh = 1e-6)
idx_glm_lasso = which(glm_lasso$beta !=0)

ela_cd_lasso <- elastic.cd(scale_X, scale_y, rep(0, 400), 0.1, alpha =1)
idx_ela_cd_lasso = which(ela_cd_lasso$beta !=0)


idx_glm_lasso
idx_lasso_cd
idx_ela_cd_lasso

glm_lasso$beta[idx_glm_lasso]
lasso_cd$beta[idx_lasso_cd]
ela_cd_lasso$beta[idx_ela_cd_lasso]
