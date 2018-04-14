authors <- function(){
  c("Bernard Chu", "Bowei Wei", "Jiarou Quan", "Ning Zhang")
}


X <- as.matrix(read.csv("X.csv", header = F, sep = ","))
y <- as.matrix(read.csv("y.csv", header = F, sep = ","))

##standardize data to have features with same scale, center data to make intercept 0
scale_X <- scale(X, center = T)
scale_y <- scale(y, center = T)


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


output_elastic.cd <- function(X,y,beta= rep(0, 400),lambdas = c(0.01,0.1,1,10),
                            tol=1e-8,maxiter=10000,quiet=FALSE, alpha = 0.95){
  ncolumn = length(lambdas)
  beta_output <- data.frame(matrix( rep(NA, 400*ncolumn), nrow = 400, ncol =ncolumn ))
  for (i in 1:length(lambdas)){
    beta_output[, i] <- elastic.cd(X,y,beta,lambda = lambdas[i],tol,maxiter,quiet, alpha)$beta
  }
  colnames(beta_output) <- paste("lambda", c(0.01,0.1,1,10), sep = "")
  return(beta_output)
}


beta_elastic_0.95 <- output_elastic.cd(scale_X, scale_y,lambdas = c(0.01,0.1,1,10))
beta_elastic_1 <- output_elastic.cd(scale_X, scale_y,lambdas = c(0.01,0.1,1,10), alpha = 1)

write.csv(beta_elastic_0.95, "b1.csv")
write.csv(beta_elastic_1, "b2.csv")


