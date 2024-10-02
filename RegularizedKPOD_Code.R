##############################################
########  main functions & examples    #######
##############################################

#packages--------------
library(MASS)
library(abind)
library(RSKC)
library(Rcpp)
library (RcppArmadillo) 
library(Rfast)
library(deldir)
library(kpodclustr)

#define functions--------------

#K-means++
kmeans_pp <- function(x, k, iter.max = 10, nstart = 1, ...){
  n <- nrow(x) # number of data points
  centers <- numeric(k) # IDs of centers
  distances <- matrix(numeric(n * (k - 1)), ncol = k - 1) # distances[i, j]: The distance between x[i,] and x[centers[j],]
  res.best <- list(tot.withinss = Inf) # the best result among <nstart> iterations
  for (rep in 1:nstart) {
    ini_centers <- qpp_arma(x, k, r = 2)$centers
    ## Perform k-means with the obtained centers
    res <- kmeans(x, ini_centers, iter.max = iter.max, nstart = 1, ...)
    res$initial.centers <- x[centers, ]
    ## Store the best result
    if (res$tot.withinss < res.best$tot.withinss) {
      res.best <- res
    }
  }
  res.best
}

#kpod
fun_kpod<-function(x, k, nstart = 1, iter.max = 100, eps = 1e-5, init_method = c("random","comp")){
  n <- nrow(x)
  p <- ncol(x)
  loss_vec <- rep(0,iter.max)
  
  miss_arrind <- which(is.na(x) == TRUE, arr.ind = TRUE)
  miss_ind <- which(is.na(x) == TRUE, arr.ind = FALSE)
  omg_mat <- matrix(1, nrow = n, ncol = p)
  omg_mat[miss_ind] <- 0
  omgc_mat <- matrix(0, nrow = n, ncol = p)
  omgc_mat[miss_ind] <- 1
  
  mu_vec <- colMeans(x, na.rm = TRUE)
  
  if(init_method[1] == "random"){
    cc_indicator <- 0
  }
  
  if(init_method[1] == "comp"){
    complete_obs <- which(Rfast::rowsums(omg_mat) == p)
    miss_obs <- which(Rfast::rowMins(omg_mat, value = TRUE) == 0)
    if(length(complete_obs) < k){
      print("There are no enough complete cases. We will use random starts")
      cc_indicator <- 0
    }else{
      x_comp <- x[-miss_obs,]#complete data
      n_comp <- nrow(x_comp)
      x_tmp <- x
      x_tmp[miss_ind] <- 0
      mu_comp <- colMeans(x_comp, na.rm = TRUE)
      cc_indicator <- 1
    }
  }
  
  umat <- matrix(0, n, k)
  cmat <- matrix(0, k, p)
  
  can_loss <- Inf
  for(t in 1:nstart){
    tmp_j <- iter.max
    #Initialization: cc_indicator=1 comp;  0 random
    #------------------------------------------------------------------------
    if(cc_indicator == 1){
      #Comp initialization
      ini_centers <- kmeans_pp(x_comp, k, iter.max = 100, nstart = 100)$centers
      cmat <- ini_centers
      
      #Compute an initial membership matrix
      dist_mat <- proxy::dist(x, cmat)
      min_val <- Rfast::rowMins(dist_mat, value = TRUE)^2
      cluster <- Rfast::rowMins(dist_mat)
      umat <- matrix(0, n, k)
      for(j in 1:k){
        umat[which(cluster==j),j] <- 1
      }
    }else{
      #Random initialization
      ini_index <- sample(1:n, k)
      x_tmp <- x
      mu_mat <- matrix(mu_vec, nrow = n, ncol = p, byrow = TRUE)
      x_tmp[which(is.na(x))] <- mu_mat[which(is.na(x))] 
      cmat <- x_tmp[ini_index,]
      if(any(duplicated(cmat)) == TRUE){
        cmat <- cmat + matrix(rnorm(k*p,0,0.01),k,p)
      }
      #Compute an initial membership matrix
      dist_mat <- proxy::dist(x, cmat)
      min_val <- Rfast::rowMins(dist_mat, value = TRUE)^2
      cluster <- Rfast::rowMins(dist_mat)
      umat <- matrix(0, n, k)
      for(j in 1:k){
        umat[which(cluster==j),j] <- 1
      }  
    }
    #------------------------------------------------------------------------
    
    #Fill-in for x
    ab_mat <- umat%*%cmat
    x_tmp <- x_tmp*omg_mat + ab_mat*omgc_mat
    
    #Main loop
    #------------------------------------------------------------------------
    for(r in 1:iter.max){
      #Update centroids and menbership matrix
      try_km <- try(res_km <- kmeans(x = x_tmp, centers = cmat, iter.max = iter.max))
      if(class(try_km) == "try-error"){
        break
      }
      cmat <- res_km$centers
      cluster <- res_km$cluster
      umat <- matrix(0, n, k)
      for(j in 1:k){
        umat[which(cluster==j),j] <- 1
      }
      
      #Fill in
      ab_mat <- umat%*%cmat
      x_tmp <- x_tmp*omg_mat + ab_mat*omgc_mat
      
      #Record loss
      loss_vec[r] <- sum( (x[-miss_ind] - ab_mat[-miss_ind])^2 )
      if(r > 1){
        if( (loss_vec[r-1] - loss_vec[r])/loss_vec[r-1] < eps){
          tmp_j <- r
          tmp_res <- list(centers = cmat, cluster = cluster, loss = loss_vec[1:tmp_j])
          break
        }
      }
    }
    #------------------------------------------------------------------------
    if(tmp_j == iter.max){
      tmp_res <- list(centers = cmat, cluster = cluster, loss = loss_vec[1:tmp_j])
    }
    if(tmp_res$loss[tmp_j] < can_loss){
      can_loss <- tmp_res$loss[tmp_j]
      res <- tmp_res
    }
  }#END for t
  return(res)
}

#rkpod
fun_rkpod_g<-function(x, k, lambda=0.1, weight_vec=weight_vec, nstart = 1, iter.max = 100, eps = 1e-5, init_method = c("random","comp")){
  n <- nrow(x)
  p <- ncol(x)
  loss_vec <- NULL
  
  miss_arrind <- which(is.na(x) == TRUE, arr.ind = TRUE)
  miss_ind <- which(is.na(x) == TRUE, arr.ind = FALSE)
  omg_mat <- matrix(1, nrow = n, ncol = p)
  omg_mat[miss_ind] <- 0
  omgc_mat <- matrix(0, nrow = n, ncol = p)
  omgc_mat[miss_ind] <- 1
  
  mu_vec <- colMeans(x, na.rm = TRUE)
  
  if(init_method[1] == "random"){
    cc_indicator <- 0
  }
  
  if(init_method[1] == "comp"){
    complete_obs <- which(Rfast::rowsums(omg_mat) == p)
    miss_obs <- which(Rfast::rowMins(omg_mat, value = TRUE) == 0)
    if(length(complete_obs) < k){
      print("There are no enough complete cases. We will use random starts")
      cc_indicator <- 0
    }else{
      x_comp <- x[-miss_obs,]#complete data
      n_comp <- nrow(x_comp)
      x_tmp <- x
      x_tmp[miss_ind] <- 0
      mu_comp <- colMeans(x_comp, na.rm = TRUE)
      cc_indicator <- 1
    }
  }
  
  umat <- matrix(0, n, k)
  cmat <- matrix(0, k, p)
  
  can_loss <- Inf
  for(t in 1:nstart){
    tmp_j <- iter.max
    #Initialization: cc_indicator=1 comp;  0 random
    #------------------------------------------------------------------------
    if(cc_indicator == 1){
      #Comp initialization
      ini_centers <- kmeans_pp(x_comp, k, iter.max = 100, nstart = 100)$centers
      cmat <- ini_centers
      
      #Compute an initial membership matrix
      dist_mat <- proxy::dist(x, cmat)
      min_val <- Rfast::rowMins(dist_mat, value = TRUE)^2
      cluster <- Rfast::rowMins(dist_mat)
      umat <- matrix(0, n, k)
      for(kk in 1:k){
        umat[which(cluster==kk),kk] <- 1
      }
    }else{
      #Random initialization
      ini_index <- sample(1:n, k)
      x_tmp <- x
      mu_mat <- matrix(mu_vec, nrow = n, ncol = p, byrow = TRUE)
      x_tmp[which(is.na(x))] <- mu_mat[which(is.na(x))] 
      cmat <- x_tmp[ini_index,]
      if(any(duplicated(cmat)) == TRUE){
        cmat <- cmat + matrix(rnorm(k*p,0,0.01),k,p)
      }
      #Compute an initial membership matrix
      dist_mat <- proxy::dist(x, cmat)
      min_val <- Rfast::rowMins(dist_mat, value = TRUE)^2
      cluster <- Rfast::rowMins(dist_mat)
      umat <- matrix(0, n, k)
      for(kk in 1:k){
        umat[which(cluster==kk),kk] <- 1
      }  
    }
    #------------------------------------------------------------------------
    
    #Fill in x
    ab_mat <- umat%*%cmat
    x_tmp <- x_tmp*omg_mat + ab_mat*omgc_mat
    
    #Main loop
    #------------------------------------------------------------------------
    for(r in 1:iter.max){
      #Update centroids and menbership matrix
      try_rkm <- try(res_rkm <- solve_rkm_g_cpp(x = x_tmp, k=k, weight_vec = weight_vec, lambda=lambda, initial_centers = cmat, itermax = iter.max, eps=eps))
      if(class(try_rkm) == "try-error"){
        break
      }
      cmat <- res_rkm$centers
      cluster <- res_rkm$cluster
      umat <- matrix(0, n, k)
      for(kk in 1:k){
        umat[which(cluster==kk),kk] <- 1
      }
      
      #Fill in x
      ab_mat <- umat%*%cmat
      x_tmp <- x_tmp*omg_mat + ab_mat*omgc_mat
      
      #Record loss
      wcss <- sum( (x[-miss_ind] - ab_mat[-miss_ind])^2 )
      penalty <- n*lambda*sum(weight_vec*apply(cmat,2,function(z) sqrt(sum(z^2)) ))
      loss_vec[r] <- wcss + penalty
      if(r > 1){
        if( (loss_vec[r-1] - loss_vec[r])/loss_vec[r-1] < eps){
          tmp_j <- r
          tmp_res <- list(centers = cmat, cluster = cluster, loss = loss_vec[1:tmp_j],wcss=wcss)
          break
        }
      }
    }
    #------------------------------------------------------------------------
    if(tmp_j == iter.max){
      tmp_res <- list(centers = cmat, cluster = cluster, loss = loss_vec[1:tmp_j],wcss=wcss)
    }
    if(tmp_res$loss[tmp_j] < can_loss){
      can_loss <- tmp_res$loss[tmp_j]
      res <- tmp_res
    }
  }#END for t
  return(res)
}

fun_rkpod_L0 <- function(x, k, lambda=0.1, nstart = 1, iter.max = 100, eps = 1e-5, init_method = c("random","comp")){
  n <- nrow(x)
  p <- ncol(x)
  loss_vec <- NULL
  
  miss_arrind <- which(is.na(x) == TRUE, arr.ind = TRUE)
  miss_ind <- which(is.na(x) == TRUE, arr.ind = FALSE)
  omg_mat <- matrix(1, nrow = n, ncol = p)
  omg_mat[miss_ind] <- 0
  omgc_mat <- matrix(0, nrow = n, ncol = p)
  omgc_mat[miss_ind] <- 1
  
  mu_vec <- colMeans(x, na.rm = TRUE)
  
  if(init_method[1] == "random"){
    cc_indicator <- 0
  }
  
  if(init_method[1] == "comp"){
    complete_obs <- which(Rfast::rowsums(omg_mat) == p)
    miss_obs <- which(Rfast::rowMins(omg_mat, value = TRUE) == 0)
    if(length(complete_obs) < k){
      print("There are no enough complete cases. We will use random starts")
      cc_indicator <- 0
    }else{
      x_comp <- x[-miss_obs,]#complete data
      n_comp <- nrow(x_comp)
      x_tmp <- x
      x_tmp[miss_ind] <- 0
      mu_comp <- colMeans(x_comp, na.rm = TRUE)
      cc_indicator <- 1
    }
  }
  
  umat <- matrix(0, n, k)
  cmat <- matrix(0, k, p)
  
  can_loss <- Inf
  for(t in 1:nstart){
    tmp_j <- iter.max
    #Initialization: cc_indicator=1 comp;  0 random
    #------------------------------------------------------------------------
    if(cc_indicator == 1){
      #Comp initialization
      ini_centers <- kmeans_pp(x_comp, k, iter.max = 100, nstart = 100)$centers
      cmat <- ini_centers
      
      #Compute an initial membership matrix
      dist_mat <- proxy::dist(x, cmat)
      min_val <- Rfast::rowMins(dist_mat, value = TRUE)^2
      cluster <- Rfast::rowMins(dist_mat)
      umat <- matrix(0, n, k)
      for(kk in 1:k){
        umat[which(cluster==kk),kk] <- 1
      }
    }else{
      #Random initialization
      ini_index <- sample(1:n, k)
      x_tmp <- x
      mu_mat <- matrix(mu_vec, nrow = n, ncol = p, byrow = TRUE)
      x_tmp[which(is.na(x))] <- mu_mat[which(is.na(x))] 
      cmat <- x_tmp[ini_index,]
      if(any(duplicated(cmat)) == TRUE){
        cmat <- cmat + matrix(rnorm(k*p,0,0.01),k,p)
      }
      #Compute an initial membership matrix
      dist_mat <- proxy::dist(x, cmat)
      min_val <- Rfast::rowMins(dist_mat, value = TRUE)^2
      cluster <- Rfast::rowMins(dist_mat)
      umat <- matrix(0, n, k)
      for(kk in 1:k){
        umat[which(cluster==kk),kk] <- 1
      }  
    }
    #------------------------------------------------------------------------
    
    #Fill in x
    ab_mat <- umat%*%cmat
    x_tmp <- x_tmp*omg_mat + ab_mat*omgc_mat
    
    #Main loop
    #------------------------------------------------------------------------
    for(r in 1:iter.max){
      #Update centroids and menbership matrix
      try_rkm <- try(res_rkm <- solve_rkm_L0_cpp(x = x_tmp, k=k , lambda=lambda, initial_centers = cmat, itermax = iter.max, eps = eps))
      if(class(try_rkm) == "try-error"){
        break
      }
      cmat <- res_rkm$centers
      cluster <- res_rkm$cluster
      umat <- matrix(0, n, k)
      for(kk in 1:k){
        umat[which(cluster==kk),kk] <- 1
      }
      
      #Fill in x
      ab_mat <- umat%*%cmat
      x_tmp <- x_tmp*omg_mat + ab_mat*omgc_mat
      
      #Record loss
      wcss <- sum( (x[-miss_ind] - ab_mat[-miss_ind])^2 )
      penalty <- n*lambda*sum(apply(cmat,2,function(z) sqrt(sum(z^2)) )>0)
      loss_vec[r] <- wcss + penalty
      if(r > 1){
        if((loss_vec[r-1] - loss_vec[r])/loss_vec[r-1] < eps){
          tmp_j <- r
          tmp_res <- list(centers = cmat, cluster = cluster, loss = loss_vec[1:tmp_j],wcss=wcss)
          break
        }
      }
    }
    #------------------------------------------------------------------------
    if(tmp_j == iter.max){
      tmp_res <- list(centers = cmat, cluster = cluster, loss = loss_vec[1:tmp_j],wcss=wcss)
    }
    if(tmp_res$loss[tmp_j] < can_loss){
      can_loss <- tmp_res$loss[tmp_j]
      res <- tmp_res
    }
  }#END for t
  return(res)
}


#Example 1-----------------

##generate data for example 1
EG1_data<-function(n){
  
  zlab <- rbinom(n = n, size = 1, prob = 1/2)
  x <- matrix(0, nrow = n, ncol = 2)
  n0 <- sum(1-zlab)
  n1 <- n - n0
  
  a <- 2
  
  x[which(zlab == 0), 2] <- rnorm(n0, mean = -a, sd = 1)
  x[which(zlab == 1), 2] <- rnorm(n1, mean = a, sd = 1)
  x[which(zlab == 0), 1] <- rnorm(n0, mean = 0, sd = 2)
  x[which(zlab == 1), 1] <- rnorm(n1, mean = 0, sd = 2)
  
  miss_mat <- cbind(rbinom(n = n, size = 1, prob = 2/3), rbinom(n = n, size = 1, prob = 1/3))
  tmp_index <- which(rowSums(miss_mat) == 0)
  miss_index <- which(rowSums(miss_mat) == 2)
  
  x_miss <- x
  x_miss[which(miss_mat == 1)] <- NA
  x_miss <- x_miss[-miss_index,]
  x_comp <- x[tmp_index,]
  
  return(list(Orig=x, Missing=x_miss, CompleteCase=x_comp, truth=zlab[-miss_index]+1, truth_orig=zlab+1))
}

k=2
p=2
d=1
n=10000

##different methods for example 1
set.seed(123)
data=EG1_data(n=10000)
#kmeans
kmres=kmeans(data$Orig,centers=k, nstart = 100, iter.max = 100)
#kpod
kpodres=fun_kpod(x=data$Missing,k=k,nstart = 100, init_method = "random")
#proposed method (L0)
rkpodres_0_rand=fun_rkpod_L0(x=data$Missing,k=k,lambda = 1, nstart=100, init_method = "random")

##draw figure 1
par(mar=c(3,3,1,1))
tesselation_km <- deldir(kmres$centers, rw = c(-7,7,-7,7))
tesselation_kpod <- deldir(kpodres$centers, rw = c(-7,7,-7,7))
tesselation_rkpod_0 <- deldir(rkpodres_0_rand$centers, rw = c(-7,7,-7,7))
#kmeans
plot(data$Orig, pch = 16, col = grey(0.8), asp = 1, cex = 0.2, xlim = c(-7,7), ylim = c(-7,7))
plot(tesselation_km, add=TRUE, wlines="tess", pch = 17, lty = 1, cmpnt_col = c(NULL,"black"), showpoints = FALSE, cex = 1, lwd = 2.5)
points(kmres$centers, pch = "+", cex = 1.5, col = "black", lwd = 1.5)
#kpod
plot(data$Orig, pch = 16, col = grey(0.8), asp = 1, cex = 0.2, xlim = c(-7,7), ylim = c(-7,7))
plot(tesselation_kpod, add=TRUE, wlines="tess", pch = 17, lty = 3, cmpnt_col = c(NULL,"red"), showpoints = FALSE, cex = 1, lwd = 2)
points(kpodres$centers, pch = 2, cex = 1, col = "red", lwd = 1.5)
#proposed method (L0)
plot(data$Orig, pch = 16, col = grey(0.8), asp = 1, cex = 0.2, xlim = c(-7,7), ylim = c(-7,7))
plot(tesselation_rkpod_0, add=TRUE, wlines="tess", pch = 17, lty = 4, cmpnt_col = c(NULL,"green4"), showpoints = FALSE, cex = 1, lwd = 2)
points(rkpodres_0_rand$centers, pch = 1, cex = 1.5, col = "green4", lwd = 1.5)





#Example 2-----------------

##generate data for example 2
EG2_data<-function(n,p=2,k=2,mustar=rbind(c(0,2),c(0,-2)),sigma=c(1,1),missing=c(1/3,2/3)){
  x=matrix(0,n,p)
  
  #true label
  truth=sample(1:k,size=n,replace = TRUE)
  #full data
  for (kk in 1:k){
    x[which(truth==kk),]=mvrnorm(n=sum(truth==kk),mu=mustar[kk,],Sigma = diag(sigma,p,p))
  }
  #missing position
  ismiss=matrix(0,n,p)
  for (j in 1:p){
    ismiss[,j]=rbinom(n,size=1,prob=missing[j])
  }
  #result
  Orig=x
  Missing=x
  Missing[ismiss==1]=NA
  CompleteCase=x[which(rowSums(ismiss)==0),]
  return(list(Orig=Orig,Missing=Missing,CompleteCase=CompleteCase,truth=truth))
  
}

n=10000
k=4
d=2
p=d+98
a=2
mu_star_EG2=cbind(matrix(c(1,1,-1,-1,1,-1,1,-1),4,2)*a,matrix(0,k,p-d))

##different methods for example 2
set.seed(11)
data=EG2_data(n=10000,p=p,k=k, mustar=mu_star_EG2,sigma=c(rep(1,d),rep(4,p-d)),missing=c(rep(0.3,d),rep(0.3,p-d)))
#kmeans
kmres=kmeans(data$Orig,centers=k, nstart = 100, iter.max = 100)
#kpod
kpodres=fun_kpod(x=data$Missing,k=k,nstart = 100, init_method = "random")
#proposed method (g)
delta=0.01
cmat_l2=apply(kpodres$centers, 2, function(z) sqrt( sum( z^2  ) )  )
weight_vec=1/( (cmat_l2>=delta)*cmat_l2 + (cmat_l2<delta)*delta  )
rkpodres_g_rand=fun_rkpod_g(x=data$Missing,k=k,lambda = 0.1, weight_vec=weight_vec, nstart=100, init_method = "random")

##draw figure 2
par(mar=c(3,3.5,1,1))
#kmeans
plot(x=c(1:p),y=apply(kmres$centers, 2, function(z) sqrt( sum( z^2  ) )  )[1:p],
     pch=3,cex=0.8,lwd=1, col="black",ylim=c(0,4),xlab="",ylab="",xaxt="n",cex.axis=0.8)
title(ylab=TeX("Estimated value of $\\|M_{(j)}\\|$"),line=1.8, cex.lab=1)
title(xlab=TeX("Index $j$") , line=2)
axis(1, at = c(1,20,40,60,80,100),cex.axis=0.8)
lines(x=c(1:p),y=rep(0,p),lty=2,col=grey(0.5))
colvec=kmres$cluster
plot(data$Orig, pch = 16, col = colvec, asp = 1, cex = 0.35, xlim = c(-7,7), ylim = c(-7,7),
     xlab="",ylab="")
#kpod
plot(x=c(1:p),y=apply(kpodres$centers, 2, function(z) sqrt( sum( z^2  ) )  )[1:p],
     pch=2,cex=0.5,lwd=1,  col="red",ylim=c(0,4),xlab="",ylab="",xaxt="n",cex.axis=0.8)
title(ylab=TeX("Estimated value of $\\|M_{(j)}\\|$"),line=1.8, cex.lab=1)
title(xlab=TeX("Index $j$") , line=2)
axis(1, at = c(1,20,40,60,80,100),cex.axis=0.8)
lines(x=c(1:p),y=rep(0,p),lty=2,col=grey(0.5))
colvec=rep(1,n)
for (i in 1:n){
  if (kpodres$cluster[i]==1){colvec[i]=4}
  if (kpodres$cluster[i]==2){colvec[i]=1}
  if (kpodres$cluster[i]==3){colvec[i]=3}
  if (kpodres$cluster[i]==4){colvec[i]=2}
}
plot(data$Orig, pch = 16, col = colvec, asp = 1, cex = 0.35, xlim = c(-7,7), ylim = c(-7,7),
     xlab="",ylab="")
#proposed method (g)
plot(x=c(1:p),y=apply(rkpodres_g_rand$centers, 2, function(z) sqrt( sum( z^2  ) )  )[1:p],
     pch=1,cex=0.8,lwd=1, col="green4",ylim=c(0,4),xlab="",ylab="",xaxt="n",cex.axis=0.8)
title(ylab=TeX("Estimated value of $\\|M_{(j)}\\|$"),line=1.8, cex.lab=1)
title(xlab=TeX("Index $j$") , line=2)
axis(1, at = c(1,20,40,60,80,100),cex.axis=0.8)
lines(x=c(1:p),y=rep(0,p),lty=2,col=grey(0.5))
colvec=rep(1,n)
for (i in 1:n){
  if (rkpodres_g_rand$cluster[i]==1){colvec[i]=1}
  if (rkpodres_g_rand$cluster[i]==2){colvec[i]=4}
  if (rkpodres_g_rand$cluster[i]==3){colvec[i]=3}
  if (rkpodres_g_rand$cluster[i]==4){colvec[i]=2}
}
plot(data$Orig, pch = 16, col = colvec, asp = 1, cex = 0.35, xlim = c(-7,7), ylim = c(-7,7),
     xlab="",ylab="")


