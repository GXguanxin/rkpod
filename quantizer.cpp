// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;


// [[Rcpp::export]]
List qpp_arma(const arma::mat xmat, const int K, const double r){
  double tmp = 0, tmp_dist = 0;
  int n = xmat.n_rows, d = xmat.n_cols;
  arma::mat tmp_cmat = arma::zeros(K,d);
  arma::vec tmp_dvec = arma::ones(n), index_vec(K), all_index = arma::linspace(0, n-1, n), tmp_index(1); 
  
  for(int k = 0; k < K; ++k){//main iteration
    tmp_index = RcppArmadillo::sample(all_index, 1, false, tmp_dvec);
    index_vec(k) = tmp_index(0);
    tmp_cmat(k,arma::span::all) = xmat(tmp_index(0),arma::span::all);
    for(int i = 0; i < n; ++i){
      tmp = sum( pow(xmat(i,arma::span::all) - tmp_cmat(k,arma::span::all), 2) );
      tmp_dist = tmp;
      if(r != 2){
        tmp_dist = pow(tmp, r/2);
      }
      if(k == 0){
        tmp_dvec(i) = tmp_dist;
      }else if(tmp_dvec(i) > tmp_dist){
        tmp_dvec(i) = tmp_dist;
      }
    }//END FOR i
  }//main iteration

  List res;
  res["centers"] = tmp_cmat;
  res["loss"] = mean(tmp_dvec);
  res["index"] = index_vec;

  return res;
}

// [[Rcpp::export]]
List wqpp_arma(const arma::mat xmat, const int K, const double r, const arma::vec wvec){
  double tmp = 0, tmp_dist = 0;
  int n = xmat.n_rows, d = xmat.n_cols;
  arma::mat tmp_cmat = arma::zeros(K,d);
  arma::vec pr_vec = wvec, index_vec(K), all_index = arma::linspace(0, n-1, n), tmp_index(1); 
  
  for(int k = 0; k < K; ++k){//main iteration
    tmp_index = RcppArmadillo::sample(all_index, 1, false, pr_vec);
    index_vec(k) = tmp_index(0);
    tmp_cmat(k,arma::span::all) = xmat(tmp_index(0),arma::span::all);
    for(int i = 0; i < n; ++i){
      tmp = sum( pow(xmat(i,arma::span::all) - tmp_cmat(k,arma::span::all), 2) );
      tmp_dist = tmp;
      if(r != 2){
        tmp_dist = pow(tmp, r/2);
      }
      if(k == 0){
        pr_vec(i) = wvec(i)*tmp_dist;
      }else if(pr_vec(i) > tmp_dist){
        pr_vec(i) = wvec(i)*tmp_dist;
      }
    }//END FOR i
  }//main iteration

  List res;
  res["centers"] = tmp_cmat;
  res["loss"] = mean(pr_vec);
  res["index"] = index_vec;

  return res;
}

// [[Rcpp::export]]
List quantizer_arma(const arma::mat xmat, const arma::mat cmat, const double r, 
                    double eta, int niter, double eps, double tol, int e_type, int prnt){

  double tmp_dist = 0, tmp = 0, tmp_r = r, dist = 0, loss = 0;
  int n = xmat.n_rows, K = cmat.n_rows, tmp_k = 0, tmp_j = niter-1;
  arma::mat tmp_umat = arma::zeros(n,K), tmp_cmat = cmat;
  arma::vec tmp_dsq(n), tmp_w(n), tmp_vec(n), tmp_loss(n);
	arma::vec loss_vec(niter);
	
  for(int j = 0; j < niter; ++j){//main iteration
    tmp_umat = arma::zeros(n,K);
    //Assignment
    for(int i = 0; i < n; ++i){
      dist = 0;
      tmp_k = 0;
      tmp_dsq(i) = 0;
      for(int k = 0; k < K; ++k){
        tmp_dist = norm(xmat(i,arma::span::all) - tmp_cmat(k,arma::span::all));
        if(k == 0){
          dist = tmp_dist;
          tmp_dsq(i) = dist;
        }else if(tmp_dist < dist){
          dist = tmp_dist;
          tmp_dsq(i) = dist;
          tmp_k = k;
        }//END IF
      }//END FOR k
      tmp_umat(i,tmp_k) = 1;
    }
    //Update centroids
    //------------------------------
    /*Compute loss and weights*/
    for(int i = 0; i < n; ++i){
      tmp_loss(i) = pow(tmp_dsq(i),tmp_r);
      if(e_type == 0){
        tmp_w(i) = pow(tmp_dsq(i)+eps,tmp_r-2);
      }else if(e_type == 1){
        if(tmp_dsq(i) < eps){
          tmp_w(i) = pow(tmp_dsq(i)+eps,tmp_r-2);
        }else{
          tmp_w(i) = pow(tmp_dsq(i),tmp_r-2);
        }
      }
    }
    /*Compute centroids*/
    for(int k = 0; k < K; ++k){
      tmp_vec = tmp_w%tmp_umat(arma::span::all,k);
      tmp = sum(tmp_vec);
      tmp_cmat(k,arma::span::all) = (tmp_vec.t()/tmp)*xmat;
    }

    //------------------------------
    loss_vec(j) = mean(tmp_loss);
    loss = loss_vec(j);
    if(prnt == 1){
      Rprintf("%ith-loss = %f \n",j,loss_vec(j));
    }
    if(j > 0){
	     if(abs(loss_vec(j) - loss_vec(j-1))/loss_vec(j-1) < tol){
          tmp_j = j;
          break;
	     }
    }
    tmp_r = eta*tmp_r;
	}//main iteration

  List res;
  res["loss"] = loss;
  res["loss_vec"] = loss_vec.subvec(0,tmp_j);
  res["centers"] = tmp_cmat;
  res["membership"] = tmp_umat;
  res["r"] = tmp_r;

  return res;
}

// [[Rcpp::export]]
List power_quantizer_arma(
      const arma::mat xmat, const arma::mat cmat, const double r, const double a, double eta_r, double eta_a, int niter, 
      int nrep, double eps, double tol, int prnt
    ){

  double tmp_dist = 0, tmp = 0, tmp_r = r, tmp_a = a, loss = 0;
  double c = a*r-2, d = a*r, f = 1 - 1/a, diff = 0;
  int n = xmat.n_rows, K = cmat.n_rows, tmp_j = niter-1;
  arma::mat tmp_wmat = arma::zeros(n,K), tmp_cmat = cmat, old_cmat = cmat;
  arma::vec tmp_vec(n), tmp_loss(n), loss_vec(niter);

  for(int j = 0; j < niter; ++j){//main iteration
    tmp_a = a;
    for(int s = 0; s < nrep; ++s){
      c = tmp_a*tmp_r-2, d = tmp_a*tmp_r, f = 1 - 1/tmp_a;
      //Compute weights
      for(int i = 0; i < n; ++i){
        tmp = 0;
        for(int k = 0; k < K; ++k){
          tmp_dist = norm(xmat(i,arma::span::all) - tmp_cmat(k,arma::span::all));
          tmp_wmat(i,k) = pow(tmp_dist+eps, c);
          tmp += pow(tmp_dist, d)/K;
        }//END FOR k
        tmp_wmat(i,arma::span::all) = tmp_wmat(i,arma::span::all)/pow(tmp,f);
        tmp_loss(i) = pow(tmp, 1/tmp_a);
      }//END FOR i
      //Update centroids
      //------------------------------
      /*Compute centroids*/
      for(int k = 0; k < K; ++k){
        tmp_vec = tmp_wmat(arma::span::all,k);
        tmp = sum(tmp_vec);
        tmp_cmat(k,arma::span::all) = (tmp_vec.t()/tmp)*xmat;
      }
      //------------------------------
      loss_vec(j) = mean(tmp_loss);
      loss = loss_vec(j);
      diff = norm(old_cmat - tmp_cmat)/norm(old_cmat);
      if(prnt == 1){
        Rprintf("%ith-diff = %f \n",j, diff );
      }
      old_cmat = tmp_cmat;
      if(s > 0){
         if(diff < tol){
            tmp_j = j;
            break;
         }
      }
      tmp_a = eta_a*tmp_a;
    }//END FOR s
    tmp_r = eta_r*tmp_r;
    if(eta_r == 1){
      break;
    }
  }//main iteration

  List res;
  res["loss"] = loss;
  res["loss_vec"] = loss_vec;//.subvec(0,tmp_j);
  res["centers"] = tmp_cmat;
  res["r"] = tmp_r;

  return res;
}


    