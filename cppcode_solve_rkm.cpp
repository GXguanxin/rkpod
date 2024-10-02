// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// [[Rcpp::export]]
List solve_rkm_g_cpp(const arma::mat x, int k, arma::vec weight_vec,
                     double lambda, const arma::mat initial_centers,
                     int itermax, double eps){
  arma::mat xmat = x, cmat_old = initial_centers; 
  int n = xmat.n_rows, p = xmat.n_cols; 
  arma::mat umat = arma::zeros(n,k), cmat = arma::zeros(k,p);
  arma::vec loss_path(itermax), clustersize(k), cmat_l2(p), cluster(n); 
  double dist = 0, tmp_dist = 0, lambdaj = 0, loss = 0;
  int kk_star = 0, final_t = 0;
  
  for(int t = 0; t < itermax; ++t){//main iteration
    //given cmat, update umat
    umat = arma::zeros(n,k);
    for(int i = 0; i < n; ++i){
      dist = 0;
      kk_star = 0;
      for(int kk = 0; kk < k; ++kk){
        tmp_dist = norm(xmat.row(i) - cmat_old.row(kk),2);
        if(kk == 0){
          dist = tmp_dist;
        }else if(tmp_dist < dist){
          dist = tmp_dist;
          kk_star = kk;
        }
      }
      umat(i,kk_star) = 1;
      cluster(i) = kk_star+1;
    }
    for(int kk = 0; kk < k; ++kk){
      clustersize(kk)=sum(umat.col(kk));
    }
    
    
    //given umat, update cmat (for each colomn of cmat, solve ridge reg)
    cmat = arma::zeros(k,p);
    cmat_l2 = arma::zeros(p);
    for(int j = 0; j < p; ++j){
      cmat_l2(j) = norm(cmat_old.col(j),2);
      lambdaj = ( n*lambda*(weight_vec(j)) )/( 2*(cmat_l2(j)) );
      for(int kk = 0; kk < k; ++kk){
        cmat(kk,j) = sum((umat.col(kk))%(xmat.col(j)))/(clustersize(kk)+lambdaj);
      }
    }
    
    //record
    loss = pow( norm(xmat - umat*cmat,"fro") , 2 )  + n*lambda*sum( weight_vec%cmat_l2 );
    loss_path(t) = loss;
    
    //stop criterion
    if(t>0){
      if( abs(loss_path(t-1)-loss_path(t))/loss_path(t-1) < eps ){
        final_t = t;
        break;
      }else{
        cmat_old = cmat;
      }
    }else{
      cmat_old = cmat;
    }
    
  }
  
  List result;
  result["cluster"] = cluster;
  result["centers"] = cmat;
  result["membership"] = umat;
  result["loss"] = loss;
  result["loss_path"] = loss_path.subvec(0,final_t);
  
    
  return result;
}




// [[Rcpp::export]]
List solve_rkm_L0_cpp(const arma::mat x, int k, 
                     double lambda, const arma::mat initial_centers,
                     int itermax, double eps){
  arma::mat xmat = x, cmat_old = initial_centers; 
  int n = xmat.n_rows, p = xmat.n_cols; 
  arma::mat umat = arma::zeros(n,k), cmat = arma::zeros(k,p), cmat_star = arma::zeros(k,p);
  arma::vec leftterm(p), loss_path(itermax), clustersize(k), cluster(n); 
  double dist = 0, tmp_dist = 0, rightterm = 0, loss = 0;
  int kk_star = 0, nbactivej = 0, final_t = 0;
  
  for(int j = 0; j < p; ++j){
    leftterm(j) = pow( norm(xmat.col(j),2), 2);
  }
  
  for(int t = 0; t < itermax; ++t){//main iteration
    //given cmat, update umat
    umat = arma::zeros(n,k);
    for(int i = 0; i < n; ++i){
      dist = 0;
      kk_star = 0;
      for(int kk = 0; kk < k; ++kk){
        tmp_dist = norm(xmat.row(i) - cmat_old.row(kk),2);
        if(kk == 0){
          dist = tmp_dist;
        }else if(tmp_dist < dist){
          dist = tmp_dist;
          kk_star = kk;
        }
      }
      umat(i,kk_star) = 1;
      cluster(i) = kk_star+1;
    }
    for(int kk = 0; kk < k; ++kk){
      clustersize(kk)=sum(umat.col(kk));
    }
    
    //given umat, update cmat (for each colomn of cmat, solve ridge reg)
    cmat = arma::zeros(k,p);
    cmat_star = arma::zeros(k,p);
    nbactivej = 0; 
    for(int j = 0; j < p; ++j){
      for(int kk = 0; kk < k; ++kk){
        cmat_star(kk,j) = sum( (umat.col(kk))%(xmat.col(j))  )/clustersize(kk);
      }
      //leftterm = pow( norm(xmat.col(j),2), 2);
      rightterm = sum( pow( xmat.col(j) - umat*(cmat_star.col(j)) , 2)) + n*lambda; 
      if(leftterm(j) > rightterm){
        cmat.col(j) = cmat_star.col(j);
        nbactivej = nbactivej + 1; 
      }
    }
    
    //record
    loss = pow( norm(xmat - umat*cmat,"fro") , 2 )  + n*lambda*nbactivej;
    loss_path(t) = loss;
    
    //stop criterion
    if(t>0){
      if( abs(loss_path(t-1)-loss_path(t))/loss_path(t-1) < eps ){
        final_t = t;
        break;
      }else{
        cmat_old = cmat;
      }
    }else{
      cmat_old = cmat;
    }
    
  }
  
  List result;
  result["cluster"] = cluster;
  result["centers"] = cmat;
  result["membership"] = umat;
  result["loss"] = loss;
  result["loss_path"] = loss_path.subvec(0,final_t);
  
  
  return result;
}