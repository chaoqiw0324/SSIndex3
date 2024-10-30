#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]] 

double kernal(double dx) 
{
  double out = 0.0;
  if ((dx <= 1.0) && (dx >= -1.0)) {
    out = 315 * (3 - 11 * dx * dx) * (1 - dx * dx)* (1 - dx * dx)* (1 - dx * dx) / 512; // Quartic/biweight
    // out = 0.75 * (1 - dx * dx); // Epanechnikov
    // out = 15 * (1 - dx * dx) * (1 - dx * dx) / 16 
  }
  // if(out<0){
  //   out = 1e-5;
  // }
  return(out) ;
}




// [[Rcpp::export]] 
double M2(const arma::vec& beta,
          const arma::mat& Z,
          const arma::vec& T2,
          const arma::vec& C2,
          const arma::vec& Y2,
          double h,
          double ht)
{
  int n2 = T2.n_elem;
  arma::mat kmat(n2, n2);
  arma::mat bzmat(n2, n2);
  arma::mat imat(n2, n2);
  
  kmat.zeros();
  bzmat.zeros();
  imat.zeros();
  
  arma::vec bz = Z * beta;
  
  
  for(int i = 0; i < n2; i++)
  {
    for(int j = 0; j < n2; j++)
    {
      imat(i,j) = T2(j) <= T2(i) && T2(i) <= C2(j);
      imat(i,j) = imat(i,j)*Y2(j);
      double x = (T2(i)-T2(j))/ht;
      if(x < 1 && x>-1)
      {
        kmat(i,j) = kernal(x)/ht*Y2(j);       
      }
      x = (bz(i)-bz(j))/h;
      if(x < 1 && x>-1)
      {
        bzmat(i,j) = kernal(x)/h;
      }
    }
    
  }
  
  
  
  arma::vec ov(n2);
  ov.ones();
  kmat.diag() = Y2%ov*0.75/ht;
  bzmat.diag() = ov*0.75/h;
  imat.diag() = Y2%ov;
  
  
  // double m1 = arma::sum(log(arma::sum(kmat % bzmat))) - arma::sum( log(arma::sum(kmat % imat)) );
  
  arma::vec num = arma::sum(kmat % bzmat,1);
  arma::vec den = arma::sum(bzmat % imat,1);
  arma::vec r = num/den;
  r.elem(arma::find(r < 0)).fill(1e-3);
  
  double m1 = arma::sum(Y2 % log(r));
  
  
  return m1 ;
}






// [[Rcpp::export]] 

double shapeFun3(int n, 
                 arma::vec& m, 
                 arma::vec& midx, 
                 arma::vec& tij, 
                 arma::vec& yi, 
                 arma::mat& xb,
                 double x, 
                 double t, 
                 double h, 
                 arma::vec& w ,
                 arma::vec& med,
                 double result
) 
{
  
  double nu = 0.0;
  double de = 0.0;
  for (int i = 0; i < n; i++) {
    for (int k = 0; k < m[i]; k++) {
      if (tij[midx[i] + k] >= t) {
        de = 0.0;
        nu = 0.0;
        nu = w[i] * kernal((x - xb[i]) / h) / h * med[midx[i] + k];
        for (int j = 0; j < n; j++) {
          for (int l = 0; l < m[j]; l++) {
            /* if (tij[midx[i] + k] == tij[midx[j] + l])   */
            /*   nu += w[j] * kernal((x[0] - xb[j]) / h[0]) / h[0];    */
            if (tij[midx[i] + k] >= tij[midx[j] + l] && tij[midx[i] + k] <= yi[j])
              de += w[j] * kernal((x - xb[j]) / h) / h * med[midx[j] + l];
          }
        }
        if(de == 0) {
          result += 0;
        } else {
          result += nu / de; 
        }
      }
    }
  }
  return result;
}



// [[Rcpp::export]] 

double size(int n, 
            arma::mat& xr,
            arma::vec& mFhat,
            arma::vec& w
) 
{
  double result = 0;
  for (int i = 0; i < n; i++) {
    // tmp = shapeFun(n, m, midx, tij, yi, xb, xb[i], yi[i]);    
    for (int j = 0; j < n; j++) {
      if (xr[i] > xr[j]) {
        result += w[i] * w[j] * mFhat[i];
      }
      if (xr[i] == xr[j] && i != j) result += 0.25 * w[i] * w[j] * (mFhat[i] + mFhat[j]);
    }
  }
  return result;
}
