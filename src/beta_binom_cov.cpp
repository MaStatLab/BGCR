
#include <math.h>
#include <RcppArmadillo.h>
#include "beta_binom_cov.h"

// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

double log_exp_x_plus_exp_y(double x, double y)
{
  double result;
  if( ( std::isinf( fabs(x) ) == 1 ) && ( std::isinf( fabs(y) ) == 0 )  )
    result = y;
  else if ( ( std::isinf( fabs(x) ) == 0 ) && ( std::isinf( fabs(y) ) == 1 )  )
    result = x;
  else if ( ( std::isinf( fabs(x) ) == 1 ) && ( std::isinf( fabs(y) ) == 1 )  )
    result = x;
  else if ( x - y >= 100 ) result = x;
  else if ( x - y <= -100 ) result = y;
  else {
    if (x > y) {
      result = y + log( 1 + exp(x - y) );
    }
    else result = x + log( 1 + exp(y - x) );
  }
  return result;
}

double afun(double nu, double yl, double yr, double theta)
{
  double output;
  output = R::digamma(theta * nu + yl) - R::digamma((1 - theta) * nu + yr) - R::digamma(theta * nu)
           + R::digamma((1-theta) * nu);
  return output;
}


double bfun(double nu, double yl, double yr, double theta)
{
  double output;
  output = R::trigamma(theta * nu + yl) + R::trigamma((1 - theta) * nu + yr) - R::trigamma(theta * nu)
           - R::trigamma((1 - theta) * nu);
  return output;
}


double eval_h(arma::vec beta, arma::mat X, arma::vec data_0, arma::vec data_1, double nu, double sigma)
{
  int n = data_0.n_elem;
  int p = X.n_cols;
  vec theta(n);
  theta = 1 / (1 + exp(-X * beta));

  double y = -0.5 * dot(beta, beta) / pow(sigma, 2) - 0.5 * p * log(2 * M_PI * pow(sigma, 2));
  for(int j = 0; j < n; j++)
  {
    y += R::lbeta(theta(j) * nu + data_0(j), (1 - theta(j)) * nu + data_1(j));
    y -= R::lbeta(theta(j) * nu, (1 - theta(j)) * nu);
  }
  return y;
}


arma::vec truncate(arma::vec x, double a)
{
  int p = x.n_elem;
  for(int i=0; i<p; i++)
  {
    if(x(i) > a)
      x(i) = a;
    else if(x(i) < -a)
      x(i) = -a;
  }

  return x;
}


arma::mat NewtonRaphson(arma::vec data_0, arma::vec data_1, arma::mat X, double nu, double sigma)
{
  int p = X.n_cols;
  int n = X.n_rows;
  vec beta0 = zeros<vec>(p);
  vec beta1 = zeros<vec>(p);
  mat W1(n,n,fill::zeros);
  mat W2(n,n,fill::zeros);
  mat id = eye<mat>(p,p);
  mat Hess(p,p,fill::zeros);
  mat result(p,p+1,fill::zeros);

  double max_step = 2;
  vec step = zeros<vec>(p);

  vec theta(n);
  vec z(n);
  theta = 1 / (1 + exp(-X * beta0)); //logit link

  double tolerance = 10E-10;
  int maxIter = 20;

  for(int i = 0; i < maxIter; i++)
  {
    for(int j = 0; j < n; j++)
    {
      double a,b,theta_temp;
      a = afun(nu, data_0(j), data_1(j), theta(j));
      b = bfun(nu, data_0(j), data_1(j), theta(j));
      theta_temp = theta(j);

      z(j) = theta_temp * (1 - theta_temp);
      W1(j,j) = a;
      W2(j,j) = - nu * b * pow(z(j), 2) - a * z(j) * (1 - 2 * theta_temp);
    }

    Hess =   nu * (X.t() * W2 * X + id / (nu * pow(sigma, 2))) ;
    step = truncate(inv(Hess) * ( X.t() * W1 * z - beta0 / (nu * pow(sigma, 2)) ) * nu, max_step) ;
    beta1 = beta0 + step ;

    if(norm(beta1-beta0, 2) < tolerance)
      i = maxIter;
    else
    {
      beta0 = beta1;
      theta =  1 / (1 + exp(-X * beta0));
    }
  }

  result.col(0) = beta1;
  result.cols(1,p) = Hess;
  return  result;
}


// [[Rcpp::export]]

double compute_m_cov(arma::vec data_0, arma::vec data_1, arma::mat X, arma::vec nu_vec, double sigma)
{
  int n_grid = nu_vec.n_elem;
  int p = X.n_cols;

  double marginal = log(0.0);

  for(int g = 0; g < n_grid; g++)
  {
    mat nr = NewtonRaphson(data_0, data_1, X, nu_vec(g), sigma);
    vec mle = nr.col(0);
    mat hess = nr.cols(1, p);

    marginal = log_exp_x_plus_exp_y(marginal, eval_h(mle, X, data_0, data_1, nu_vec(g), sigma) +
              0.5 * p * log(2.0 * M_PI) - 0.5 * log(det(hess)) - log(n_grid));
  }

  return marginal;
}







