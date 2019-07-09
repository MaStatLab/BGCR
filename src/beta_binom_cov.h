
#ifndef HELPERS_H
#define HELPERS_H

#include <math.h>
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

double log_exp_x_plus_exp_y(double x, double y);

double afun(double nu, double yl, double yr, double theta);

double bfun(double nu, double yl, double yr, double theta);

double eval_h(arma::vec beta, arma::mat X, arma::vec data_0, arma::vec data_1, double nu, double sigma);

arma::mat NewtonRaphson(arma::vec data_0, arma::vec data_1, arma::mat X, double nu, double sigma );

double compute_m_cov(arma::vec data_0, arma::vec data_1, arma::mat X, vec nu_vec, double sigma);

#endif
