

####---X is a n*p matrix with each row is the pre-experiment covariates of a subject
####---Y is a n*2 matrix. The first column of y represent the successes.


beta_binom_cov = function(n1, n2, X_null, X_alt, nu = 10 ^ (seq(-1, 4, 0.25)), sigma = 4){
  Y = rbind(n1, n2)
  data_0 = Y[, 1]
  data_1 = Y[, 2]
  res1 = 0
  res2 = 0
  res1 = compute_m_cov(data_0, data_1, X_null, nu, sigma)
  res2 = compute_m_cov(data_0, data_1, X_alt, nu, sigma)
  return(exp(res2 - res1) )
}

 

 
 
