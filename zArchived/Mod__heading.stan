// generated with brms 2.22.0
functions {
  /* softplus link function inverse to 'log1p_exp'
   * Args:
   *   x: a positive scalar
   * Returns:
   *   a scalar in (-Inf, Inf)
   */
   real log_expm1(real x) {
     return log(expm1(x));
   }
  /* softplus link function inverse to 'log1p_exp' (vectorized)
   * Args:
   *   x: a positive vector
   * Returns:
   *   a vector in (-Inf, Inf)
   */
   vector log_expm1(vector x) {
     return log(expm1(x));
   }
 /* compute correlated group-level effects
  * Args:
  *   z: matrix of unscaled group-level effects
  *   SD: vector of standard deviation parameters
  *   L: cholesky factor correlation matrix
  * Returns:
  *   matrix of scaled group-level effects
  */
  matrix scale_r_cor(matrix z, vector SD, matrix L) {
    // r is stored in another dimension order than z
    return transpose(diag_pre_multiply(SD, L) * z);
  }
  
    real unwrap_von_mises_lpdf(real y, real mu, real kappa) {
      return von_mises_lpdf(y | mod_circular(mu), kappa);
    }
    real unwrap_von_mises_rng(real mu, real kappa) {
      return von_mises_rng( mod_circular(mu) , kappa);
    }
    real unwrap_von_mises_vect_lpdf(vector y, real mu, real kappa) {
    real tmp = 0;
    for(i in 1:size(y))
    {
    tmp = tmp + unwrap_von_mises_lpdf(y[i] | mu, kappa);
    }
      return tmp;
    }
  
  
    real mod_circular(real y) {
      return fmod(y + pi(), 2*pi()) - pi();
    }
  
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K_fmu;  // number of population-level effects
  matrix[N, K_fmu] X_fmu;  // population-level design matrix
  int<lower=1> K_zmu;  // number of population-level effects
  matrix[N, K_zmu] X_zmu;  // population-level design matrix
  int<lower=1> K_kappa;  // number of population-level effects
  matrix[N, K_kappa] X_kappa;  // population-level design matrix
  int<lower=1> Kc_kappa;  // number of population-level effects after centering
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  array[N] int<lower=1> J_1;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_kappa_1;
  vector[N] Z_1_kappa_2;
  int<lower=1> NC_1;  // number of group-level correlations
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  matrix[N, Kc_kappa] Xc_kappa;  // centered version of X_kappa without an intercept
  vector[Kc_kappa] means_X_kappa;  // column means of X_kappa before centering
  for (i in 2:K_kappa) {
    means_X_kappa[i - 1] = mean(X_kappa[, i]);
    Xc_kappa[, i - 1] = X_kappa[, i] - means_X_kappa[i - 1];
  }
}
parameters {
  vector[K_fmu] b_fmu;  // regression coefficients
  vector[K_zmu] b_zmu;  // regression coefficients
  vector[Kc_kappa] b_kappa;  // regression coefficients
  real Intercept_kappa;  // temporary intercept for centered predictors
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  matrix[M_1, N_1] z_1;  // standardized group-level effects
  cholesky_factor_corr[M_1] L_1;  // cholesky factor of correlation matrix
  
real kappamu;
                           
  
real deltakappamu;
                           
}
transformed parameters {
  matrix[N_1, M_1] r_1;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_1] r_1_kappa_1;
  vector[N_1] r_1_kappa_2;
  real lprior = 0;  // prior contributions to the log posterior
  // compute actual group-level effects
  r_1 = scale_r_cor(z_1, sd_1, L_1);
  r_1_kappa_1 = r_1[, 1];
  r_1_kappa_2 = r_1[, 2];
  lprior += normal_lpdf(b_fmu[1] | 0, pi()/1);
  lprior += normal_lpdf(b_fmu[2] | 0, pi()/2);
  lprior += normal_lpdf(b_kappa[1] |  0.0, 3.0);
  lprior += normal_lpdf(Intercept_kappa |  3.0, 3.0);
  lprior += student_t_lpdf(sd_1[1] | 3, 0, 3.0)
    - 1 * student_t_lccdf(0 | 3, 0, 3.0);
  lprior += student_t_lpdf(sd_1[2] | 3, 0, 3.0)
    - 1 * student_t_lccdf(0 | 3, 0, 3.0);
  lprior += lkj_corr_cholesky_lpdf(L_1 | 1);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] nlp_fmu = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] nlp_zmu = rep_vector(0.0, N);
    // initialize non-linear predictor term
    vector[N] mu;
    // initialize linear predictor term
    vector[N] kappa = rep_vector(0.0, N);
    nlp_fmu += X_fmu * b_fmu;
    nlp_zmu += X_zmu * b_zmu;
    kappa += Intercept_kappa + Xc_kappa * b_kappa;
    for (n in 1:N) {
      // add more terms to the linear predictor
      kappa[n] += r_1_kappa_1[J_1[n]] * Z_1_kappa_1[n] + r_1_kappa_2[J_1[n]] * Z_1_kappa_2[n];
    }
    for (n in 1:N) {
      // compute non-linear predictor values
      mu[n] = (nlp_fmu[n] + nlp_zmu[n]);
    }
    kappa = log1p_exp(kappa);
    for (n in 1:N) {
      target += unwrap_von_mises_lpdf(Y[n] | mu[n], kappa[n]);
    }
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(to_vector(z_1));
  target += unwrap_von_mises_vect_lpdf(b_zmu[1:5] | 0, log1p_exp(kappamu)) + normal_lpdf(b_zmu[1:5] | 0, 2*pi());
  target += unwrap_von_mises_vect_lpdf(b_zmu[6:10] | 0, log1p_exp(kappamu+deltakappamu)) + normal_lpdf(b_zmu[6:10] | 0, 2*pi());
  target += normal_lpdf(kappamu | 3.0, 3.0);
  target += normal_lpdf(deltakappamu | 0.0, 3.0);
}
generated quantities {
  // actual population-level intercept
  real b_kappa_Intercept = Intercept_kappa - dot_product(means_X_kappa, b_kappa);
  // compute group-level correlations
  corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1);
  vector<lower=-1,upper=1>[NC_1] cor_1;
  
real kappa_mu = log1p_exp(kappamu);
          
  
  vector [K_fmu] mod_mu;
  for(i in 1:size(mod_mu))
  {
    mod_mu[i] = mod_circular(b_fmu[i]);
  }
  vector [K_zmu] mod_zmu;
  for(i in 1:size(mod_zmu))
  {
    mod_zmu[i] = mod_circular(b_zmu[i]);
  }
          
  
real kappa_mu_delta = log1p_exp(kappamu + deltakappamu);
          
  // extract upper diagonal of correlation matrix
  for (k in 1:M_1) {
    for (j in 1:(k - 1)) {
      cor_1[choose(k - 1, 2) + j] = Cor_1[j, k];
    }
  }
}

