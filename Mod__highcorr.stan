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
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  array[N] int<lower=1> J_1;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_kappa_1;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
}
parameters {
  vector[K_fmu] b_fmu;  // regression coefficients
  vector[K_zmu] b_zmu;  // regression coefficients
  real Intercept_kappa;  // temporary intercept for centered predictors
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  array[M_1] vector[N_1] z_1;  // standardized group-level effects
  
real kappamu;
                           
}
transformed parameters {
  vector[N_1] r_1_kappa_1;  // actual group-level effects
  real lprior = 0;  // prior contributions to the log posterior
  r_1_kappa_1 = (sd_1[1] * (z_1[1]));
  lprior += normal_lpdf(b_fmu | 0,pi()/2);
  lprior += unwrap_von_mises_vect_lpdf(b_zmu | 0, log1p_exp(kappamu));
  lprior += normal_lpdf(Intercept_kappa | 3,2);
  lprior += student_t_lpdf(sd_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
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
    kappa += Intercept_kappa;
    for (n in 1:N) {
      // add more terms to the linear predictor
      kappa[n] += r_1_kappa_1[J_1[n]] * Z_1_kappa_1[n];
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
  target += std_normal_lpdf(z_1[1]);
  target += normal_lpdf(kappamu | 3, 3);
}
generated quantities {
  // actual population-level intercept
  real b_kappa_Intercept = Intercept_kappa;
  
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
          
}

