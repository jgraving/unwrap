require(brms)

#set up required Stan functions
#circular modulo in Stan code
modulo_circular_fun = stanvar(scode = "
  real modulo_circular(real y) {
    return fmod(y + pi(), 2*pi()) - pi();
  }
",
                           block = 'functions')

#custom likelihood function using the shifted modulo link
stan_unwrap_fun = stanvar(scode = "
  real unwrap_von_mises_lpdf(real y, real mu, real kappa) {
    return von_mises_lpdf(y | modulo_circular(mu), kappa);
  }
  real unwrap_von_mises_rng(real mu, real kappa) {
    return von_mises_rng( modulo_circular(mu) , kappa);
  }
",
                          block = 'functions') 
#define the custom family
unwrap_von_mises = custom_family(
  name = "unwrap_von_mises",
  dpars = c("mu",
            "kappa"),
  links = c('identity',#N.B. the link function is defined via the LPD function
            "softplus"), 
  lb = c(-pi,#lower bound of mu should be -pi
         0),
  ub = c(pi,#upper bound of mu should be pi
         NA),
  type = "real",#takes continuous response data
)
