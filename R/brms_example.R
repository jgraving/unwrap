require(circular)#load packages
require(brms)
source("unwrap_family.R")

dt_vonmises = read.table('vonmises_data.csv',#load data
                         sep = ',',
                         header = TRUE)
                       
form_vonmises = bf(y~x, kappa~x, family = unwrap_von_mises)#formula for a circular regression
prior_vonmises =  prior('unwrap_von_mises(0,1e-16)', class = 'Intercept', dpar = 'mu') +
                     prior('normal(0.0,pi()/2)', class = 'b', dpar = 'mu') +
                     prior('normal(log1p_exp(2.0),1.0)', class = 'Intercept', dpar = 'kappa') +
                     prior('normal(0.0,0.5)', class = 'b', dpar = 'kappa')#unbiased priors
                     
model_vonmises = brm(formula = form_vonmises,
                     data = dt_vonmises,
                     family = unwrap_von_mises,
                     stanvars = stan_unwrap_fun + modulo_circular_fun,
                     prior = prior_vonmises)#fit the model