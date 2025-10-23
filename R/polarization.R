# A model for the effect of Degree of Linear Polarization on Dung Beetle Orientation Accuracy -----
#   
### Details
#   
# AUTHOR: James Foster 2025 10 23
# 
# MODIFIED: James Foster 2025 10 23
# 
# DESCRIPTION: Fit a hierarchical maximum-likelihood von Mises to dung beetle exit angles under a rotatable polarized light stimulus.
# Modified from [beetles.ipynb](https://github.com/jgraving/unwrap/notebooks/)
# 
# 
# INPUTS:   `'DLdata201611.csv'`
# 
# OUTPUTS:  Plots and test statistics
# 
# CHANGES: 
#   - 
# 
# REFERENCES:
#   - Foster, J.J., $et~al$.  (2019). 
# Orienting to polarized light at night – matching lunar skylight to performance in a nocturnal beetle.
# Journal of Experimental Biology 222, jeb188532. 
# DOI:[10.1242/jeb.188532](https://doi.org/10.1242/jeb.188532)
# 
# - Sayin S, ... Graving JM, $et~al$. (2025)
# The behavioral mechanisms governing collective motion in swarming locusts.
# Science 387,995-1000
# DOI:[10.1126/science.adq7832](https://doi.org/10.1126/science.adq7832)
# 
# - Graving JM & Foster JJ (in preparation)
# Unwrapping Circular Statistics: Bayesian Linear Models for Circular Data
# 
# ---
#   # TODO
#   -

# Set up workspace --------------------------------------------------------
angle_unit = 'degrees'
angle_rot = 'clock'

#Check the operating system and assign a logical flag (T or F)
sys_win = Sys.info()[['sysname']] == 'Windows'
#find file location
here_path = tryCatch(expr = #look in the folder containing this file: sys.frame(1)$ofile
                       {file.path(dirname(sys.frame(1)$ofile))},
                     error = function(e)
                     {#if that fails, try to find the "Documents" folder
                       file.path(Sys.getenv('HOME'))
                     }
)

tryCatch(expr = #try to load functions from the folder containing this file
           {
             source(file = file.path(here_path,
                                     'unwrap_functions.R',
                                     fsep = if(sys_win){'\\'}else{'/'}
             )
             )
           },
         error = function(e)
         {#if that fails, ask the user
           
           path_functions =  
             if(sys_win){#choose.files is only available on Windows
               message('\n\nPlease select the "unwrap_functions.R" file\n\n')
               Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
               choose.files(
                 default = file.path(here_path, "*.R"),#For some reason this is not possible in the "root" user
                 caption = 'Please select the "unwrap_functions.R" file'
               )
             }else{
               message('\n\nPlease select the "unwrap_functions.R" file\n\n')
               Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
               file.choose(new=F)
             }
           #show the user the path they have selected
           if(is.null(path_functions) | !length(path_functions))
           {stop('No file selected.')}else
           {print(path_functions)}
           source(path_functions)
         }
)

path_file = file.path(here_path, "Notebooks/DLdata201611.csv")
if(file.exists(path_file))
{
  print(path_file)
}else
{
  # set path to file
  if(sys_win){#choose.files is only available on Windows
    message('\n\nPlease select the ".csv" file\n\n')
    Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
    path_file = choose.files(
      default = file.path(here_path, "*.csv"),#For some reason this is not possible in the "root" user
      caption = 'Please select the "DLdata201611.csv" file'
    )
  }else{
    message('\n\nPlease select the "DLdata201611.csv" file\n\n')
    Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
    path_file = file.choose(new=F)
  }
  #show the user the path they have selected
  if(is.null(path_file) | !length(path_file))
  {stop('No file selected.')}else
  {print(path_file)}
}
# Custom family --------------------------------------------------------
#In the "unwrap" family, all variables 
unwrap_von_mises = custom_family(
  "unwrap_von_mises", dpars = c("mu", "kappa"),
  links = c('identity',#brms cannot accept custom link functions, do via nl instead
            "softplus"), 
  lb = c(-pi, 0),
  ub = c(pi, NA),
  type = "real"
)

### Stan variables ---------------------------------------------------------


### Functions ------------------------------------------------------------

stan_unwrap_fun = stanvar(scode = "
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
",
                          block = 'functions') + 
  stanvar(scode = "
  real mod_circular(real y) {
    return fmod(y + pi(), 2*pi()) - pi();
  }

  real inv_mod_circular(real y) {
    return mod_circular(y);
  }
",
          block = 'functions') 


### Parameters -----------------------------------------------------------


#concentration of random effects on fixed effect on mean angle
kappamu_param = stanvar(scode = "
real kappamu1;
real kappamu2;
                           ",
                        block = "parameters") + 
  stanvar(scode = "
real kappa_mu1 = log1p_exp(kappamu1);
real kappa_mu2 = log1p_exp(kappamu1+kappamu2);
          ", 
          block = 'genquant')
#concentration of random effects on fixed effect on mean angle
uni_kappamu_param = stanvar(scode = "
real kappamu1;
                           ",
                            block = "parameters") + 
  stanvar(scode = "
real kappa_mu1 = log1p_exp(kappamu1);
          ", 
          block = 'genquant')


stanvars_nopred = stan_unwrap_fun + kappamu_param
stanvars_nopred_uni = stan_unwrap_fun + uni_kappamu_param


# Input Variables ----------------------------------------------------------
all_plots = FALSE # to speed up

# Read in the data and format ---------------------------------------------

#select the reorganised data
cd = read.table(file = path_file, 
                header = T, 
                sep  = ',')
if(all_plots)
{
  View(cd)
}
#night,trial1,trial2,cond1,cond2,deg.of.pol,Beetle
cd = within(cd,
            {
              ID = as.factor(Beetle) # beedance identifier as a factor
              date = as.factor(night) # date as a factor
            }
)

cd_long = with(cd, 
               {
                 data.frame(
                   bearing = c(trial1, trial2),
                   condition = c(cond1, cond2),
                   ID = c(ID, ID),
                   trial = rep(1:2,times = length(trial1)),
                   date = c(date, date),
                   DoLP = c(deg.of.pol, deg.of.pol)
                 )
               }
)

## Specify condition types by rotation angle
cd_long = within(cd_long,
                 {
                  lDoLP = log10(DoLP)#intention: 0.99 -> -0.01; 0.02 -> -0.98
                  signed_angle = Mod360.180(bearing)  # bearing between -180 and 180
                  angle = as.numeric(rad(signed_angle)) # circular format suppresses later warnings
                  pol_bearing = sapply(X = condition,
                                     switch,
                                     PolNorth = 0,
                                     UnpolNorth = 0,
                                     PolEast = 90,
                                     UnpolEast = 90,
                                     NA
                                     )
                  pol_angle = rad(pol_bearing)
                 }
)
u_id = with(cd_long, unique(ID)) # unique beetles
length(u_id)#340 individuals


# Fit model ---------------------------------------------------------------

## Simpler variable names ------------------------------------------------
mod_dt = within(cd_long,
            {
              LD = lDoLP
              PA = pol_angle
              DT = date
              TR = trial
              CD = ifelse(test = pol_angle == 0,
                          yes = 0,
                          no = 1
                          ) # set to 1 for East, 0 for North
            }
)



## Formula ---------------------------------------------------------------

#unimodal model

formula_uni = bf(#modulus may not be necessary, included in lpd function
  formula = mod_circular(angle*2) ~ mu,
  nlf(mu ~ fmu1 + zmu1), 
  fmu1 ~ 1 + CD, # all conditions
  zmu1 ~ 0 + ID, # mean angle combines fixed and random effects
  nlf(kappa ~ k1),
  k1 ~ 1 + LD + (1 |ID), #for kappa this occurs in linear space, and BRMS can set it up automatically
  family = unwrap_von_mises,#mod mu, kappa via the softplus
  nl = TRUE)#to accept user-defined extra parameters (zmu) we need to treat the formula as nonlinear

formula_mix = bf(#modulus may not be necessary, included in lpd function
  formula = angle ~ mu,
  nlf(mu1 ~ fmu1 + zmu1), 
  nlf(mu2 ~ fmu1 + fmu2 + zmu1), #assume it is the same, but with some small increment
  fmu1 ~ 1 + CD, # all conditions
  zmu1 ~ 0 + ID, # mean angle combines fixed and random effects
  # fmu2 ~ 1, #just UV-dim
  # zmu2 ~ 0 + ID, # mean angle combines fixed and random effects
  nlf(kappa1 ~ k1),
  nlf(kappa2 ~ k1), #for 2nd mean, this is specific to individual and condition combination
  k1 ~ 1 + LD + (1 |ID), #for kappa this occurs in linear space, and BRMS can set it up automatically
  theta1 ~ 1 + (1 |ID), #for mixture weighting, this is specific to individual and condition combination
  family = mixture(unwrap_von_mises,#mod mu, kappa via the softplus
                   unwrap_von_mises
  ),#mixture, specified by the log ratio of theta1 : theta2
  nl = TRUE)#to accept user-defined extra parameters (zmu) we need to treat the formula as nonlinear


## Priors ----------------------------------------------------------------

#STRATEGY
#bias mu1 to near 0° to help find any population mean
#force theta to near 0
#could bias effect of condition to 90° (strong expectation)

#priors for mu
pr_mu_uni = 
  prior(normal(0,pi()/2), class = b,  nlpar = 'fmu1', coef = 'Intercept') + # closer to 0°
  prior(normal(pi()/2,pi()/2), class = b,  nlpar = 'fmu1', coef = 'CD') + # closer to 90°
  set_prior(paste("target +=", 
                  'unwrap_von_mises_vect_lpdf(b_zmu1 | 0, log1p_exp(kappamu1))',
                  '+ normal_lpdf(b_zmu1 | 0, 2*pi())'# additional prior to keep estimates from walking around the circle
  ),
  check = FALSE) +
  set_prior("target += normal_lpdf(kappamu1 | 3.0, 3.0)", #prior to higher values, indiv differences should be small
            check = FALSE) 
pr_mu_mix = 
  pr_mu_uni + #same as unimodal, plus priors for 2nd mean
  prior(normal(-pi(),pi()/3), class = b, nlpar = 'fmu2', coef = 'Intercept') + # closer to 180°
  set_prior(paste("target +=", 
                  'unwrap_von_mises_vect_lpdf(b_zmu1 | 0, log1p_exp(kappamu1))',
                  '+ normal_lpdf(b_zmu1 | 0, 2*pi())'# additional prior to keep estimates from walking around the circle
  ),
  check = FALSE) #+
  # set_prior("target += normal_lpdf(kappamu1 | 3.0, 3.0)", #prior to lower values, indiv differences could be large
  #           check = FALSE) + 
  # set_prior(paste("target +=", 
  #                 'unwrap_von_mises_vect_lpdf(b_zmu2 | 0, log1p_exp(kappamu1 + kappamu2))',
  #                 '+ normal_lpdf(b_zmu2 | 0, 2*pi())'# additional prior to keep estimates from walking around the circle
  # ),
  # check = FALSE) +
  # set_prior("target += normal_lpdf(kappamu2 | 0.0, 1.0)", #prior to higher values, indiv differences should be small
  #           check = FALSE)
#priors for kappa
pr_kappa_uni = 
  prior(normal(3.0,3.0), class = b, nlpar = 'k1', coef = 'Intercept') + # bias to oriented
  prior(normal(0.0,3.0), class = b, nlpar = 'k1', coef = 'LD') + # bias to oriented
  prior(student_t(3, 0, 3.0), class = sd, nlpar = 'k1')  # weak bias no differences
pr_kappa_mix = pr_kappa_uni  #identical
#priors for theta (mixture weight)
pr_theta_mix = 
  # prior(normal(1.0,0.25), class = Intercept, dpar = 'theta1') + # force to mu1 as primary #20251007 approx 30% of data are at mu2, so setting to 1.0
  prior(normal(0.0,0.5), class = Intercept, dpar = 'theta1') + # force to mu1 as primary #20251007 approx 30% of data are at mu2, so setting to 1.0
  prior(student_t(3, 0, 0.5), class = sd, dpar = 'theta1') #should this be so small?
# prior(student_t(3, 0, 3.0), class = sd, dpar = 'theta1') #should this be so small?

#all unimodal priors
pr_uni = pr_mu_uni + pr_kappa_uni
#all mixture priors
pr_mix = pr_mu_mix + pr_kappa_mix + pr_theta_mix



## Run model -------------------------------------------------------------
wup = 500
sam = 200


### Unimodal -------------------------------------------------------------
#very long compile time
system.time(
  {
    
    np_uni = brm(formula = formula_uni,
                 data = mod_dt,
                 prior = pr_uni,
                 stanvars =stanvars_nopred_uni,
                 warmup = wup,#may be necessary 
                 iter = wup+sam, #doesn't take a lot of runs
                 chains = 4, # 4 chains in parallel
                 cores = 4, # on 4 CPUs
                 control = list(adapt_delta = 0.90),
                 refresh = 100, # echo chain progress every n iterations
                 silent = 1, # echo some Stan messages
                 backend = 'cmdstanr')
    
  }
)

save(np_uni,
     file = file.path(dirname(path_file),
                      'polarization_uni.RData')
)

sm_uni = summary(np_uni, robust = TRUE)
print(sm_uni$fixed[c('fmu1_Intercept',
                     'fmu1_CD',
                     'k1_Intercept',
                     'k1_LD'),],
      digits = 3)

UnwrapRhats(np_uni, variable = 'fmu')
summary(UnwrapRhats(np_uni), variable = 'zmu')

plot(x = sort(unique(mod_dt$DoLP)),
     y = A1(
         softplus( 
           sort(unique(mod_dt$LD))*
           sm_uni$fixed['k1_LD','Estimate']+
           sm_uni$fixed['k1_Intercept','Estimate']
         )
       ),
     xlab = 'degree of polarization',
     ylab = 'predicted kappa')
abline(h = c(0,1), v = c(0,1))


### Bimodal -------------------------------------------------------------
#warning, takes nearly 15 minutes for 700 iterations!
system.time(
  {
    
    np_mix = brm(formula = formula_mix,
                 data = mod_dt,
                 prior = pr_mix,
                 stanvars =stanvars_nopred,
                 warmup = wup,#may be necessary 
                 iter = wup+sam, #doesn't take a lot of runs
                 chains = 4, # 4 chains in parallel
                 cores = 4, # on 4 CPUs
                 control = list(adapt_delta = 0.90),
                 refresh = 100, # echo chain progress every n iterations
                 silent = 1, # echo some Stan messages
                 backend = 'cmdstanr')
    
  }
)

save(np_mix,
     file = file.path(dirname(path_file),
                      'polarization_mix.RData')
)

sm_mix = summary(np_mix, robust = TRUE)
print(sm_mix$fixed[c('fmu1_Intercept',
                     'fmu1_CD',
                     'k1_Intercept',
                     'k1_LD'),],
      digits = 3)

UnwrapRhats(np_mix, variable = 'fmu')
summary(UnwrapRhats(np_mix), variable = 'zmu')


#  Plot model coefficients ------------------------------------------------
plot(np_mix)


## Model comparison ------------------------------------------------------
#better predictions should justify fitting a mixture model 
loo_uni = loo(np_uni)
loo_mix = loo(np_mix)
lc_unimix = loo_compare(loo_uni, loo_mix)
print(lc_unimix)
#the mixture model has higher predictive power
pnorm(q = lc_unimix[2,1], sd = lc_unimix[2,2],lower.tail = FALSE)


### Plot model comparison ------------------------------------------------
lc_plot = data.frame(elpd = c(loo_uni$estimates['elpd_loo','Estimate'],
                              loo_mix$estimates['elpd_loo','Estimate'],
                              loo_uni$estimates['elpd_loo','Estimate'] - lc_unimix[2,'elpd_diff']),
                     se = c(loo_uni$estimates['elpd_loo','SE'],
                            loo_mix$estimates['elpd_loo','SE'],
                            lc_unimix[2,'se_diff'])
)

par(mar = c(0,4,0,4),
    mfrow = c(1,1))
plot(x = 1:dim(lc_plot)[1],
     y = lc_plot$elpd,
     xlab = '',
     ylab = 'expected log predictive density',
     xlim = c(1,dim(lc_plot)[1]) + c(-1,1)*0.5,
     ylim = with(lc_plot, {range(elpd+se%*% t(c(-2,2)))}), #within 2sigma of all estimates
     pch = 19,
     col = c('darkred', 'darkblue', 'gray35'),
     cex = 2,
     axes = FALSE)
with(lc_plot,
     {
       arrows(x0 = 1:dim(lc_plot)[1],
              x1 = 1:dim(lc_plot)[1],
              y0 = elpd - se,
              y1 = elpd + se,
              code = 3,
              angle = 90,
              length = 0.1,
              lwd = 3,
              col =   c('darkred', 'darkblue', 'gray35')
       )
     }
)
axis(2,
     at = pretty(c(lc_plot$elpd,
                   with(lc_plot, {range(elpd+se%*% t(c(-2,2)))}))
     )
)
axis(4,
     at = with(lc_plot,
               {
                 seq(from = elpd[1], to = elpd[dim(lc_plot)[1]] + se[dim(lc_plot)[1]]*4, by = 100)
               }
     ),
     labels = with(lc_plot,
                   {
                     seq(from = 0, to =  elpd[dim(lc_plot)[1]] + se[dim(lc_plot)[1]]*4 - elpd[1],  by = 100)
                   }
     )
)
abline(h = lc_plot$elpd[1],
       col = 'gray')
mtext(text = 'ELPD difference',
      side = 4,
      line = 3
)
mtext(side = 1,
      line = -1,
      at = 1:dim(lc_plot)[1],
      text = c('unimodal\nmodel',
               'bimodal\nmodel',
               'difference'),
      col = c('darkred', 'darkblue', 'gray35')
)


## Plot coefficients -----------------------------------------------------

#check divergences
# nt_mix = nuts_params(np_mix)
# bayesplot::mcmc_scatter(full_mix,
#                       np = nt_mix,
#                       regex_pars = c('Intercept_theta1', '^b_fmu1_Intercept')
#                       )

if(all_plots)
{
  
  #main effects means
  #weighting (logistic scaled)
  plot(np_mix,
       variable = 'Intercept_theta1',
       transform = plogis)
  plot(np_mix,
       variable = '^b_theta1',
       regex = TRUE) 
  #primary mu
  plot(np_mix,
       variable = '^b_fmu1', # all effects have similar names
       regex = TRUE,
       nvariables = 4,
       transform = unwrap_circular_deg) 
  #secondary mu
  plot(np_mix,
       variable = '^b_fmu2', # all effects have similar names
       regex = TRUE,
       nvariables = 5,
       transform = unwrap_circular_deg) 
  #primary kappa
  plot(np_mix,
       variable = '^b_k1',
       regex = TRUE)#main effects means converge well
  plot(np_mix,
       variable = '^sd_ID__theta1',
       regex = TRUE)
  plot(np_mix,
       variable = '^sd_ID__k1',
       regex = TRUE)
  plot(np_mix,
       variable = '^kappamu1',
       regex = TRUE)
  # plot(np_mix,
  #      variable = '^kappamu2',
  #      regex = TRUE)
  
  #Many random effects of mu
  plot(np_mix,
       variable = '^b_zmu1',
       regex = TRUE,
       ask = FALSE,
       transform = unwrap_circular_deg)
  # plot(np_mix,
  #      variable = '^b_zmu2',
  #      regex = TRUE,
  #      ask = FALSE,
  #      transform = unwrap_circular_deg)
}

# Plot conditional effects ----
# for calculating marginal effects/conditional expectations
posterior_epred_unwrap_von_mises =   function(draws,component="all") {

        mu <- brms::get_dpar(draws, "mu")

        kappa <- brms::get_dpar(draws, "kappa")

        kappa = softplus(kappa)

  }
cond_mix = conditional_effects(np_mix, dpar = 'kappa1')
# conditional_effects(np_mix)

