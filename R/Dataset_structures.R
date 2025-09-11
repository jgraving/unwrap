#FOR A 'CLEAN' RUN, PRESS ctrl+shift+F10 to RESTART Rstudio
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2025 09 11
#     MODIFIED:	James Foster              DATE: 2025 09 11
#
#  DESCRIPTION: Attempt to run a two-way interaction model on the GUV dances data
#               using the circular modulo modelling method devel. by Jake Graving.
#               Includes a secondary mean to account for bimodality for UV-dim stimulus.
#               Modified from GUV_CircMod_v2.R
#               
#       INPUTS: 
#               
#      OUTPUTS: Plots and test statistics
#
#	   CHANGES: - 
#
#   REFERENCES: Sayin S, Couzin-Fuchs E, Petelski I, Günzel Y, Salahshour M, 
#               Lee CY, Graving JM, Li L, Deussen O, Sword GA, et al. (2025)
#               The behavioral mechanisms governing collective motion in swarming locusts. 
#               Science. 387(6737):995–791
#               
#               Gabry J, Češnovar R, Johnson A (2022). 
#               cmdstanr: R Interface to 'CmdStan'.
#               https://mc-stan.org/cmdstanr/
# 
#               Bürkner, P.-C. (2018). 
#               Advanced Bayesian Multilevel Modeling with the R Package brms. 
#               The R Journal 10, 395–411.
# 
#               Carpenter, B., Gelman, A., Hoffman, M. D., Lee, D., Goodrich, B., 
#               Betancourt, M., Brubaker, M., Guo, J., Li, P. and Riddell, A. (2017). 
#               Stan: A Probabilistic Programming Language. 
#               Journal of Statistical Software 76 doi: 10.18637/jss.v076.i01
# 
#
#    EXAMPLES:  
#
# 
#TODO   ---------------------------------------------
#TODO   
#- Separate functions
#- Rmarkdown version

set.seed(0120810506)#ISBN Batschelet, 1981


# Set up plot colours -----------------------------------------------------

col_kappa = "#1E78B5"
col_rho = '#F08024'
col_sd = '#E74A29'
col_sd2 = '#E57461'
col_pdf = adjustcolor(col = '#21A885',
                      alpha.f = 0.7)
col_obs = '#3E1F51'

# Set up functions --------------------------------------------------------


require(circular)
require(brms)
#function for estimating sd from kappa
MardiaSD = function(k,
                    method = 'analytical',
                    n = 1e4)
{
  switch(method,
         simulation = 
           circular::sd.circular(
             circular::rvonmises(n = n, 
                                 mu = circular::circular(pi, # this is directly in the middle of (0,2pi)
                                                         template = 'none'),
                                 kappa = k)
           ),
         analytical = sqrt( -2* log( circular::A1(k))),
         sqrt( -2* log( circular::A1(k)))
  )
}

#estimate normal SD from kappa
NormalSD = function(k,
                    n = 1e4)
{
  sd(
    as.numeric(
      circular::rvonmises(n = n, 
                          mu = circular::circular(pi, # this is directly in the middle of (0,2pi)
                                                  template = 'none'),
                          kappa = k)
    )
  )
}

#convert angles to signed angles in (-180, 180)
Mod360.180 = function(x)
{#use atan2 to convert any angle to the range (-180,180)
  deg(
    atan2(y = sin(rad(x)),
          x = cos(rad(x))
    )
  )
}


#generic mean angle simulator
MeanRvm = function(n, #representative sample size
                   mu = circular(0), #mean (defaults to 0rad)
                   kappa, #kappa required
                   au = 'degrees', #units
                   ar = 'clock') #rotation direction
{
  mean.circular(rvonmises(n = n, 
                          mu = circular(mu, units = au, rotation = ar), 
                          kappa = kappa,
                          control.circular = list(units = au, rotation = ar)))
}


#Simulate confidence intervals for a unimodal or bimodal distribution
#fitted to a vector of "angles"
CI_vM = function(angles, #vector of angles fitted (used for sample size)
                 m1, #primary mean
                 k1, #primary concentration
                 m2 = NA, #secondary mean (ignored if NULL or NA)
                 k2 = NA, #secondary kappa
                 w1 = 1, #weighting of primary mean
                 force_mu = FALSE, #force median at true mu?
                 n = 1e4, #number of simulations
                 au = 'degrees', 
                 ar = 'clock',
                 calc_q = TRUE,
                 alternative = 'one.sided', #two.sided less conservative
                 interval = 0.95, #confidence interval to calculate
                 speedup_parallel = TRUE
)
{
  if(speedup_parallel) #3x faster
  {
    cl = parallel::makePSOCKcluster(parallel::detectCores()-1)
    parallel::clusterExport(cl = cl, 
                            varlist = c('mean.circular',
                                        'circular',
                                        'rvonmises'),
                            envir = .GlobalEnv
    )
    parallel::clusterExport(cl = cl, 
                            varlist = c('MeanRvm',
                                        'angles',
                                        'm1',
                                        'k1',
                                        'm2',
                                        'k2',
                                        'w1',
                                        'n',
                                        'au',
                                        'ar'),
                            envir = environment()
    )
    #simulate primary mean
    m1_est = 
      parallel::parSapply(cl = cl,
                          X = 1:n,
                          FUN = function(i)
                          {
                            eval.parent(
                              {
                                MeanRvm(n = round(length(angles)*w1), #estimate number of observations at primary mean
                                        mu = m1, 
                                        kappa = k1,
                                        au = au,
                                        ar = ar)
                              }
                            )
                          },
                          simplify = 'array' #return an array of simulated angles
      )
    if(!is.na(m2)) #if there is a valid secondary mean
    {
      m2_est = 
        parallel::parSapply(cl = cl,
                            X = 1:n,
                            FUN = function(i)
                            {
                              eval.parent(
                                {
                                  MeanRvm(n = round(length(angles)*(1-w1)), #estimate number of observations at secondary mean
                                          mu = m2, 
                                          kappa = k2,
                                          au = au,
                                          ar = ar)
                                }
                              )
                            },
                            simplify = 'array' #return an array of simulated angles
        )
    }
    parallel::stopCluster(cl)
  }else
  { #if not using parallel, use the slower version via replicate()
    m1_est = replicate(n = n, 
                       MeanRvm(n = round(length(angles)*w1), 
                               mu = m1, 
                               kappa = k1,
                               au = au,
                               ar = ar)
    )
    if(!is.na(m2))
    {
      m2_est = replicate(n = n, 
                         MeanRvm(n = round(length(angles)*(1-w1)), 
                                 mu = m2, 
                                 kappa = k2,
                                 au = au,
                                 ar = ar)
      )
    }
  }
  return(
    if(calc_q) #calculate quantiles only if requested
    {
      #either two-sided, symmetrical around mean change
      #or one-sided, from zero change towards mean change
      probs1 = switch(alternative,
                      two.sided = sort(c(c(0,1)+c(1,-1)*(1-interval)/2, 0.5)),
                      one.sided = sort(c(c(0,1)+
                                           (if(Mod360.180(m1)>0) #N.B. quantile.circular counts anticlockwise
                                           {c(1,0)}else
                                           {c(0,-1)}
                                           )*(1-interval), 0.5)),
                      sort(c(c(0,1)+ #default to one-sided
                               (if(Mod360.180(m1)>0)
                               {c(1,0)}else
                               {c(0,-1)}
                               )*(1-interval), 0.5))
      )
      if(is.na(m2))
      {
        if(force_mu)
        {
          Mod360.180(
            quantile( Mod360.180(as.numeric(m1_est) - m1),
                      probs = probs1) + m1
          )
        }else
        {
          Mod360.180(
            quantile.circular(x = circular(x = m1_est,
                                           units = au,
                                           rotation = ar),
                              probs = probs1)
          )
        }
      }else
      {
        probs2 = switch(alternative,
                        two.sided = sort(c(c(0,1)+c(1,-1)*(1-interval)/2, 0.5)),
                        one.sided = sort(c(c(0,1)+
                                             (if(Mod360.180(m2)>0)
                                             {c(1,0)}else
                                             {c(0,-1)}
                                             )*(1-interval), 0.5)),
                        sort(c(c(0,1)+ #default to one-sided
                                 (if(Mod360.180(m2)<0)
                                 {c(1,0)}else
                                 {c(0,-1)}
                                 )*(1-interval), 0.5))
        )
        list(m1 = if(force_mu)
        {
          Mod360.180(
            quantile( Mod360.180(as.numeric(m1_est) - m1),
                      probs = probs1) + m1
          )
        }else
        {
          Mod360.180(
            quantile.circular(x = circular(x = m1_est,
                                           units = au,
                                           rotation = ar),
                              probs = probs1)
          )
        },
        m2 = if(force_mu)
        {
          Mod360.180(
            quantile( Mod360.180(as.numeric(m2_est) - m2),
                      probs = probs2) + m2
          )
        }else
        {
          Mod360.180(
            quantile.circular(x = circular(x = m2_est,
                                           units = au,
                                           rotation = ar),
                              probs = probs2)
          )
        }
        )
      }
    }else
    { #if quantiles not requested, return the simulations (mainly for troubleshooting)
      if(is.na(m2))
      {
        m1_est = 
          sapply(X = m1_est, FUN = Mod360.180)
      }else
      {
        list(
          m1_est = 
            sapply(X = m1_est, FUN = Mod360.180),     
          m2_est = 
            sapply(X = m2_est, FUN = Mod360.180),
        )
      }
    }
  )
}


PlotCI_vM = function(ci_vec,
                     col = 'salmon',
                     lwd = 2,
                     radius = 0.95,
                     ...)#passed to lines()
{
  ci_vec = as.numeric(ci_vec)#remove circular formatting!
  #changed on 20250815, plotting issues near median
  angle_seq1.1 = 
    seq(from = ci_vec[1], #lower
        to = ci_vec[1] +
          Mod360.180(ci_vec[2]-ci_vec[1]), #median
        length.out =1e2/2)
  angle_seq1.2 = 
    seq(from = ci_vec[2], #median
        to = ci_vec[2] +
          Mod360.180(ci_vec[3]-ci_vec[2]) , #upper
        length.out =1e2/2)
  lines(x = radius*sin( rad(angle_seq1.1) ),
        y = radius*cos( rad(angle_seq1.1) ),
        col = col,
        lwd = lwd,
        lend = 'butt',
        ...
  )
  lines(x = radius*sin( rad(angle_seq1.2) ),
        y = radius*cos( rad(angle_seq1.2) ),
        col = col,
        lwd = lwd,
        lend = 'butt',
        ...
  )
  if(!is.na(ci_vec[4]))
  {
    #changed on 20250815
    angle_seq2.1 = 
      seq(from = ci_vec[1+3],
          to = ci_vec[1+3] +
            Mod360.180(ci_vec[2+3]-ci_vec[1+3]),
          length.out =1e2/2)
    
    angle_seq2.2 = 
      seq(from = ci_vec[2+3],
          to = ci_vec[2+3] +
            Mod360.180(ci_vec[3+3]-ci_vec[2+3]) ,
          length.out =1e2/2)
    lines(x = radius*sin( rad(angle_seq2.1) ),
          y = radius*cos( rad(angle_seq2.1) ),
          col = col,
          lwd = lwd,
          lend = 'butt',
          ....)
    lines(x = radius*sin( rad(angle_seq2.1) ),
          y = radius*cos( rad(angle_seq2.1) ),
          col = col,
          lwd = lwd,
          lend = 'butt',
          ...)
  }
}


DescriptCplot = function(k,
                         m = 0,#in degrees
                         ndata = 20,
                         lw = 3,
                         pcol = col_obs,
                         mvcol = col_rho,
                         sdcol = col_sd,
                         cicol = col_rho,
                         denscol = col_pdf,
                         bins = 360/5-1,
                         stack = TRUE,
                         sep = 0.05,
                         shrink = 1.25,
                         refline = NA,
                         save_sample = FALSE,
                         seed = 20250815,
                         ...#passed to points.circular
)
{
  #convert mu to circular class
  cm = circular(x = m,
                units = 'degrees',
                rotation = 'clock',
                zero = pi/2)
  #kappa = 0 v. random, make repeatable
  set.seed(seed)
  #generate dataset
  cd = rvonmises(mu = cm,
                 kappa = k,
                 n = ndata)
  #Generate a vector of angles
  ra = seq(from = -pi, to = pi, by = 0.01)
  #plot dataset
  plot.circular(x = circular(NA,
                             units = 'degrees',
                             rotation = 'clock',
                             zero = pi/2),
                stack = stack,
                bins = bins,
                shrink = shrink,
                sep = sep,
                col = pcol,
                axes = FALSE
  )
  #add reference line
  if(!is.na(refline))
  {
    arrows.circular(x = circular(refline,
                                 units = 'degrees',
                                 rotation = 'clock',
                                 zero = pi/2),
                    col = 'gray',
                    length = 0,
                    lend = 'butt',
                    lwd = lw)
  }
  #add probability density
  lines.circular(x = ra,
                 y = dvonmises(x = circular(ra,
                                            units = 'radians',
                                            rotation = 'clock',
                                            zero = pi),
                               mu = cm,
                               kappa = k)/shrink,
                 col = denscol,
                 lwd = lw)
  points.circular(cd,
                  stack = stack,
                  bins = bins,
                  shrink = shrink,
                  sep = sep,
                  col = pcol,
                  ...
  )
  #add the mean vector
  if(k>0)
  {
    arrows.circular(x = cm,
                    y = A1(k),
                    col = mvcol,
                    lwd = lw,
                    length = 0.1/shrink#,shrink = 1/shrink
    )
  }else
  {
    points(x = 0,
           y = 0, 
           col = mvcol,
           pch = 19,
           lwd = lw)
  }
  #add SD as lines
  #estimate Mardia SD
  msd = MardiaSD(k = k)
  arrows.circular(x = cm+circular(deg(msd),
                                  units = 'degrees',
                                  rotation = 'clock',
                                  zero = pi/2),
                  y = 1+5*sep*shrink,
                  col = sdcol,
                  lwd = lw,
                  lty = 3,
                  length = 0)
  arrows.circular(x = cm-circular(deg(msd),
                                  units = 'degrees',
                                  rotation = 'clock',
                                  zero = pi/2),
                  y = 1+5*sep*shrink,
                  col = sdcol,
                  lwd = lw,
                  lty = 3,
                  length = 0)
  # lines.circular(x = circular(rep(0-deg(msd), 2),
  #                             units = 'degrees',
  #                             rotation = 'clock',
  #                             zero = pi/2),
  #                 y = c(2,3)*sep*shrink,
  #                col = sdcol,
  #                lwd = lw)
  #add CI of the mean
  #calculate vector of estimates
  civ = CI_vM(angles = cd,
              m1 = cm,
              k1 = k,
              alternative = 'two.sided',
              force_mu = if(k == 0 ){TRUE}else{FALSE})
  #plot quantiles of estimates
  PlotCI_vM(ci_vec = civ,
            col = cicol, lwd = lw,
            radius = 1+5*sep*shrink)
  #add a title
  mtext(text = paste0('κ = ', k),
        side = 1,
        line = -3)
  if(save_sample)
  {return(cd)}
}

PCfun = function(angles,
                 col,
                 shrink = 1.5,
                 title = '',
                 plot_rho = TRUE,
                 side = 1,
                 ...)
{
  ca = circular(x = angles,
                units = 'degrees',
                rotation = 'clock')
  plot.circular(x = ca,
                col = col,
                stack = TRUE,
                bins = 355/5,
                units = 'degrees',
                rotation = 'clock',
                zero = pi/2,
                shrink = shrink,
                ...)
  mtext(text = title,
        side = side,
        line = -2)
  lines(x = c(0,0),
        y = c(-1,1),
        col = 'gray')
  if(plot_rho)
  {
  arrows.circular(x = mean.circular(ca),
                  y = rho.circular(ca),
                  zero = pi/2,
                  rotation = 'clock',
                  col = col,
                  length =0.1)
  }
}

#histograms on a vertical axis
#any data, but plotted as a histogram on a vertical rather than horizontal axis
VertHist = function(data, # numerical data vector
                    breaks = 1e2,
                    ylab = 'data',
                    xlab = 'density',
                    ylim = NULL,
                    main = '',
                    col = 'gray',
                    border = NA,
                    ...)
{
  hst = hist(x = data, # calculate the histogram but don't plot it
             breaks = breaks, # user defined breaks
             plot = FALSE)
  with(hst,
       {
         plot(x = NULL, #open an empty plot
              xlim = c(0, max(density)),
              ylim = if(is.null(ylim)){range(mids)}else{ylim},
              xlab = xlab,
              ylab = ylab,
              main = main)
         #plot each bar
         for(i in 1:length(mids))
         {
           rect(xleft = 0,
                xright = density[i],
                ybottom = breaks[i], 
                ytop = breaks[i + 1],
                col = col,
                border = border,
                ...
           )
         }
       }
  )
}


#invert the softplus link
#https://en.wikipedia.org/wiki/Softplus
#we are using this as our _inverse_ link function for kappa,
#maps almost 1:1 but keeps values >0 for low estimates
softplus = function(x)
{
  log(exp(x)+1) 
}
#this would return our kappa estimates back to the original scale
inv_softplus = function(x)
{
  log(exp(x)-1) 
}

#function for credible interval of the intercept
CI_unwrap = function(data,
                     formula = bf(y~1),
                     est_kappa = TRUE,
                     predictors = TRUE,
                     return_model = FALSE,
                     interval = 0.95,
                     prior = NULL,
                     ...#passed to brms
)
{
  if(is.null(prior))
  {
    if(predictors)
    {
      prior =  prior('normal(0,pi())', class = 'Intercept', dpar = 'mu') +
        prior('normal(0,pi()/3)', class = 'b', dpar = 'mu') +
        prior('normal(3,3)', class = 'Intercept', dpar = 'kappa')
      prior('normal(0,3)', class = 'b', dpar = 'kappa')
    }else{
      prior =  prior('normal(0,pi())', class = 'Intercept', dpar = 'mu') +
        prior('normal(3,3)', class = 'kappa')
    }
  }
  #set up required inverse link
  softplus = function(x)
  {
    log(exp(x)+1) 
  }
  #set up required Stan functions
  mod_circular_fun = stanvar(scode = "
    real mod_circular(real y) {
      return fmod(y + pi(), 2*pi()) - pi();
    }
  ",
                             block = 'functions')
  unwrap_von_mises = custom_family(
    "unwrap_von_mises", dpars = c("mu", "kappa"),
    links = c('identity',#brms cannot accept custom link functions, do via nl instead
              "softplus"), 
    lb = c(-pi, 0), ub = c(pi, NA),
    type = "real",
  )
  
  stan_unwrap_fun = stanvar(scode = "
    real unwrap_von_mises_lpdf(real y, real mu, real kappa) {
      return von_mises_lpdf(y | mod_circular(mu), kappa);
    }
    real unwrap_von_mises_rng(real mu, real kappa) {
      return von_mises_rng( mod_circular(mu) , kappa);
    }
  ",
                            block = 'functions') 
  
  bmod = brm(formula = formula,
             data = data,
             family = unwrap_von_mises,
             stanvars = stan_unwrap_fun + mod_circular_fun,
             prior = prior,
             ...
  )
  all_draws = as_draws_df(bmod)
  # mu_draws = as_draws_df(bmod, 
  #                        variable = 'Intercept',
  #                        dpar = 'mu')
  pbs =  c(0,0.5,1)+
    c(1,0,-1)*
    (1-interval)*0.5
  mu_q = with(all_draws,
              {
                quantile.circular(x = circular(Intercept,
                                               template = 'none'),
                                  probs = pbs)
              }
  )
  if(predictors)
  {
    mux_q = with(all_draws,
                 {
                   quantile.circular(x = circular(Intercept+b_x,
                                                  template = 'none'),
                                     probs = pbs)
                 }
    )
  }
  
  # kappa_draws = as_draws_df(bmod, 
  #                           variable = 'kappa',
  #                           regex = TRUE)
  if(est_kappa)
  {
    if(predictors)
    {
      kappa_q = with(all_draws,
                     {
                       softplus(x = quantile(x = Intercept_kappa,
                                             probs = pbs)
                       )
                     }
      )
      kappax_q = with(all_draws,
                      {
                        quantile(x = softplus(Intercept_kappa + b_kappa_x),
                                 probs = pbs)
                      }
      )
    }else{
      kappa_q = with(all_draws,
                     {
                       softplus(x = quantile(x = kappa,
                                             probs = pbs)
                       )
                     }
      )
    }
  }
  
  if(predictors)
  {
    rlst = list(mu = mu_q,
                mu_x = mux_q)
    rlst$kappa = if(est_kappa){kappa_q}else{NULL}
    rlst$kappa_x = if(est_kappa){kappax_q}else{NULL}
    rlst$model = if(return_model){bmod}else{NULL}
  }else{
    rlst = list(mu = mu_q)
    rlst$kappa = if(est_kappa){kappa_q}else{NULL}
    rlst$model = if(return_model){bmod}else{NULL}
  }
  return(rlst)
}


# Divergence from home direction -------------------------------------------
set.seed(0120810506)#ISBN Batschelet, 1981

par(pty = 's')
par(mar = c(0,0,0,0))

cd_divergence = DescriptCplot(m = -15,
                              k = 10,
                              refline = 0, 
                              ndata = 20,
                              sdcol = NA,
                              denscol = NA,
                              save_sample = TRUE)

# circular zero
c0 = circular(x = 0,
              units = 'degrees',
              rotation = 'clock',
              zero = pi/2)

#v-test finds significant orientation in 0°, but we know the true direction is different
rayleigh.test(cd_divergence, mu = c0)
ci_divergence = CI_vM(angles = cd_divergence,
                      m1 = -15,
                      k1 = 10,
                      alternative = 'two.sided')
mtext(text = paste0('(',paste(signif(ci_divergence[-2], 2), collapse = ' '), ')'),
      side = 1,
      line = -1)


## model version --------------------------------------------------------

#fit a generic unwrap model
ci_uw = CI_unwrap(data = data.frame(y = rad(cd_divergence)), 
                  est_kappa = TRUE,
                  return_model = TRUE,
                  predictors = FALSE,
                  backend = 'cmdstan'#faster and more reliable
)

draws_divergence = with(ci_uw, as_draws_df(model))


par(pty = 's')
par(mar = c(0,0,0,0),
    mfrow = c(1,2))
PCfun(cd_divergence,
      col = col_obs,
      sep = 0.05,
      shrink = 1.25,
      plot_rho = FALSE)
arrows.circular(x = circular(-15,
                             units = 'degrees',
                             rotation = 'clock',
                             zero = pi/2),
                y = A1(10),
                col = col_rho,
                lwd = 5,
                length = 0.1/1.25
)
with(draws_divergence,
     points(x = sin(Intercept)*A1(softplus(kappa)),
            y = cos(Intercept)*A1(softplus(kappa)),
            col = adjustcolor(col = col_sd,
                              alpha.f = 1/255),
            pch = 19)
            )
with(ci_uw,
     arrows.circular(x = circular(mu[2],
                                  rotation = 'clock',
                                  zero = pi/2),
                     y = A1(kappa[2]),
                     lwd = 2,
                     length = 0.1/1.25,
                     col = adjustcolor(col_sd, alpha.f = 200/255))
)
with(draws_divergence,
VertHist(data = deg(Intercept), 
         main = 'mean angle',
         ylim = c(-30, 15),
         col = adjustcolor(col_sd, alpha.f = 100/255),
         cex.axis = 0.7))
abline(h = 0,
       col = 'gray',
       lwd = 7)

# Reduced concentration -------------------
set.seed(0120810506)#ISBN Batschelet, 1981

par(pty = 's')
par(mar = c(0,0,0,0))
par(mfrow = c(1,2))

cd_3 = DescriptCplot(m = 0,
                     k = 3,
                     ndata = 20,
                     sdcol = NA,
                     denscol = NA,
                     refline = c0,
                     save_sample = TRUE)

cd_0.5 = DescriptCplot(m = 0,
                       k = 0.5,
                       ndata = 20,
                       sdcol = NA,
                       denscol = NA,
                       refline = c0,
                       save_sample = TRUE)

## Model version ---------------------------------------------------------


ci_kappa = CI_unwrap(data = data.frame(y = rad(c(cd_3, cd_0.5)),
                                       x = c(rep(0, length(cd_3)),
                                             rep(1, length(cd_0.5)))
                                      ),
                    formula = bf(y~x,
                                 kappa~x),
                    est_kappa = TRUE,
                    return_model = TRUE,
                    backend = 'cmdstan'#faster and more reliable
                    )

draws_kappa = with(ci_kappa, as_draws_df(model))


par(pty = 's')
par(mar = c(0,0,0,0),
    mfrow = c(1,3))
PCfun(cd_3,
      col = col_obs,
      sep = 0.05,
      shrink = 1.25,
      plot_rho = FALSE)
arrows.circular(x = circular(0,
                             units = 'degrees',
                             rotation = 'clock',
                             zero = pi/2),
                y = A1(3),
                col = col_rho,
                lwd = 5,
                length = 0.1/1.25
)
with(draws_kappa,
     points(x = sin(Intercept)*A1(softplus(Intercept_kappa)),
            y = cos(Intercept)*A1(softplus(Intercept_kappa)),
            col = adjustcolor(col = col_sd,
                              alpha.f = 1/255),
            pch = 19)
)
with(ci_kappa,
     arrows.circular(x = circular(mu[2],
                                  rotation = 'clock',
                                  zero = pi/2),
                     y = A1(kappa[2]),
                     lwd = 2,
                     length = 0.1/1.25,
                     col = adjustcolor(col_sd, alpha.f = 200/255))
)
PCfun(cd_0.5,
      col = col_obs,
      sep = 0.05,
      shrink = 1.25,
      plot_rho = FALSE)
arrows.circular(x = circular(0,
                             units = 'degrees',
                             rotation = 'clock',
                             zero = pi/2),
                y = A1(0.5),
                col = col_rho,
                lwd = 5,
                length = 0.1/1.25
)
with(draws_kappa,
     points(x = sin(Intercept+b_x)*A1(softplus(Intercept_kappa+b_kappa_x)),
            y = cos(Intercept+b_x)*A1(softplus(Intercept_kappa+b_kappa_x)),
            col = adjustcolor(col = col_sd,
                              alpha.f = 1/255),
            pch = 19)
)
with(ci_kappa,
     arrows.circular(x = circular(mu_x[2],
                                  rotation = 'clock',
                                  zero = pi/2),
                     y = A1(kappa_x[2]),
                     lwd = 2,
                     length = 0.1/1.25,
                     col = adjustcolor(col_sd, alpha.f = 200/255))
)
with(draws_kappa,
     VertHist(data = b_kappa_x, 
              main = 'change in kappa',
              ylim = c(-10, 1),
              col = adjustcolor(col_sd, alpha.f = 100/255),
              cex.axis = 0.7))
abline(h = 0,
       col = 'gray',
       lwd = 7)

with(draws_kappa, mean(b_kappa_x < 0)*100 ) #nearly all estimates suggest a reduction in kappa

watson.two.test(cd_3, cd_0.5)#no difference detected
