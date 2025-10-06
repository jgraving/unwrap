
#  Load packages ----------------------------------------------------------
#package for circular data
require(circular)
#package for Bayesian modelling in Stan
require(cmdstanr)
#package for building Stan models
require(brms)

# Set up functions --------------------------------------------------------
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

#the radian equivalent
mod_circular = function(x)
{
  atan2(y = sin(x),
        x = cos(x))
}

#the unwrapped (no discontinuities)
unwrap_circular = function(x)
{
  mux = mean.circular(x = circular(x = x, template = 'none'))
  centx = atan2(y = sin(x - mux),
                x = cos(x  - mux))
  unwrx = centx + mux
}

#degree version
unwrap_circular_deg = function(x)
{
  mux = mean.circular(x = circular(x = x, template = 'none'))
  centx = atan2(y = sin(x - mux),
                x = cos(x  - mux))
  unwrx = centx + mux
  return(deg(unwrx))
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
                         calc_ci = FALSE,
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
  if(!is.na(denscol))
  {
  lines.circular(x = ra,
                 y = dvonmises(x = circular(ra,
                                            units = 'radians',
                                            rotation = 'clock',
                                            zero = pi),
                               mu = cm,
                               kappa = k)/shrink,
                 col = denscol,
                 lwd = lw)
  }
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
  if(!is.na(sdcol))
  {
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
  }
  if(calc_ci)
  {
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
  }
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


## Modelling functions ---------------------------------------------------


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
                    axes = TRUE,
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
              main = main,
              axes = axes)
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


#add contours
Draws2Cont = function(draws,
                      palette = 'Heat 2',
                      nlevels = 20,
                      x_string = 'sin(Intercept)*A1(softplus(Intercept_kappa))',
                      y_string = 'cos(Intercept)*A1(softplus(Intercept_kappa))',
                      alpha = 200/255,
                      ngrid = 25, # defaults to a 25x25 grid
                      cropc = FALSE, #crop region outside circle
                      denstype = 'relative' # 'normalised' or 'relative' (normalised fails to plot low densities)
)
{
  kdc = with(draws,
             {
               MASS::kde2d(x = eval(str2lang(x_string)),
                           y = eval(str2lang(y_string)),
                           n = ngrid)
             }
  )
  if(cropc)
  {
    xy = with(kdc, expand.grid(x = x,y= y))#find coordinates of z variable
    idc = with(xy, x^2+y^2 < 1.0)
    #crop edge of circle
    kdc = within(kdc,
                 {z[!idc] = 0})
  }
  with(kdc,
       {
         .filled.contour(x = x,
                         y = y,
                         z = z,
                         levels = (1:nlevels)*
                           switch(EXPR = denstype,
                                  relative = max(z/nlevels),
                                  normalised = sum(z/nlevels),#warning, fails to plot low densities
                                  max(z/nlevels)
                           ),
                         col = hcl.colors(n = nlevels,
                                          palette = palette,
                                          rev = TRUE,
                                          alpha = alpha)
         )
       }
  )
}

# Draws2CircCont = function(draws,
#                           palette = 'Heat 2',
#                           nlevels = 20,
#                           x_string = 'Intercept',
#                           y_string = 'A1(softplus(kappa))',
#                           alpha = 200/255,
#                           denstype = 'relative' # 'normalised' or 'relative'
# )
# {
#   with(draws,
#        {
#          xx = eval(str2lang(x_string))
#          yy = eval(str2lang(y_string))
#          xtile = c(xx - 2*pi, xx, xx + 2*pi)
#          ytile = c(yy, yy, yy)
#          kdc = MASS::kde2d(x = xtile,
#                            y = ytile, 
#                            n = 30)
#          ln = length(kdc$x)/3
#          midspan = ln + 1:ln #indices of the middle span
#          midspanx = with(kdc,order(sin(x[midspan])*y[midspan]) )#indices of the middle span
#          midspany = with(kdc, order(cos(x[midspan])*y[midspan]) ) #indices of the middle span
#          with(kdc,
#               {
#                 .filled.contour(x = sin(x[midspanx])*y[midspanx],
#                                 y = cos(x[midspany])*y[midspany],
#                                 z = z[midspanx, midspany],
#                                 levels = (1:nlevels)*
#                                             switch(EXPR = denstype,
#                                                     relative = max(z[midspanx, midspany]/nlevels),
#                                                     normalised = sum(z[midspanx, midspany])/nlevels,
#                                                     relative = max(z[midspanx, midspany]/nlevels)
#                                                    ),
#                                 col = hcl.colors(n = nlevels,
#                                                  palette = palette,
#                                                  rev = TRUE,
#                                                  alpha = alpha)
#                 )
#               }
#          )
#        }
#   )
# }


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

#TODO replace and archive
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

#calculate rhat for circular variables with extreme ranges
Rhat_unwrap = function(x){rhat(unwrap_circular(x))}
UnwrapRhats = function(uwmod,
                       variable = '^b_zmu',
                       regex = TRUE,
                       digits = 5,
                       ...)
{
  rh =   
    apply(X = as_draws_df(uwmod,
                          variable = variable, 
                          regex = regex,
                          ...),
          MARGIN = 2,
          FUN = Rhat_unwrap)
  nrh = names(rh)
  rh = rh[!( nrh %in% c(".chain", ".iteration", ".draw") )]
  return( round(rh,digits = digits) )
}

#copy of log_lik_von_mises, with helper functions included & removal of circular formatting
log_lik_unwrap_von_mises <- function(i, prep) {
  #remove circular formatting?
  prep$data$Y = as.numeric(prep$data$Y)
  
  args <- list(
    mu = get_dpar(prep, "mu", i),
    kappa = get_dpar(prep, "kappa", i = i)
  )
  
  # ----------- log_lik helper-functions -----------
  # compute (possibly censored) log_lik values
  # @param dist name of a distribution for which the functions
  #   d<dist> (pdf) and p<dist> (cdf) are available
  # @param args additional arguments passed to pdf and cdf
  # @param prep a brmsprep object
  # @return vector of log_lik values
  log_lik_censor <- function(dist, args, i, prep) {
    pdf <- get(paste0("d", dist), mode = "function")
    cdf <- get(paste0("p", dist), mode = "function")
    y <- prep$data$Y[i]
    cens <- prep$data$cens[i]
    if (is.null(cens) || cens == 0) {
      x <- do_call(pdf, c(y, args, log = TRUE))
    } else if (cens == 1) {
      x <- do_call(cdf, c(y, args, lower.tail = FALSE, log.p = TRUE))
    } else if (cens == -1) {
      x <- do_call(cdf, c(y, args, log.p = TRUE))
    } else if (cens == 2) {
      rcens <- prep$data$rcens[i]
      x <- log(do_call(cdf, c(rcens, args)) - do_call(cdf, c(y, args)))
    }
    x
  }
  
  # adjust log_lik in truncated models
  # @param x vector of log_lik values
  # @param cdf a cumulative distribution function
  # @param args arguments passed to cdf
  # @param i observation number
  # @param prep a brmsprep object
  # @return vector of log_lik values
  log_lik_truncate <- function(x, cdf, args, i, prep) {
    lb <- prep$data[["lb"]][i]
    ub <- prep$data[["ub"]][i]
    if (is.null(lb) && is.null(ub)) {
      return(x)
    }
    if (!is.null(lb)) {
      log_cdf_lb <- do_call(cdf, c(lb, args, log.p = TRUE))
    } else {
      log_cdf_lb <- rep(-Inf, length(x))
    }
    if (!is.null(ub)) {
      log_cdf_ub <- do_call(cdf, c(ub, args, log.p = TRUE))
    } else {
      log_cdf_ub <- rep(0, length(x))
    }
    x - log_diff_exp(log_cdf_ub, log_cdf_lb)
  }
  
  # weight log_lik values according to defined weights
  # @param x vector of log_lik values
  # @param i observation number
  # @param prep a brmsprep object
  # @return vector of log_lik values
  log_lik_weight <- function(x, i, prep) {
    weight <- prep$data$weights[i]
    if (!is.null(weight)) {
      x <- x * weight
    }
    x
  }
  
  
  out <- log_lik_censor(
    dist = "von_mises", args = args, i = i, prep = prep
  )
  out <- log_lik_truncate(
    out, cdf = pvon_mises, args = args, i = i, prep = prep
  )
  log_lik_weight(out, i = i, prep = prep)
}
