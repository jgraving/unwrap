Before running
make sure you have installed:
 - https://github.com/brian-lau/MatlabStan
 - https://github.com/brian-lau/MatlabProcessManager/
 - cmdstant (recommended https://mc-stan.org/cmdstanr/)

 hard code the cmdstan paths and versions in:
 - stan_home.m (e.g. d = '/Users/username/.cmdstan/cmdstan-2.37.0';)
 - StanModel.m
BY SETTING

       function ver = stan_version(self)
             ver = [2, 37, 0]%for cmdstan-2.37.0
       end

AND
           if isempty(p.Results.stan_version)
             self.stan_version = [2 37 0];
 %             self.stan_version = self.get_stan_version();
          else
             self.stan_version = p.Results.stan_version;
          end

TO CORRESPOND WITH YOUR cmdstan VERSION
