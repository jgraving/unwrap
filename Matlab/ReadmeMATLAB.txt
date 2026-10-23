Native MATLAB examples
----------------------
MatlabOnly.m is a small circular regression example.
CircularDatasetStructures.m analyses the bundled circular_datasets/*.csv
using helper functions in unwrap_functions/.
Both require the Statistics and Machine Learning Toolbox (hmcSampler).
They do not require Stan. Open the chosen script in MATLAB and run it.

Stan example: UnwrapStan.m
-------------------------
The following setup instructions apply only to UnwrapStan.m.
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
