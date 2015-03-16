function [my_noise] = add_noise(my_mean, my_sigma, noise_type);
% function [my_noise] = add_noise(my_mean, my_sigma, noise_type);
% 
%  my_mean    - vector of the mean values 
%  my_sigma   - vector of corresponding standard deviations
%  noise_type - type of noise (either 'normal' or 'log_normal')
% 
%  my_noise   - my_mean + noise
% 



my_sigma = my_sigma.*my_mean;

my_var   = my_sigma.^2;
my_noise = [];

  if(strcmp('normal', noise_type))
    my_noise = normrnd(my_mean, my_sigma);
  elseif(strcmp('log_normal', noise_type))
    % 
    % First we have to log transform the mean and standard deviation
    % 
    %lt_mean = log(my_mean);
    lt_mean = log(my_mean.^2./sqrt(my_var + my_mean.^2));
    lt_sigma= sqrt(log(my_var./(my_mean.^2)+1));

    % 
    % Now we sample from the log transformed distribution 
    % 
    my_noise =  lognrnd(lt_mean , lt_sigma);
  end


