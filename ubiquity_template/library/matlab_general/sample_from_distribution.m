function [my_parameters] = sample_from_distribution(my_parameters)
%
%  function [my_parameters] = sample_from_distribution(my_parameters)
%
%  For a list of parameters with a mean, standard deviation, and distribution
%  type (normal or lognormal), a random sample is generated. This is used when
%  running monty carlo simulations. The information about the parameter
%  distributions is specified in the distribution_information field:
% 
%  my_parameters.distribution_information =  ...
%          [ {'CL'}         .2       .03          {'lognormal'}
%            {'Vc'}         50       9.4          {'normal'} ];
%           % parameter    mean     standard      distribution
%           % name         value    deviation     type
%                  
%  The first column represents the parameter name, the second is the mean
%  followed by the standard deviation (third). The fourth column specifies the
%  type of distribution ('lognormal' or 'normal') this applies to.
%                  
%  The field 'my_sample' is created. This is a column vector where each row
%  corresponds to the row in the distribution_information cell array. For the
%  example above it may look something like:
%  my_parameters.my_sample = [.1
%                              53.2];
%                  


[num_parameters, ncols] = size(my_parameters.distribution_information);

my_parameters.my_sample = [];


for parameter_idx=1:num_parameters
  % sigma = standard deviation
  % variance = standard_deviation ^2
  tmp_mean  = cell2mat( my_parameters.distribution_information(parameter_idx, 2));
  tmp_sigma = cell2mat( my_parameters.distribution_information(parameter_idx, 3));
  tmp_var   = tmp_sigma.^2;
  if(strcmp('normal', char(my_parameters.distribution_information(parameter_idx, 4))))
    my_parameters.my_sample = [my_parameters.my_sample; normrnd(tmp_mean, tmp_sigma)];
  elseif(strcmp('lognormal', char(my_parameters.distribution_information(parameter_idx, 4))))
    % 
    % First we have to log transform the mean and standard deviation
    % 
    lt_mean =  log((tmp_mean.^2)./sqrt(tmp_var+tmp_mean.^2));
    lt_sigma= sqrt(log(tmp_var./(tmp_mean.^2)+1));
    % 
    % Now we sample from the log transformed distribution 
    % 
    my_parameters.my_sample = [my_parameters.my_sample; lognrnd(lt_mean , lt_sigma)];
  else
      disp('Warning: the distribution type %s was not recognized',  char(my_parameters.distribution_information(parameter_idx, 4)));
  end
end
