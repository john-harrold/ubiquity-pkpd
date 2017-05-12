function [cfg]=system_set_covariate(cfg, covariate, times, values)
% [cfg]=system_set_covariate(cfg, covariate, times, values)
%
%  covariate = string containing the covariate name
%  times = vector of times when the covariate changes
%  values = corresponding vector of values
%
  isgood = true;
  if((length(times) ~= length(values)))
    vp(cfg, 'The times and values have differnt lengths');
    isgood = false;
  end

  if(sum(strcmp(cfg.options.inputs.covariate_names, covariate)) ~= 1)
    vp(cfg, sprintf('The covariate name %s could not be found', covariate));
    isgood = false;
  end
  
  if(isgood)
    eval(sprintf('cfg.options.inputs.covariates.%s.times.values  = times;', covariate));
    eval(sprintf('cfg.options.inputs.covariates.%s.values.values = values;', covariate));
  else
    vp(cfg, 'Something went wrong and the covariate')
    vp(cfg, 'was not set, see the messages above.');
  end
