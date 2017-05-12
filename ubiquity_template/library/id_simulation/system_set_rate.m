function [cfg]=system_set_rate(cfg, rate, times, levels)
% [cfg]=system_set_rate(cfg, rate, times, levels)
%
%  rate = string containing the rate name
%  times = vector of times when the rate changes
%  levels = corresponding vector of infusion rates
%
  isgood = true;
  if((length(times) ~= length(levels)))
    vp(cfg, 'The times and levels have different lengths');
    isgood = false;
  end

  if(sum(strcmp(cfg.options.inputs.infusion_rate_names, rate)) ~= 1)
    vp(cfg, sprintf('The rate name %s could not be found', rate));
    isgood = false;
  end
  
  if(isgood)
    eval(sprintf('cfg.options.inputs.infusion_rates.%s.times.values  = times;', rate));
    eval(sprintf('cfg.options.inputs.infusion_rates.%s.levels.values = levels;', rate));
  else
    vp(cfg, 'Something went wrong and the infusion rate')
    vp(cfg, 'was not set, see the messages above.');
  end
