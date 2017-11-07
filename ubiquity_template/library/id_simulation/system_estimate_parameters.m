function [pest]=system_estimate_parameters(cfg, flowctl, analysis_name)
% function [pest]=system_estimate_parameters(cfg, flowctl, analysis_name)
% system_estimate_parameters - controls the estimation process
% flowctl - controls the flow of the estimation process with the
%           following options:
%           'plot previous estimate'
%           'previous estimate as guess'
%           'estimate'
%           'plot guess'
% analysis_name - name used for archinving results


if(strcmp(flowctl , 'estimate') | strcmp(flowctl, 'previous estimate as guess'))
  % checking the analysis_name
  name_check = ubiquity_name_check(analysis_name);
  if(~name_check.isgood)
    vp(cfg, sprintf('Error: the analyssis name >%s< is invalid', analysis_name  ))
    vp(cfg, sprintf('Problems: %s', name_check.msg))
    analysis_name   = 'analysis';
    vp(cfg, sprintf('Instead Using: %s', analysis_name  ))
  end
  % Loading the previous estimate
  if(strcmp(flowctl, 'previous estimate as guess'))
    eval(sprintf('load output%sanalysis_%s.mat pest',filesep, analysis_name))
    cfg.estimation.parameters.guess = pest;
  end

  % Performing the parameter estimation: 
  pest = estimate_parameters(cfg);

  % Saving the results:
  eval(sprintf('save output%sanalysis_%s.mat pest',filesep, analysis_name))
  archive_estimation(analysis_name, cfg);
elseif(strcmp(flowctl , 'plot guess'))
  pest   = system_fetch_guess(cfg);
else 
  % Loading the previously saved results:
  eval(sprintf('load output%sanalysis_%s.mat pest',filesep, analysis_name))
end

