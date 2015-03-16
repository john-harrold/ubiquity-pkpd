function [simulation_options]=process_covariates(cfg, simulation_options)

% checking to make sure there are covariates
if(isfield(cfg.options, 'inputs'))
if(isfield(cfg.options.inputs, 'covariates'))
for cov_idx = 1:length(cfg.options.inputs.covariate_names)

  mycov_str  = cfg.options.inputs.covariate_names{cov_idx};

  % pulling the covariate out of cfg
  mycov = getfield(cfg.options.inputs.covariates, mycov_str);

  % constructing time varying inputs from
  % the covariates 
  if(strcmp(mycov.cv_type, 'linear'))
    simulation_options.mytimevarying(cov_idx).name        = mycov_str;          
    simulation_options.mytimevarying(cov_idx).value(1,:)  = mycov.times.values;
    simulation_options.mytimevarying(cov_idx).value(2,:)  = mycov.values.values;
  elseif(strcmp(mycov.cv_type, 'step'))
    simulation_options.mytimevarying(cov_idx).name        = mycov_str;          
    simulation_options.mytimevarying(cov_idx).value       = make_step([mycov.times.values; mycov.values.values], simulation_options.output_times)';
  else
    vp(cfg, sprintf('Unknown covariate type (%s) for covariate (%s)', mycov.cv_type, mycov_str));
  end
end
end
end

function [step_val]=make_step(levels,output_times)

switch_times  = levels(1,:);
magnitude     = levels(2,:);

delta         = 500*eps;

for i = 1:numel(switch_times)
    if (1==i)
      step_val = [switch_times(i) magnitude(i)];
    else
      step_val = [step_val
                   switch_times(i)          magnitude(i-1)
                  (switch_times(i)+delta)   magnitude(i)];
    end
end

% extending the last infusion rate to the end of 
% simulation time
if(max(output_times) > max(step_val(:,1)))
  step_val = [step_val
              max(output_times) step_val(end,2)];
end
