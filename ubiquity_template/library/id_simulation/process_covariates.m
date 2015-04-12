function [simulation_options]=process_covariates(cfg, simulation_options)

% checking to make sure there are covariates
if(isfield(cfg.options, 'inputs'))
if(isfield(cfg.options.inputs, 'covariates'))
tvcounter = 1;
for cov_idx = 1:length(cfg.options.inputs.covariate_names)

  mycov_str  = cfg.options.inputs.covariate_names{cov_idx};

  % pulling the covariate out of cfg
  mycov = getfield(cfg.options.inputs.covariates, mycov_str);

  % constructing time varying inputs from
  % the covariates 
  if(strcmp(mycov.cv_type, 'linear'))
    simulation_options.mytimevarying(tvcounter).name        = mycov_str;          
    simulation_options.mytimevarying(tvcounter).value(1,:)  = mycov.times.values;
    simulation_options.mytimevarying(tvcounter).value(2,:)  = mycov.values.values;
  elseif(strcmp(mycov.cv_type, 'step'))
    simulation_options.mytimevarying(tvcounter).name        = mycov_str;          
    simulation_options.mytimevarying(tvcounter).value       = make_step([mycov.times.values; mycov.values.values], simulation_options.output_times)';
  else
    vp(cfg, sprintf('Unknown covariate type (%s) for covariate (%s)', mycov.cv_type, mycov_str));
  end
  tvcounter = tvcounter +1;

  % 
  %  Next we create an initial condition for the covariate so it can be used
  %  in static secondary parameters 
  % 
  simulation_options.mytimevarying(tvcounter).name        = sprintf('SIMINT_CVIC_%s', mycov_str);          
  simulation_options.mytimevarying(tvcounter).value       = make_step([mycov.times.values(1); mycov.values.values(1)], simulation_options.output_times)';
  tvcounter = tvcounter +1;
end
end
end

function [step_val]=make_step(levels,output_times)

switch_times  = levels(1,:);
magnitude     = levels(2,:);


for i = 1:numel(switch_times)
    if (1==i)
      step_val = [switch_times(i) magnitude(i)];
    else
      if(switch_times(i) == 0)
        delta         = 250*eps;
      else
        delta         = 250*eps*switch_times(i);
      end
      step_val = [step_val
                  (switch_times(i)-delta)   magnitude(i-1)
                  (switch_times(i)+delta)   magnitude(i)];
    end
end

% extending the last infusion rate to the end of 
% simulation time
if(max(output_times) > max(step_val(:,1)))
  step_val = [step_val
              max(output_times) step_val(end,2)];
end
