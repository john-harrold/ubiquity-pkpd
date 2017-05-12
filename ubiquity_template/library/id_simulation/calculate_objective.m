function [value, details]=calculate_objective(parameters,cfg)
%
%  function [value, details]=calculate_objective(parameters,cfg)
%   parameters    -- set of parameters about which to calculate the objective
%   cfg           -- configureation variable passed through to function calls
%                    with the following required fields
%
%   cfg.estimation.observation_function
%       name of matlab function which is called using the following format:
%
%       [observations]=observation_function(parameters,cfg);
%
%       parameters  - is the vector of system and variance parameters
%       cfg         - is a data structure with the following fields specific
%                     to the objective function
%
%       cfg.estimation.objective_type
%           type of objective function to use with the following
%           acceptable values:
%               'wls'    - weighted least squares
%               'ml'     - maximum likelihood
%       
%       The variable 'observations' is the data structure with an field for
%       each subject, and each subject can have multiple outputs 1 ... l. For
%       each subject/output combination a matrix with the following structure
%       needs to be defined:
%       
%       observations.subject1.output1   = [t1   y(t1), yhat(t1), var(t1)
%                                          t2   y(t2), yhat(t2), var(t2)
%                                          .      .       .        .    
%                                          .      .       .        .    
%                                          .      .       .        .    
%                                          tn   y(tn), yhat(tn), var(tn)];
%
%        t    - time 
%        y    - measured value at time (ti)
%        yhat - model prediction at time (ti)
%        var  - variance at time or 1/weight at (ti)
%
%

%
% Initializing outputs
%
value                = 0;
num_measurements.all = 0;
details.value        = 0;
details.outputs      = {};
errorflag = 0;

%
% checking the cfg file to make sure it's complete
%
if(~strcmp('good', check_cfg(cfg)))
    error('There seems to be some problem with the ''cfg'' variable');
end


% bounding the parameters
if(isfield(cfg.estimation.parameters, 'lower_bound') & ...
   isfield(cfg.estimation.parameters, 'upper_bound'))
  [parameters, bperror, objmult] = bound_parameters(parameters, ...
                                    cfg.estimation.parameters.lower_bound, ...
                                    cfg.estimation.parameters.upper_bound);
  if(bperror)
    disp(' -> bound_parameters.m failed');
    disp(' -> ERROR: LB should be less than UB');
    errorflag = 1;
  end
else
  disp(sprintf(' -> unable to retrieve parameter bounds'));
  errorflag = 1;
end


% checking for timeout
% if timeout is reached then an objective function value of infinity is
% returned

if(isfield(cfg, 'estimation'))
if(isfield(cfg.estimation, 'timeout'))
    try
      if( toc > cfg.estimation.timeout)
        value     = inf;
        return;
      end
    catch
     error(' -> To use timeout option, you must also run "tic" to start the timer');
    end
end
end


try 
  % obtaining observation informaiotn by
  % evaluating the user defined function
  eval(sprintf('[observations]=%s(parameters,cfg);',cfg.estimation.observation_function));

  % getting all subjects
  subjects = fieldnames(observations);

catch TRYERROR
  vp(cfg, sprintf(' -> unable to retrieve observations'));
  vp(cfg, sprintf(' -> possible causes:'));
  vp(cfg, sprintf('      o cfg.estimation.observation_function is not defined'));
  vp(cfg, sprintf('      o odd parameter combinations sent to the'));
  vp(cfg, sprintf('        objective function during estimation '));
  vp(cfg, sprintf('        is causing problems '));
  system_log_tryerror(cfg, TRYERROR);
  errorflag = 1;
  value = inf;
end


if(~errorflag)
  try
  %
  % calculating the objective and storing it in 'value'
  %
  
  % looping through each subject
  for  subidx =1:length(subjects)
    % pulling the individual assocaited with subid
    individual = getfield(observations, char(subjects(subidx)));
    % getting that individuals outputs
    outputs  = fieldnames(individual);
    % now looping through each output in that subject
    for  outidx =1:length(outputs)
      output     = getfield(individual, char(outputs(outidx)));
      % calcualting objective function for each case
      %
      % This is taken from section 3.7 in the
      % ADAPT5 Users Guide
      %
      num_measurements.all = num_measurements.all + length(output(:,3));
      yobs  = output(:,2);
      ypred = output(:,3);
      yvar  = output(:,4);
      if(strcmp('wls', cfg.estimation.objective_type))
          component = sum(((yobs-ypred).^2).*(1./yvar));
      elseif(strcmp('ml', cfg.estimation.objective_type))
          component = 1/2.*sum(((yobs-ypred).^2).*(1./yvar)) + sum(log(yvar));
      end
      value   =  value + component;
      % adding component to appropriate output
      % making sure a field exists for this output in the details outputs
      % structure
      if(~isfield(details.outputs, char(outputs(outidx))))
        details.outputs = setfield(details.outputs, char(outputs(outidx)), 0);
      end
        details.outputs = setfield(details.outputs, char(outputs(outidx)), ...
          (component + getfield(details.outputs, char(outputs(outidx)))));
    end
  end
  catch
    value = inf;
  end
end

% adding the constant portion to the negative log likelihood 
% objective though this really shouldn't do anything
if(strcmp('ml', cfg.estimation.objective_type))
  value = value + num_measurements.all*log(2*pi())/2;
end

% Fminsearch does not have an explicit constraint for parameter bounds, so
% these are handled by  increasing the objective function due to the parameter
% bounds being crossed.
if(strcmp(cfg.estimation.optimizer, 'fminsearch'))
  value = value*(1 + objmult);
end

details.value = value;
if(errorflag)
 error(' -> calculate_objective.m failed, see above');
end
