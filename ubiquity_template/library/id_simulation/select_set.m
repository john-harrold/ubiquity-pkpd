function [cfg]=select_set(cfg, set_name, parameter_names)

disp('------------------------------------------------------');
disp('Warning: function select_set has been depreciated');
disp('use system_select_set instead');
disp('------------------------------------------------------');

if(exist('parameter_names', 'var'))
  cfg = system_select_set(cfg, set_name, parameter_names);
else
  cfg = system_select_set(cfg, set_name);
end


% % function [cfg]=select_set(cfg, set_name, parameter_names)
% % 
% % % Used to select a parameter set and define a subset of parameters to be
% % % considered for estimation purposes 
% %
% % % After generating a system you can create the default cfg variable using
% % % the following command:
% %
% % cfg = auto_fetch_system_information
% %
% % % This creates the cfg.parameters data structure with the 'values' set to
% % % the default parameter set. To overwrite these with the set named 'mouse'
% % % (assuming this has been defined in system.txt), you would use the
% % % following:
% %
% % cfg=select_set(cfg, 'mouse');
% %
% % % This overwrites cfg.parameters.values with the values in the parameter
% % % set 'mouse'. It also creates the data structure
% % % cfg.parameters.estimation.parameters 
% % % 
% % % with the names, bounds, and values of mouse.
% % % 
% % % Next say you only want to estimate a subset of the parameters. This can
% % % be done by passing a cell array of the parameter names:
% %
% % cfg=select_set(cfg, 'mouse', {'Vp', 'CL'});
% %
% % % In this case we want to estimate Vp and CL
% % % 
% % % The following fields are also created:
% %
% % cfg.estimation.to_estimate  - index of parameters being estimated
% % cfg.estimation.set_name     - source set name used 
% % cfg.estimation.mi           - data structure to map parameter names to their
% %                               estimation subset index (see below)
% %
% % cfg.estimation.parameters   - parameter information about the subset of
% %                               parameters to be estimated
% % cfg.estimation.parameters.guess - vector of parameter values used to start
% %                                   the estimation
% %
% % % to pull out the initial guess and overwrite the default 
% % % clearance with a  value of 0.1;
% %
% % parameters_guess =  cfg.estimation.parameters.guess;
% % parameters_guess(cfg.estimation.mi.CL) = 0.1;
% %
% %
% % % If the parameter subset to be estimated contains variance parameters, the
% % % objective type will default to maximum likelihood (ml). If there are no
% % % variance parameters then it will default to weighted least squares (wls).
% 
% 
% %
% % First we make sure the specified set exists
% %
%   % finding the indices of the parameters to estimate
%   % and storing them in cfg.estimation.to_estimate:
% if(~isfield(cfg.parameters.sets, set_name))
%   disp(sprintf('Error: parameter set >%s< was not defined in cfg', set_name));
%   disp(sprintf('       using the ''default'' parameter set      '));
%   set_name = 'default';
% end
% 
% % overwriting the values with those for the selected set:
% cfg.parameters.values = getfield(getfield(cfg.parameters.sets, set_name), 'values');
% 
% 
% % first we gather the indices of all of the parameters to be estimated
% % both system and variance
% tmp_to_estimate_all  = [];
% 
% % Now we populate that array depending on what has been specified by the user.
% % Either use the cell array of parameter names that were specified, or default
% % to all of the parameters
% cfg.estimation.set_name    = set_name;
% if(exist('parameter_names', 'var'))
%   % parameter_names was defined, but we still check to see if 
%   % it was empty and default to all parameters if it was
%   if(length(parameter_names) > 0)
%     for parameter_name_idx = 1:length(parameter_names)
%       tmp_to_estimate_all = [ tmp_to_estimate_all
%       getfield(cfg.options.mi.parameters, parameter_names{parameter_name_idx})];
%     end
%   else
%     disp(sprintf('Warning: parameter_names specified but emtpy '))
%     disp(sprintf('         defaulting to all parameters        '))
%     tmp_to_estimate_all = [1:length(cfg.parameters.values)];
%   end
% else
%   tmp_to_estimate_all = [1:length(cfg.parameters.values)];
% end
% 
% % Now we need to split the system and variance
% % parameters.  These will hold the numeric indices 
% % of each of the parameter types:
% tmp_to_estimate_system   = [];
% tmp_to_estimate_variance = [];
% 
% 
% for ctr  = 1:length(tmp_to_estimate_all)
%   % if it's a variance parameter we add that parameter index
%   % to the tmp variance vector
%   if(strcmp(cfg.parameters.ptype(tmp_to_estimate_all(ctr)), 'variance'))
%     tmp_to_estimate_variance = [tmp_to_estimate_variance 
%                                 tmp_to_estimate_all(ctr)];
%   else
%     tmp_to_estimate_system   = [tmp_to_estimate_system   
%                                 tmp_to_estimate_all(ctr)];
%   end
% end
% 
% 
% %now we stick them together
% cfg.estimation.to_estimate = [tmp_to_estimate_system; tmp_to_estimate_variance];
% 
% if(length(tmp_to_estimate_variance) > 0)
%    cfg.estimation.parameters.system = length(tmp_to_estimate_system);
%    cfg.estimation.objective_type    = 'ml';
%  else
%    cfg.estimation.objective_type    = 'wls';
% end
% 
% % setting default values that can be overwritten 
% cfg.estimation.options = optimset('Display', 'iter', 'maxfuneval', 1500);
% cfg.estimation.effort  =  1; 
% 
% 
% % Now we need to push any variance parameters to the end of the 
% % parameters vector so the estimation routines will do well
% 
% % now we're creating the cfg.estimation.parameters data 
% % structure to hold the subset of parameters we want to estimate 
% cfg.estimation.parameters.names       = cfg.parameters.names       (cfg.estimation.to_estimate);
% cfg.estimation.parameters.guess       = cfg.parameters.values      (cfg.estimation.to_estimate);
% cfg.estimation.parameters.lower_bound = cfg.parameters.lower_bound (cfg.estimation.to_estimate);
% cfg.estimation.parameters.upper_bound = cfg.parameters.upper_bound (cfg.estimation.to_estimate);
% cfg.estimation.parameters.units       = cfg.parameters.units       (cfg.estimation.to_estimate);
% cfg.estimation.parameters.editable    = cfg.parameters.editable    (cfg.estimation.to_estimate);
% cfg.estimation.parameters.ptype       = cfg.parameters.ptype       (cfg.estimation.to_estimate);
% cfg.estimation.parameters.type        = cfg.parameters.type        (cfg.estimation.to_estimate);
% 
% % storing the map between parameter name and index 
% % in the estimation structure:
% for p_idx = 1:length(cfg.estimation.parameters.names)
%  eval(sprintf('cfg.estimation.mi.%s = p_idx;', cfg.estimation.parameters.names{p_idx}));
% end
% 
% 
% 
% % applying set specific covariates 
% 
% % checking to make sure there are covariates
% % if there are, we assign the set values to them
% if(isfield(cfg.options, 'inputs'))
% if(isfield(cfg.options.inputs, 'covariates'))
% for cov_idx = 1:length(cfg.options.inputs.covariate_names)
%   mycov_str  = cfg.options.inputs.covariate_names{cov_idx};
%   mycov = getfield(cfg.options.inputs.covariates, mycov_str);
%   if(isfield(mycov.parameter_sets, set_name))
%     cov_set_values = getfield(mycov.parameter_sets, set_name);
%   else
%     cov_set_values = getfield(mycov.parameter_sets, 'default');
%   end
% 
%   mycov.times.values   = cov_set_values.times;
%   mycov.values.values  = cov_set_values.values;
%   cfg.options.inputs.covariates = setfield(cfg.options.inputs.covariates, mycov_str, mycov);
%   
% end
% end
% end
% 
