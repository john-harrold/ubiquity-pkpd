function [simout_mapped] = run_simulation(parameters,cfg)
% function [simout_mapped] = run_simulation(parameters,cfg)
%
% This runs a simulation for a model created in the system.txt format
%
% % compile   the system to make sure the 
% % latest changes have been committed. 
% build_system
%
% See the generated file:
%   transient/auto_simulation_driver.m 
% for examples on how to control different aspects of the simulation. 
%

% Setting up some default simulation options:
simulation_options.model_name             = 'ode_simulation';
simulation_options.default_simopts.Solver = 'ode23s';
simulation_options.output_times           = linspace(0,10,1000)';

% others will be added and defaults will be overwritten here:
if(isfield(cfg.options, 'simulation_options'))
  SIMINT_fields = fieldnames(cfg.options.simulation_options);
  for SIMINT_idx = 1:length(SIMINT_fields)
     eval(sprintf('simulation_options.%s = cfg.options.simulation_options.%s;', SIMINT_fields{SIMINT_idx}, SIMINT_fields{SIMINT_idx}));
  end
end

% The for loop will create a variable for each parameter automatically in this workspace
for parameter_idx = 1:length(parameters)
  % this defines the parameter
  eval(sprintf('%s = parameters(parameter_idx);',cfg.parameters.names{parameter_idx}));
  % this stores that paramter in the SIMINT_META data structure
  % so that it can be passed along with the simulation
  % output for validation purposes
  eval(sprintf('SIMINT_META.parameters.%s = parameters(parameter_idx);',cfg.parameters.names{parameter_idx}));
end

% covariates are evaluated at the first value specified 
% so they can be used in the calculation of static secondary parameters
if(isfield(cfg.options, 'inputs'))
if(isfield(cfg.options.inputs, 'covariates'))
  SIMINT_start_time = min(simulation_options.output_times);
  for SIMINT_covariate_idx = 1:length(cfg.options.inputs.covariate_names)
   SIMINT_covariate_info = getfield(cfg.options.inputs.covariates, cfg.options.inputs.covariate_names{SIMINT_covariate_idx});
   % defining the covariate
   eval(sprintf('%s = SIMINT_covariate_info.values.values(1);',cfg.options.inputs.covariate_names{SIMINT_covariate_idx}));
   % storing the covariate value used for static secondary parameter
   eval(sprintf('SIMINT_META.covariate_IC.%s = SIMINT_covariate_info.values.values(1);',cfg.options.inputs.covariate_names{SIMINT_covariate_idx}));
  end
end
end


% Here are the static secondary parameters are created
% in an automated fashion:  
if(isfield(cfg.options.misc, 'static_secondary_parameters'))
  static_secondary_parameters = fieldnames(cfg.options.misc.static_secondary_parameters);
  for parameter_idx = 1:length(static_secondary_parameters)
    SIMINT_defined_flag = 0;

    % first we see if there are any static secondary parameters that are
    % defined by conditional statements:
    if(isfield(cfg.options.misc, 'static_secondary_parameters_conditional'))
      % Now we check this specific parameter:
      if(isfield(cfg.options.misc.static_secondary_parameters_conditional, static_secondary_parameters{parameter_idx}))
        % if the secondary parameter is defined in terms of a conditional
        % statement, then it's evaluated here
        eval(getfield(cfg.options.misc.static_secondary_parameters_conditional, static_secondary_parameters{parameter_idx}));

        % we flip the defined flag so we don't define it below using a regular
        % assignment statement
        SIMINT_defined_flag = 1;
      end
    end

    % If the defined flag is 0 then we didn't defien this parameter as 
    % a conditional and need to define it here
    if(SIMINT_defined_flag == 0)
      eval(sprintf('%s = %s;', ...
                    static_secondary_parameters{parameter_idx}, ...
           getfield(cfg.options.misc.static_secondary_parameters,static_secondary_parameters{parameter_idx})));
    end

    % this stores the secondary parameter in the SIMINT_META data structure so
    % that it can be passed along with the simulation output for validation
    % purposes
    eval(sprintf('SIMINT_META.secondary_parameters.%s = %s;', ...
             static_secondary_parameters{parameter_idx}, ...
             static_secondary_parameters{parameter_idx}))
  end 
end 

% First set the default initial condition to be zero for all states 
simulation_options.initialstate  = zeros(1,length(fieldnames(cfg.options.mi.states)));
 
% Next we overwrite those states that are nonzero 
% The following 'if' and 'for' chunks will define the initial conditions
% components of the simulation_options for the states that have non-zero
% initial conditions. This is done in an  automated fashion.  Similar 
% to the parameters, the explicit declarations are made below 
% (again for reference)
if(isfield(cfg.options, 'initial_conditions'))
  initial_conditions          = fieldnames(cfg.options.initial_conditions);
  for state_idx = 1:length(initial_conditions)
    eval(sprintf('simulation_options.initialstate(cfg.options.mi.states.%s) = %s;', ...
             char(initial_conditions(state_idx)), ...
         getfield(cfg.options.initial_conditions,char(initial_conditions(state_idx)))));
  end 
end 



% we look and see if there are any inputs 
% and if there are we look for each type
% and add them to the simulation_options 
% data structure
if(isfield(cfg.options, 'inputs'))
  %
  % Bolus inputs
  %
  if(isfield(cfg.options.inputs, 'bolus'))
    % first we setup the first row which is the bolus times
    eval(sprintf('SIMINT_bolus_times = cfg.options.inputs.bolus.times.values.*%s;', cfg.options.inputs.bolus.times.scale));
    simulation_options.bolus_inputs(1,:) = SIMINT_bolus_times;

    SIMINT_fields = fieldnames(cfg.options.inputs.bolus.species);
    for SIMINT_idx = 1:length(SIMINT_fields)
      SIMINT_species = SIMINT_fields{SIMINT_idx};
      
      eval(sprintf('SIMINT_bolus_values = cfg.options.inputs.bolus.species.%s.values;', SIMINT_species));
      % JMH moved scaling to a different function remove lines here and just
      % below later (2015.10)
      % eval(sprintf('SIMINT_bolus_scale  = cfg.options.inputs.bolus.species.%s.scale;', SIMINT_species)); 
      eval(sprintf('SIMINT_state_idx    = cfg.options.mi.states.%s;', SIMINT_species));

      
      % this creates the bolus magnitudes:
      % JMH (2015.10)
      % eval(sprintf('simulation_options.bolus_inputs(end+1,:) = SIMINT_bolus_values.*%s ;', SIMINT_bolus_scale));
      eval(sprintf('simulation_options.bolus_inputs(end+1,:) = SIMINT_bolus_values ;'));
      % and here we specify the state id associated with that bolus value:
                    simulation_options.bolus_inputs(end+1,:) = SIMINT_state_idx*ones(size(SIMINT_bolus_times));
    end
  end
  %
  % Infusion rates
  %
  if(isfield(cfg.options.inputs, 'infusion_rate_names'))
    SIMINT_fields = fieldnames(cfg.options.inputs.infusion_rates);
    SIMINT_fields = cfg.options.inputs.infusion_rate_names;
    for SIMINT_idx = 1:length(SIMINT_fields)
      SIMINT_rate    = SIMINT_fields{SIMINT_idx};
      eval(sprintf('SIMINT_rate_details = cfg.options.inputs.infusion_rates.%s;', SIMINT_rate));

      eval(sprintf('simulation_options.infusion_rates(%d).name         = ''%s'';',SIMINT_idx, SIMINT_rate) );
      eval(sprintf('simulation_options.infusion_rates(%d).levels(1,:)  = SIMINT_rate_details.times.values.*%s;', SIMINT_idx, SIMINT_rate_details.times.scale));
      eval(sprintf('simulation_options.infusion_rates(%d).levels(2,:)  = SIMINT_rate_details.levels.values.*%s;',SIMINT_idx, SIMINT_rate_details.levels.scale));
    end
  end

  %         
  % Covariates         
  %
  if(isfield(cfg.options.inputs, 'covariates'))
     simulation_options = process_covariates(cfg, simulation_options);
  end

end

% running the simulation
tstart   = cputime;
[simout]=run_simulation_generic(parameters, simulation_options, cfg);
tstop    = cputime;
SIMINT_META.timing.simulation = tstop-tstart;

% mapping the outputs and states
simout_mapped = auto_map_simulation_output(simout, cfg);

% adding error to the output
simout_mapped = add_observation_errors(simout_mapped, parameters, cfg);


simout_mapped.meta = SIMINT_META;


function [som]=add_observation_errors(som, parameters, cfg)

% pulling out the error models
em = fieldnames(cfg.ve);

% now we loop through them and add them one at a time
for(emidx = 1:length(em))
   current_em = getfield(cfg.ve, em{emidx});

   som =  ...
   output_add_error(som,                         ... % simulation output without error
                    em{emidx},                   ... % output
                    getfield(cfg.ve, em{emidx}), ... % error model
                    parameters,                  ... % current parameter values
                    cfg);

end


function [SIMINT_som]=output_add_error(SIMINT_som, SIMINT_output, SIMINT_em, SIMINT_parameters, SIMINT_cfg)
% function [SIMINT_som]=output_add_error(SIMINT_som, SIMINT_output, SIMINT_em, SIMINT_parameters, SIMINT_cfg)
% This function takes the current simulation output and adds the error
% specified by SIMINT_em to the utput SIMINT_output for the given set of
% parameters SIMINT_parameters
%
%
% SIMINT_som                  mapped simulation output without error
% SIMINT_output               output name
% SIMINT_em                   error model
% SIMINT_parameters           vector of parameters
% SIMINT_cfg                  cfg data structure describing the system

% Defining the time
SIMINT_TIME = SIMINT_som.times.sim_time;


% defining the PRED values locally 
eval(sprintf('PRED  = SIMINT_som.outputs.%s;', SIMINT_output));


% pulling out the variance parameter names
SIMINT_VP_NAMES = SIMINT_cfg.parameters.names(strcmp(SIMINT_cfg.parameters.ptype, 'variance'));

% vectorize the equation in em
SIMINT_em = strrep(SIMINT_em, '^', '.^');
SIMINT_em = strrep(SIMINT_em, '/', './');
SIMINT_em = strrep(SIMINT_em, '*', '.*');


for(SIMINT_IDX = 1:length(SIMINT_VP_NAMES))
    % defining the variance parameters locally
    eval(sprintf('%s = system_fetch_parameter(SIMINT_cfg, SIMINT_parameters, ''%s'');', ...
              SIMINT_VP_NAMES{SIMINT_IDX}, ...
              SIMINT_VP_NAMES{SIMINT_IDX}))

end

% calculating the error model value
SIMINT_em_VARIANCE  =  eval(SIMINT_em);

% normrnd takes as the second argument the standard deviation
% which is the sqrt of the variance
eval(sprintf('SIMINT_ERR = normrnd(zeros(size(PRED)), sqrt(SIMINT_em_VARIANCE));', ...
          SIMINT_VP_NAMES{SIMINT_IDX}, ...
          SIMINT_VP_NAMES{SIMINT_IDX}))

% Adding the error to the prediction and storing it in som
eval(sprintf('SIMINT_som.outputs.SIOE_%s = PRED + SIMINT_ERR;', SIMINT_output));


