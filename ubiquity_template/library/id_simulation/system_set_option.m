function [cfg]=system_set_option(cfg, group, option, value)
%
% function [cfg]=system_set_option(cfg, group, option, value)
%
%
%
%
%=======================================================================================
%
% Group: 'estimation'
%
%   The estimation options are initialized when a parameter set has been
%   selected. With a list of parameters to be estiatimated. For example if we
%   only want to estimate Vp and kel and we are using the default parameter
%   set, we would start by doing the following:
%
%   to_estimate = {'Vp'
%                  'kel'}';
%   
%   % get the system information:
%   cfg = auto_fetch_system_information;
%   
%   % define the parameter set and the parameters
%   % we want to include in the estimation
%   cfg = select_set(cfg, 'default', to_estimate);
%   
%   % if you want to estimate all  the parameters then 
%   % select the parameter set in the following way:
%   cfg = select_set(cfg, 'default');
%
%   The estimation can then be controlled using the following options:
%
%   Option: 'objective_type'
%   ------------------------
%   Values: 'ml', 'wls'
%
%   If there are variance parameters being estimated then the objective type
%   will be set to maximum likelihood ('ml'). Otherwise it will default to weighted
%   least squares ('wls'). 
%
%   Default:
%   cfg = system_set_option(cfg, 'estimation', 'objective_type', 'wls');
%
%   Option: 'observation_function'
%   ------------------------------
%   Values: string of m-file name without the .m extension
%
%   The function observation_details.m  must take two inputs: parameters and
%   cfg (in that order) --- parameters will be the subset of parameters being
%   estimated. The first value returned is a data structure referred to here as
%   od. This data structure should have a field for each group: group1, group2,
%   ....  Each group/output combination should have a matrix containing
%   observation information:
%
%     observations.group1.output1 = [observation_matrix];
%     observations.group1.output2 = [observation_matrix];
%     observations.group2.output2 = [observation_matrix];
%  
%   The format of the observation_matrix is a numerical array and
%   should have the following format:
%         
%                  --                          --    
%                  | t(1)  y(1)   yhat(1) var(1)|
%                  | t(2)  y(2)   yhat(2) var(2)|
%                  |  .     .        .      .   |
%                  |  .     .        .      .   |
%                  |  .     .        .      .   |
%                  |  .     .        .      .   |
%                  |  .     .        .      .   |
%                  | t(m)  y(m)   yhat(m) var(m)|
%                  --                          --
%  
%   Each output for each subject can have m different measurements (y)
%   at time (t). The model predictions are given by yhat, and var is
%   the variance of y. 
%
%
%   Calculating var:
%   
%   For WLS:
%   If you just want to do sum-squared error
%   
%     var = 1;
%   
%   The inverse of the prediction squared would be defined in the
%   following manner:
%  
%      var(i) = (yhat(i))^2;
%  
%   For ML:
%   Modeling the variance, where the variance parameters are stored in the 
%   parameter vector:
%   
%      slope      =  parameters(cfg.options.mi.parameters.slope);
%      intercept  =  parameters(cfg.options.mi.parameters.intercept);
%  
%      var(i) = (intercept + slope*yhat(i))^2;
%  
%   NOTE: while the first output must contain the observation details
%   described above, more than one output can be returned if needed. The format
%   of additional outputs can be defined by the user.
%
%   Default:
%   cfg = system_set_option(cfg, 'estimation', 'observation_function', 'system_od_general');
%
%   Option: 'effort'
%   ----------------
%   Values: integer >= 1
%
%   The 'effort' field tells the estimation routine to try "harder" to find a
%   good parameter estimate. By default it will be assumed to be 1, and will just 
%   estimate the parameters like normal. If you set it to a value greater than
%   1, that many iterations will be used to move around the parameter space
%   and try to find a better solution.
%
%   NOTE: this is only valid for the fminsearch optimizer.
%
%   cfg = system_set_option(cfg, 'estimation', 'effort', 1);
%
%   Option: 'optimizer'            
%   -------------------            
%   Values: 'fminsearch', 'ga'
%
%   Matlab provides different optimization functions that can be used to
%   perform parameter estimation. The 'fminsearch' function is a general
%   gradient based optimization routine. If you are having difficulty finding
%   a good solution because of local minima the genetic algorithm funciton
%   ('ga') can be used if the Global Optimization Toolbox is installed.
%
%   Default:
%   cfg = system_set_option(cfg, 'estimation', 'optimizer', 'fminsearch');
%
%   Option: 'optimization_options'
%   ----------------------------------------
%   Values: optimset data structure
%
%   The optimization routines have different options. For 'fminsearch' use
%   'optimset', for 'ga' use 'gaoptimset'. By default the options for
%   'fminsearch' will be used. You can overwrite these in the following
%   manner:
%
%   cfg = system_set_option(cfg, 'estimation', 'optimization_options', ...
%                optimset('Display',   'iter', ...
%                         'TolFun',     1e-3,  ...
%                         'MaxIter',    3000,  ...
%                         'MaxFunEval', 3000))
%
%   Those are a few examples. See the documentation for optimset for a
%   detailed list of options that can be passed down.
%
%   Default:
%   cfg = system_set_option(cfg, 'estimation', 'optimization_options', optimset());
%
%   Option: 'monitor_status_function'
%   ---------------------------------
%   Values: string of m-file name without the .m extension
%
%   To view the progress of the estimation graphically, you need to specify an
%   output function. The function library/id_simulation/estimation_status.m
%   can be used to do this. Just specify the following option 'estimation_status'.
%
%   To make your own custom function copy estimation_status.m to the main
%   template directory and rename it something, say mystatus.m. Then you can 
%   use the value 'mystatus' to specify that the new function should be used.
%
%   NOTE: This will only be used when effort is set to 1, and is only valid for
%   the fminsearch optimizer.
%
%   Default:
%   cfg = system_set_option(cfg, 'estimation', 'monitor_status_function', 'estimation_status');
%
%   To disable set the monitor status funciton to ''
%   cfg = system_set_option(cfg, 'estimation', 'monitor_status_function', '');
%
%   Option: 'monitor_exit_when_stable'
%   ----------------------------------
%   Values: 'yes', 'no'
%
%   Description
%   With monitoring enabled the estimation_status function allows the
%   estimation to terminate when a "stable" solution has been found. Set to
%   'yes' to enable this option. See monitor_iteration_history and
%   monitor_slope_tolerance below to control this.
%
%   Default:
%   cfg = system_set_option(cfg, 'estimation', 'monitor_exit_when_stable', 'no');
%
%   Options: 'monitor_iteration_history' & 'monitor_slope_tolerance'
%   ----------------------------------------------------------------
%   Values: integer > 1 & positive number > 0
%
%   This will tell the optimizer to finish up when the parameters and
%   objective function have 'stabilized'. This is accomplished by fitting a
%   line to these values over a specified iteration history and comparing the
%   largest slope to a tolerance. 
%
%   Default:
%   cfg = system_set_option(cfg, 'estimation', 'monitor_iteration_history', 100);
%   cfg = system_set_option(cfg, 'estimation', 'monitor_slope_tolerance',   0.001);
%
%=======================================================================================
%
% Group: 'logging'   
%
%   Option: 'enabled'
%   -----------------
%   Values: 'yes', 'no'
%
%   Default:
%   cfg = system_set_option(cfg, 'logging', 'enabled', 'yes');
%
%   Description: 
%      - no   -> Disabled logging
%
%      - yes  -> Enable Logging   
%    
%  
%   Option: 'file'     
%   -----------------
%   Values: string with log file name
%
%   Default:
%   cfg = system_set_option(cfg, 'logging', 'file', sprintf('transient%subquity_log.txt', filesep)); 
%
% 
%   Option: 'timestamp'
%   -----------------
%   Values: 'yes', 'no'
%
%   Default:
%   cfg = system_set_option(cfg, 'logging', 'timestamp',  'yes');                                          
%
%   Description: 
%      - no   -> Just log the string provided
%
%      - yes  -> Prepend a timestamp to the string provided
%
%   Option: 'ts_str'
%   -----------------
%   Values: string describing the format of the time stamp
%
%   For more information on this format, see the help for the datestr function
%   in Matlab.
%
%   Default:
%   cfg = system_set_option(cfg, 'logging', 'ts_str',     'yyyy-mm-dd HH:MM:SS)';                       
%
%  
%=======================================================================================
%
% Group: 'simulation'
%
%   Option: 'integrate_with'
%   -----------------------
%   Values: 'm-file', 'simulink'
%
%
%   Default:
%   cfg = system_set_option(cfg, 'simulation', 'integrate_with', 'm-file');
%
%   Description: 
%      - m-file   -> Uses matlabs m-file scripting interface for running 
%                    simulations (default method)
%    
%      - simulink -> uses compiled C and is faster, but you may not have a
%                    license for it
%
%   Option: 'include_important_output_times'
%   ----------------------------------------
%   Values: 'yes', 'no'
%
%   When simulating the system for generating plots, short infusions can be
%   skipped over during simulation, etc.  it may be important to force
%   sampling at important output times (such as bolus events, times when the
%   infusion rate changes, etc). Set this value to 'yes' to accomplish this
%   (default is 'yes'). To only sample at specified times change this value to
%   'no'
%
%   Default:
%   cfg = system_set_option(cfg, 'simulation', 'include_important_output_times', 'yes');
%
%
%   Option: 'output_times'
%   ----------------------
%   Values: Vector of increasing time values in the units of the simulation
%   time
%
%   This is vector of times where states and model outputs are evaluated. 
%
%   Default:
%   cfg = system_set_option(cfg, 'simulation', 'output_times', linspace(0,100));
%
%=======================================================================================
%
% Group: 'solver'
%
%   Option: 'solver'
%   ----------------
%   Values: 'ode15s', 'ode23s', 'ode45', 'ode23', 'ode113', 'ode23t', 'ode23tb', 'ode15i'
%
%   This is where the ode solver is specified. For nonstiff systems ode45
%   works well. For stiff sets of equations ode15s and ode23s generally
%   perform well.
%
%   Default:
%   cfg = system_set_option(cfg, 'solver', 'solver', 'ode23s');
%
%
%   Option: other options 
%   ---------------------
%   Any name/value paring from odeset can be used. See 'help odeset' for a
%   list of name/value parings. For example to set the relative and absolute
%   error tolerances to 1e-5 and 1e-9 respectively, the following would be
%   used:
%
%   cfg = system_set_option(cfg, 'solver', 'RelTol',  1e-5);
%   cfg = system_set_option(cfg, 'solver', 'AbsTol',  1e-9);
%
%=======================================================================================
%
% Group: 'stochastic'
%
%   Option: 'ci'
%   ----------------------------------------
%   Values: integer >= 1
%
%   This is the desired confidence (or prediction) interval in percentage that
%   we want to use when analyzing the responses for each state and model output.
%
%   Default:
%   cfg = system_set_option(cfg, 'stochastic', 'ci', 95);
%
%   Option: 'nsub'
%   ----------------------------------------
%   Values: integer >= 1 
%
%   This specifies the number of subjects to generate.
%
%   Default:
%   cfg = system_set_option(cfg, 'stochastic', 'nsub', 100);
%
%   Option: 'seed'
%   ----------------------------------------
%   Values: integer >= 1 
%
%   This is the seed the random number generator will use.
%
%   Default:
%   cfg = system_set_option(cfg, 'stochastic', 'seed', 8675309);
%
%   Option: 'ponly'
%   ----------------------------------------
%   Values:  true or false
%
%   By default simulating subjects will simulate the parameter sets for each
%   subject and their response to the specified dosing. If you simply want the
%   parameter values set 'ponly' to true.
%
%   Default:
%   cfg = system_set_option(cfg, 'stochastic', 'ponly', false);
%
%   Option: 'states'
%   ----------------------------------------
%   Values:  cell array of state names
%
%   By default all states will be returned but because stochastic simulations
%   can take up a lot of memory, it may be useful to specify a subset of
%   states
%
%   Default:
%   cfg = system_set_option(cfg, 'stochastic', 'states', fieldnames(cfg.options.mi.states));
%
%   Option: 'outputs'
%   ----------------------------------------
%   Values:  cell array of state names
%
%   By default all outputs will be returned but because stochastic simulations
%   can take up a lot of memory, it may be useful to specify a subset of
%   outputs
%
%   Default:
%   cfg = system_set_option(cfg, 'stochastic', 'outputs', fieldnames(cfg.options.mi.outputs));




% list of good group specifications
groups = {'estimation', 'logging', 'solver', 'stochastic', 'simulation'};

isgood = true;
error_msg = '';


if(sum(strcmp(groups, group)))
  % setting stochastic options
  if(strcmp(group, 'stochastic'))

    % if we're setting states and outputs we check to make
    % sure they exist
    if(strcmp(option, 'states') | strcmp(option, 'outputs'))
      valid = fieldnames(getfield(cfg.options.mi, option));
      for(validx = 1:length(value))
        if(sum(strcmp(valid, value{validx})) == 0)
          isgood = false;
          error_msg = sprintf('%s %s', error_msg, value{validx});
        end
      end
    end


    if(isgood)
      eval(sprintf('cfg.options.stochastic.%s = value;', option));
    else
      print_banner_start(cfg);
      vp(cfg, sprintf('The following %s are not valid', option));
      vp(cfg, error_msg);
      vp(cfg, sprintf(' %s were not set', option));
      print_banner_stop(cfg);
    end
  end

  if(strcmp(group, 'logging'))
    eval(sprintf('cfg.options.logging.%s = value;', option));
  end
  
  % setting simulation options
  if(strcmp(group, 'simulation'))
    eval(sprintf('cfg.options.simulation_options.%s = value;', option));
  end
  
  % setting solver options
  if(strcmp(group, 'solver'))
    eval(sprintf('cfg.options.simulation_options.solver_opts.%s = value;', option));
  end

  % setting estimation options
  if(strcmp(group, 'estimation'))

    % allowable options
    options.objective_type.cfgname               = 'cfg.estimation.objective_type';
    options.observation_function.cfgname         = 'cfg.estimation.observation_function';
    options.optimization_options.cfgname         = 'cfg.estimation.options';
    options.effort.cfgname                       = 'cfg.estimation.effort' ;  
    options.monitor_status_function.cfgname      = 'cfg.estimation.monitor.status_function';  
    options.monitor_exit_when_stable.cfgname     = 'cfg.estimation.monitor.exit_when_stable';  
    options.monitor_iteration_history.cfgname    = 'cfg.estimation.monitor.iteration_history';  
    options.monitor_slope_tolerance.cfgname      = 'cfg.estimation.monitor.slope_tolerance';  
    options.optimizer.cfgname                    = 'cfg.estimation.optimizer' ;  


    % format of the options
    options.objective_type.type                  = 'string';
    options.observation_function.type            = 'free';
    options.optimization_options.type            = 'free';
    options.effort.type                          = 'integer';
    options.monitor_status_function.type         = 'free';
    options.monitor_exit_when_stable.type        = 'string';
    options.monitor_iteration_history.type       = 'integer';
    options.monitor_slope_tolerance.type         = 'float';
    options.optimizer.type                       = 'string';

    % allowed values
    options.objective_type.allowed               = {'ml', 'wls'};        
    options.observation_function.allowed         = {};        
    options.optimization_options.allowed         = {};        
    options.effort.allowed                       = {};        
    options.monitor_status_function.allowed      = {};        
    options.monitor_exit_when_stable.allowed     = {'yes', 'no'};        
    options.monitor_iteration_history.allowed    = {};        
    options.monitor_slope_tolerance.allowed      = {};        
    options.optimizer.allowed                    = {'fminsearch', 'ga'};
                                                      
     
  if(sum(strcmp(fieldnames(options), option)))
    option_type    = getfield(getfield(options, option), 'type');
    option_allowed = getfield(getfield(options, option), 'allowed');
    option_cfgname = getfield(getfield(options, option), 'cfgname');

    option_ok = true;
    if(strcmp(option_type, 'string'))
      if(not( sum(strcmp(option_allowed, value))))
        option_ok = false;
        error_msg = sprintf('group (%s), option (%s), should be a character string (%s)', group, option,  strjoin(option_allowed, ', '))
      end
    elseif(strcmp(option_type, 'integer'))
      if(not(value == floor(value)))
      value
        option_ok = false;
        error_msg = sprintf('group (%s), option (%s), should be an integer', group, option)
      end
    elseif(strcmp(option_type, 'float'))
      if(not(isfloat(value)))
        option_ok = false;
        error_msg = sprintf('group (%s), option (%s), should be an number', group, option)
      end
    end

    % Making sure the appropriate toolbox is available 
    if(strcmp(group,  'estimation') & ...)
       strcmp(option, 'optimizer')  & ...
       strcmp(value, 'ga'))
      if(~license('checkout', 'gads_toolbox') | (exist('ga') ~= 2))
        vp(cfg, sprintf('Cannot find the Global Optimization Toolbox'));
        vp(cfg, sprintf('This is required to use the Genetic Algorithm'));
        vp(cfg, sprintf('for parameter estimation'));
        error_msg = sprintf('group (%s), option (%s) not set', group, option);
        option_ok = false;
      end

    end

    if(option_ok)
      % if everything checks out then we assign the option
      eval(sprintf('%s = value;', option_cfgname));
    else
      % if something happened above, then we print an error here:
      print_banner_start(cfg);
      vp(cfg, error_msg);
      print_banner_stop(cfg);
    end
  else
    print_banner_start(cfg);
    vp(cfg, sprintf('For the estimation option group the '));
    vp(cfg, sprintf('specified option (%s) is invalid', option));
    print_banner_stop(cfg);
  end


  end

else

  print_banner_start(cfg);
  vp(cfg, sprintf('The specified group (%s) is invalid', group))
  vp(cfg, sprintf('Valid groups are: '))

  for(gidx = 1:length(groups))
    vp(cfg, sprintf(' -> %s', groups{gidx}));
  end
  print_banner_stop(cfg);
end


function []=print_banner_start(cfg)
  vp(cfg, sprintf('------------------------------------')) 
  vp(cfg, sprintf('system_set_option()                 ')) 

function []=print_banner_stop(cfg)
  vp(cfg, sprintf('------------------------------------')) 

