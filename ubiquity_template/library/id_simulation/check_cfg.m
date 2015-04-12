function [status]=check_cfg(cfg)
% function [status]=check_cfg(cfg)
%
%  Checks the 'cfg' configuration to make sure the correct fields are defined
%
%
%  Most of these fields are created automatically for you. After building a
%  system from a system.txt file, the base cfg data structure can be created 
%  using:
%
%    cfg = auto_fetch_system_information;
%
%  Next the estimation information can be specified using select_set:
%
%    cfg = select_set(cfg, 'default');
%
% inputs:
%    cfg    configuration variable used in simulation and parameter estimation
%
%    Required fields:
%       cfg.estimation.observation_function
%                     .objective_type
%
%       cfg.estimation.parameters.names  
%                                .guess
%                                .lower_bound
%                                .upper_bound
%                                .units       
%                                .editable    
%                                .variance    
%                                .system    (ML only)
%
%      see  help select_set for details on how to create this
%
%    Optional fields:
%       cfg.estimation.timeout
%
%    Options
%        cfg.options.verbose 
%          Set this to a value of 'yes' to see progress and other information 
%          displayed 
%
%
%    User defined optional fields:
%        These are fields that you can define and use internally in the
%        observation_function or within your analysis:
%        cfg.data
%        cfg.options
%
%
%    Detailed description:
%
%       cfg.estimation.observation_function
%           Name of the m-file which takes as its first input the parameter
%           vector (system and variance) and it's second input is cfg. The
%           m-file should return a data structure of observations. Each filed
%           in observations represents a subject and each subject can have a
%           field for observations corresponding to each output. For example
%           if we have one subject with measurements for outputs 1 and 2, and
%           if we have another subject with measuremnets for output 2 the
%           observations data stucture will look like:
%
%             observations.subject1.output1 = [observation_matrix];
%             observations.subject1.output2 = [observation_matrix];
%             observations.subject2.output2 = [observation_matrix];
%
%           The format of the observation_matrix is a numerical array and
%           should have the following format:
%                 
%                          --                          --    
%                          | t(1)  y(1)   yhat(1) var(1)|
%                          | t(2)  y(2)   yhat(2) var(2)|
%                          |  .     .        .      .   |
%                          |  .     .        .      .   |
%                          |  .     .        .      .   |
%                          |  .     .        .      .   |
%                          |  .     .        .      .   |
%                          | t(m)  y(m)   yhat(m) var(m)|
%                          --                          --
%
%           Each output for each subject can have m different measurements (y)
%           at time (t). The model predictions are given by yhat, and var is
%           the variance of y. 
%           
%           Calcualting var:
%           
%           For WLS:
%           If you just want to do sum-squared error
%           
%             var = 1;
%           
%           The inverse of the prediction squared would be defined in the
%           following manner:
%
%              var(i) = (yhat(i))^2;
%
%           For ML:
%           Modeling the variance, where the variance parameters are stored in
%           elements 8 and 9 of your parameter vector:
%           
%              slope      = parameters(8);
%              intercept  = parameters(9);
%
%              var(i) = (intercept + slope*yhat(i))^2;
%
%
%       cfg.estimation.objective_type
%           This field should be a string indiciating the type of objective
%           function to use. The possible options are
%               
%           cfg.estimation.objective_type = 'ml';  % maximum likelihood
%           cfg.estimation.objective_type = 'wls'; % weighted least squares
%
%
%       cfg.estimation.timeout
%           This is a numerical value representing the time (in seconds)
%           after which the estimation will terminate. This will abandon the
%           estimation after 5 mintutes:
%
%           cfg.estimation.timeout = 600; 
%
%           Note: in order for this to work, you must first issue the 'tic'
%           command before calling an optimization routine. For example, if
%           you were using fminsearch to find parameter estimates, your call
%           might look something like:
%
%           tic;
%           [pest] = fminsearch(....);
%
%       cfg.estimation.parameters.guess
%           Vector containting the initial guess for parameters (system and
%           variance). For example if we had a once compartment system with
%           initial estimates of kel=.4 and Vc=10 witha proportional variance
%           model (sigma =.1) .
%       
%           cfg.estimation.parameters.guess  = [.4, 10, .1];
%       
%       cfg.estimation.parameters.names  
%           Cell array of strings containing the parameter names. The
%           dimensionality should be the same as cfg.estimation.parameters.guess. For
%           example:
%
%           cfg.estimation.parameters.names  = [{'kel'}, {'Vc'}, {'sigma'}];
%
%       cfg.estimation.parameters.lower_bound
%           Lower bound of parameter estimates. This elements in this vector
%           correspond to those in the 'guess' above, so the length is the
%           same. To specifiy positive bounds, we can set the lower bound to
%           the machine precition 'eps':
%
%           cfg.estimation.parameters.lower_bound = [eps eps eps];
%
%       cfg.estimation.parameters.upper_bound
%           Like the lower bound, this just represents and upper bound on the
%           parameter estimates. To leave them unbounded, these values can
%           just be set to 'inf' (infinity):
%
%           cfg.estimation.parameters.upper_bound = [inf inf inf];
%
%       cfg.estimation.parameters.system
%           This field is only relevant when using the 'ml' objective type.
%           The parameters vector contains both the system and variance
%           parameters. The first 'p' elements representing the system
%           parameters and the last 'q' elements representing those of the
%           variance. This field in the cfg.estimation.parameters structure
%           indicates the number of system parameters (p). Continuing with the
%           one compartment model with a proportional error model we have:
%
%           cfg.estimation.parameters.system = 2;
%
%
% outputs:
%    status 
%       Character string with a value of 'good' if everything is dandy and
%       'bad' if an error has been encountered. More verbose information
%       should be printed to the screen.
%

status     = 'good';


if(isfield(cfg, 'estimation'))
    if(~isfield(cfg.estimation,'observation_function'))
      status  = 'bad' ;
      disp('Missing Field: cfg.estimation.observation_function is empty');
    end 
    if(isfield(cfg.estimation,'objective_type'))
       if(~(strcmp(cfg.estimation.objective_type, 'ml') | strcmp(cfg.estimation.objective_type, 'wls')))
         status  = 'bad' ;
         disp('Format: cfg.estimation.objective_type can be either ''ml'' or ''wls'' is empty');
       end
    else
      status  = 'bad' ;
      disp('Missing Field: cfg.estimation.objective_type is empty');
    end 
else
  status  = 'bad' ;
  disp(' --> no information for controlling estimation has been provided');
  disp('Missing Field: cfg.estimation.observation_function');
  disp('Missing Field: cfg.estimation.objective_type');
end

if(isfield(cfg,'parameters'))
    num_params = 0;
    %
    % making sure 'guess' has been specified and that the format
    % is correct
    %
    if(isfield(cfg.estimation.parameters,'guess'))
      num_params = length(cfg.estimation.parameters.guess);
      if (num_params < 1)
        status  = 'bad' ;
        disp('Format: cfg.estimation.parameters.guess is empty');
      end
    else
      status  = 'bad' ;
      disp('Missing Field: cfg.estimation.parameters.guess');
    end

    %
    % checking for the presence and the dimensionality of 
    % the other parameter fields
    %
    parameter_fields = [{'names'}, {'lower_bound'}, {'upper_bound'}];
    
    for i=1:length(parameter_fields)
      if(isfield(cfg.estimation.parameters,char(parameter_fields(i))))
        if(num_params > 0)
          if( num_params ~= length(getfield(cfg.estimation.parameters,char(parameter_fields(i)))));
            status  = 'bad' ;
            disp(sprintf('Format: cfg.estimation.parameters.guess and cfg.estimation.parameters.%s have different lengths', char(parameter_fields(i))));
          end
        end
      else
        status  = 'bad' ;
        disp(sprintf('Missing Field: cfg.estimation.parameters.%s',char(parameter_fields(i))  ));
      end
    end




else
  status  = 'bad' ;
  disp(' --> no information about the parameters has been provided');
  disp('Missing Field: cfg.estimation.parameters.names');
  disp('Missing Field: cfg.estimation.parameters.guess');
  disp('Missing Field: cfg.estimation.parameters.lower_bound');
  disp('Missing Field: cfg.estimation.parameters.upper_bound');
end

