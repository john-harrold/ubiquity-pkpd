function [parameters_est, statistics_est]=estimate_parameters(cfg)
% function [parameters_est, statistics_est]=estimate_parameters(cfg)
% 
% % the following files are generated and placed in the
% % output directory 
%
%   parameters_all.csv % csv file with all parameters (estimated and fixed)
%   parameters_est.csv % csv file with only the estimated parameters 
%   report.txt         % summary report with estimates, cis and
%                      % variance/covariance matrix
% 
%
%
% For options on how to control parameter estimation see
%
% help system_set_option


% % assuming you've created the system.txt file to 
% % describe the system, this makes sure all the 
% % files are up to date:
% build_system;
% 
% % here we define the parameters to estimate:
% to_estimate = {'vp'
%                'kel'}';
%
% % get the system information:
% cfg = auto_fetch_system_information;
%
% % define the parameter set and the parameters
% % we want to include in the estimation
% cfg = select_set(cfg, 'default', to_estimate);
%
% % if you want to estimate all  the parameters then 
% % select the parameter set in the following way:
% cfg = select_set(cfg, 'default')
% 
% % If there are variance parameters being estimated then the objective type
% % will be set to maximum likelihood. Otherwise it will default to weighted
% % least squares. If you wish to overwrite them the following options can be
% % set:
%
% cfg.estimation.objective_type = 'ml';  % maximum likelihood
% cfg.estimation.objective_type = 'wls'; % weighted least squares
%
% 
% % You have to create the observation_details.m file: 
% cfg.estimation.observation_function = 'observation_details';
% 
% %  The function observation_details.m  must take two inputs: parameters and
% %  cfg (in that order) --- parameters will be the subset of parameters being
% %  estimated. The first value returned is a data structure referred to here as
% %  od. This data structure should have a field for each group: group1, group2,
% %  ....  Each group/output combination should have a matrix containing
% %  observation information:
% %
% %   observations.group1.output1 = [observation_matrix];
% %   observations.group1.output2 = [observation_matrix];
% %   observations.group2.output2 = [observation_matrix];
% %
% % The format of the observation_matrix is a numerical array and
% % should have the following format:
% %       
% %                --                          --    
% %                | t(1)  y(1)   yhat(1) var(1)|
% %                | t(2)  y(2)   yhat(2) var(2)|
% %                |  .     .        .      .   |
% %                |  .     .        .      .   |
% %                |  .     .        .      .   |
% %                |  .     .        .      .   |
% %                |  .     .        .      .   |
% %                | t(m)  y(m)   yhat(m) var(m)|
% %                --                          --
% %
% % Each output for each subject can have m different measurements (y)
% % at time (t). The model predictions are given by yhat, and var is
% % the variance of y. 
% % 
% % Calculating var:
% % 
% % For WLS:
% % If you just want to do sum-squared error
% % 
% %   var = 1;
% % 
% % The inverse of the prediction squared would be defined in the
% % following manner:
% %
% %    var(i) = (yhat(i))^2;
% %
% % For ML:
% % Modeling the variance, where the variance parameters are stored in
% % the parameter vector:
% % 
% %    slope      =  parameters(cfg.options.mi.parameters.slope)
% %    intercept  =  parameters(cfg.options.mi.parameters.intercept)
% %
% %    var(i) = (intercept + slope*yhat(i))^2;
% % NOTE: while the first output must contain the observation details
% % described above, more than one output can be returned if needed. The format
% % of additional outputs can be defined by the user.
% %
%
% % This will display the final estimates and
% % text to update the system.txt file
% % set it to 'no' if you don't want to see anything
% cfg.options.verbose                 = 'yes';
%
% % The 'effort' field tells the estimation routine to try "harder" to find a
% % good parameter estimate. By default it will be assumed to be 1, and will just 
% % estimate the parameters like normal. 
% cfg.estimation.effort = 1;
%
% % Increasing this number will implement a sort of simulated annealing where 
% % moving around the parameter space adding noise to the current estimates in 
% % an attempt to find a better solution. The larger this number, the longer
% % it will take.
%
% % The optimization routines have different options. By default the options
% % for fminsearch will be used. You can overwrite these in the folloiwng
% % manner:
% cfg.estimation.options = optimset('Display',   'iter', ...
%                                   'TolFun',     1e-3,  ...
%                                   'MaxIter',    3000,  ...
%                                   'MaxFunEval', 3000);
%
% % Those are a few examples. See the documentation for optimset for a
% % detailed list of options that can be passed down.
%
% % To view the progress of the estimation graphically, you need to specify an
% % output function. The function library/id_simulation/estimation_status.m
% % can be used to do this. Just specify the following option:
% cfg.estimation.monitor.status_function = 'estimation_status';
%
% % With monitoring enabled the estimation_status function allows the
% % following options:
% % Termination at stabilization 
% cfg.estimation.monitor.criteria.stability.exit_when_stable = 'yes';
%
% % This will tell the optimizer to finish up when the parameters and
% % objective function have 'stabilized'. This is accomplished by fitting a
% % line to these values over a specified iteration history and comparing the
% % largest slope to a tolerance. These parameters can be specified in the
% % following manner (defaults shown);
% cfg.estimation.monitor.criteria.stability.iteration_history = 100;
% cfg.estimation.monitor.criteria.stability.slope_tolerance   = 0.001;
%
% % To make your own custom function copy estimation_status.m to the main
% % template directory and rename it something, say mystatus.m. Then you can 
% % use the following to specify that the new function should be used.
% cfg.estimation.monitor.status_function = 'mystatus';
%
% % Then you can modify mystatus.m to make it do provide useful feedback at each
% % iteration. Note: When you make changes to this new function it may break
% % the stopping criteria specified above.

  
  if(not(isdir('output')))
    vp(cfg, '-----------------------------------');
    vp(cfg, 'Creating the "output" directory    ');
    vp(cfg, '-----------------------------------');
    mkdir('output');
  end

  % Setting things up for fminsearch
  if(strcmp(cfg.estimation.optimizer, 'fminsearch'))
    % Pulling the defaults for fminsearch
    options = optimset('fminsearch');
    % if the user has specified an estimation status function, then we
    % add that to the optimset options.
    if(isfield(cfg.estimation, 'monitor'))
      if(isfield(cfg.estimation.monitor, 'status_function'))
        if(~strcmp(cfg.estimation.monitor.status_function, ''))
          eval(sprintf('options = optimset(options, ''OutputFcn'', @(optimValues,state,varargin)%s(optimValues,state,cfg));', cfg.estimation.monitor.status_function));
        end
      end
    end
  end

  if(strcmp(cfg.estimation.optimizer, 'ga'))
    % Pulling the defaults for ga
    options = gaoptimset()
  end


  if( isfield(cfg.estimation, 'options'))
    %merging the user specified options with the default values
    options = optimset(options, cfg.estimation.options);
  end



  % making sure there is a default  value for effort
  if(~isfield(cfg.estimation, 'effort'))
    cfg.estimation.effort = 1;
  end

  % Defaulting to run the estimation then we run some 
  % checks to make sure everything is kosher. 
  run_estimation = 'yes';

  % First we look at the cfg variable to make sure it's all good:
  if(~strcmp(check_cfg(cfg), 'good'))
    disp('#-> Check your cfg data structure');
    run_estimation = 'no';
  end

  % Next we check the objective function to make sure it's giving something
  % reasonable back at the initial guess:
  if(~isfinite( calculate_objective(cfg.estimation.parameters.guess, cfg)))
    disp('#-> The objective function evaluates as either Inf or NaN');
    disp('#-> for the initial guess. Check your initial parameter guess');
    disp('#-> or the observation_details function');
    run_estimation = 'no';
  end

  % disabling logging in the simulation 
  cfg.options.simulation_options.logging = 'no';

  % Now we compare the objective type to the system 
  % information. We can then let the user know if
  % something doesn't make sense
  if(strcmp(cfg.estimation.objective_type, 'wls') && sum(strcmp(cfg.estimation.parameters.ptype, 'variance')))
    vp(cfg, '-----------------------------------');
    vp(cfg, 'Weighted Least Squares Objective   ');
    vp(cfg, 'selected _AND_ variance parameters ');
    vp(cfg, 'are being estimated. Consider      ');
    vp(cfg, 'switching to maximum likelihood:   ' );
    vp(cfg, 'cfg.estimation.objective_type = ''ml''');
    vp(cfg, '-----------------------------------');
  end

  if(strcmp(cfg.estimation.objective_type, 'ml') && (sum(strcmp(cfg.estimation.parameters.ptype, 'variance'))<1))
    vp(cfg, '-----------------------------------');
    vp(cfg, 'Maximum Likelihood  Objective   ');
    vp(cfg, 'selected _AND_ no variance parameters ');
    vp(cfg, 'are being estimated. Consider      ');
    vp(cfg, 'switching to weighted least squares' );
    vp(cfg, 'cfg.estimation.objective_type = ''wls''');
    vp(cfg, '-----------------------------------');
  end

  % Now we perform the estimation
  if(strcmp(run_estimation, 'yes'))

    % Here the optimizer field is used to determine the 
    % function used to estimate parameters
    vp(cfg,  '-----------------------------------');
    vp(cfg,  'Beginning parameter estiamtion');
    if(strcmp(cfg.estimation.objective_type, 'wls'))
      vp(cfg,  '   Obj: Weighted Least Squares');
    elseif(strcmp(cfg.estimation.objective_type, 'ml'))
      vp(cfg,  '   Obj: Maximum Likelihood');
    end

    % Using the genetic algorithm:
    if(strcmp(cfg.estimation.optimizer, 'ga'))

      vp(cfg,  '   Method: Genetic Algorithm (ga)');
      vp(cfg,  '-----------------------------------');
     
      initial_guess = cfg.estimation.parameters.guess ;

      if(iscolumn(initial_guess))
        initial_guess = initial_guess';
      end

      options = gaoptimset(options, 'InitialPopulation', initial_guess);
      fitnessfun = @(pvalues)calculate_objective(pvalues, cfg);  
      parameters_est = ...
      ga(fitnessfun, ...
         length(initial_guess), ...
         [], [], [], [],  ...
         cfg.estimation.parameters.lower_bound, ...
         cfg.estimation.parameters.upper_bound, ...
         [], [],  ...
         options);
    end


    % Using fminsearch:
    if(strcmp(cfg.estimation.optimizer, 'fminsearch'))
      % initializing the best objective function value to 
      % infinity, this will be updated below
      best_obj      = inf;
      best_params   = cfg.estimation.parameters.guess ;
      initial_guess = best_params;


      if(isrow(initial_guess))
        initial_guess = initial_guess';
      end

      vp(cfg,  '   Method: Nelder-Mead (fminsearch)');
     
      if(cfg.estimation.effort >1)
        vp(cfg,  'Pseudo simulated annealing selected:');
        vp(cfg,  '   Sit back and have a cup of coffee;');
        vp(cfg,  '   this is going to take a while');
     
        % seeding the random number generator
        rng(8675309)
      end
      vp(cfg,  '-----------------------------------');
     
      for est_ctr = 1:cfg.estimation.effort
         if(cfg.estimation.effort >1)
           if(est_ctr/cfg.estimation.effort <= .2)
             options =  optimset('Display',   'none', ...
                                 'maxfuneval', 300, ...
                                 'maxiter',    300);
             sigma     = .5;
             stage_ctr = 1;
           elseif(est_ctr/cfg.estimation.effort <= .5)
             options =  optimset('Display',   'none', ...
                                 'maxfuneval', 300, ...
                                 'maxiter',    300);
             sigma     = .4;
             stage_ctr = 2;
           elseif(est_ctr/cfg.estimation.effort <= .7)
             options =  optimset('Display',   'none', ...
                                 'maxfuneval', 200, ...
                                 'maxiter',    200);
             sigma     = .20;
             stage_ctr = 3;
           elseif(est_ctr/cfg.estimation.effort < 1)
             options =  optimset('Display',   'none', ...
                                 'maxfuneval', 200, ...
                                 'maxiter',    200);
             sigma     = .10 ;
             stage_ctr = 3;
           else
             options =  optimset('Display',   'Iter', ...
                                 'maxfuneval', 500, ...
                                 'maxiter',    500);
             sigma     = .0001;
             stage_ctr = 4;
           end
           
           vp(cfg,sprintf('Stage %d (sigma = %.1e) iteration %d of %d (%s)', ...
                           stage_ctr,...
                           sigma,...
                           est_ctr,...
                           cfg.estimation.effort, ...
                           datestr(now)));
             
         else
             sigma     = .0001;
         end
     
     
        
         % not trying too hard, just running fminsearch:
         [current_est, current_obj ] = ...
         fminsearch('calculate_objective', ...
                    initial_guess,         ...
                    options,               ...
                    cfg);
         
     
         % disp(sprintf('Intermediate Obj: %.6e', current_obj));
     
         % rebounding the parameters
         current_est    = bound_parameters(current_est   , ...
                                           cfg.estimation.parameters.lower_bound, ...
                                           cfg.estimation.parameters.upper_bound);
         
         % if we found a solution better than the previous best,
         % then we update our best parameter set and objective value
         if(current_obj < best_obj)
           best_params = current_est;
           best_obj    = current_obj;
           disp(sprintf('New Objective found: %.6e', best_obj));
         end
     
         % Adding noise to the best estimate so far  but only if 
         % we're not at the last stage that one starts with the
         % best guess
         if(est_ctr < cfg.estimation.effort - 1)
     
     
           initial_guess  = add_noise(best_params,  sigma.*ones(size(initial_guess)), 'log_normal');
           % rebounding guess
           initial_guess  = bound_parameters(initial_guess , ...
                                             cfg.estimation.parameters.lower_bound, ...
                                             cfg.estimation.parameters.upper_bound);
         else
           % On the last stage we start with whatever our 
           % best solution is so far, so best_params is stored
           % in initial_guess
           vp(cfg, 'On the last stage');
           initial_guess = best_params;
         end
     
      end
     
      % If we're tracking the estimation status, we'll try to 
      % dump that figure to a file as well. This will fail if 
      % user specified a different status function that doesn't 
      % create this figure with the appropriate name:
      if(isfield(cfg.estimation, 'monitor'))
      if(isfield(cfg.estimation.monitor, 'status_function'))
      if(~strcmp(cfg.estimation.monitor.status_function, ''))
        try
          figure(findobj('type', 'figure', 'name', 'Estimation Status'));
          dump_figure(sprintf('output%smonitor_estimation_progress', filesep));
          vp(cfg, 'Estimation Status figure saved to:');
          vp(cfg, sprintf('output%smonitor_estimation_progress', filesep));
        catch
          vp(cfg, 'Estimation status function was selected, but ');
          vp(cfg, 'we were unable to find a figure named ');
          vp(cfg, '''Estimation Status''');
        end
      end
      end
      end
      parameters_est = best_params;
    end




    % creating statistics_est so that if 
    % the calculation fails below, this empty 
    % variable will be returned.
    statistics_est = struct();
    try
      % We're done with the estimation, now we try to calculate the solution
      % statistics. This can be numerically tricky, so we have a fail
      % condition below:
      %calculate_objective(parameters_est, cfg)
      vp(cfg, '-----------------------------------');
      vp(cfg, 'Calculating solution statistics. Be');
      vp(cfg, 'patient this can take a while when ');
      vp(cfg, 'there are many parameters          ');
      vp(cfg, '-----------------------------------');
      statistics_est = solution_statistics(parameters_est, cfg);
      files = generate_report(parameters_est, statistics_est, cfg);
      if(strcmp(cfg.options.verbose, 'yes'))
        type(files.report_file);
      end
   
      if(strcmp(cfg.options.verbose, 'yes'))
        vp(cfg, 'If you''re happy with the results, the following');
        vp(cfg, 'can be used to update system.txt file. Just copy, ');
        vp(cfg, 'paste, and delete the previous entries');
        for parameter_idx =1:length(parameters_est)
          ptmp.set_name  = cfg.estimation.set_name;
          ptmp.value     = var2string(parameters_est(parameter_idx), 12, 5);
          ptmp.name      = cfg.estimation.parameters.names{parameter_idx};
          ptmp.units     = cfg.estimation.parameters.units{parameter_idx};
          ptmp.type      = cfg.estimation.parameters.type{parameter_idx};
          ptmp.ptype     = cfg.estimation.parameters.ptype{parameter_idx};
          %ptmp.ptype     = cfg.parameters.matrix{(strcmp(cfg.parameters.matrix(:,1) , cfg.estimation.parameters.names(parameter_idx))), end};
          ptmp.lb_number = cfg.estimation.parameters.lower_bound(parameter_idx);
          ptmp.editable  = cfg.estimation.parameters.editable{parameter_idx};
          ptmp.ub_number = cfg.estimation.parameters.upper_bound(parameter_idx);
        
          ptmp.ub_text   = num2str(ptmp.ub_number);
      
          if(ptmp.lb_number == eps)
            ptmp.lb_text   = 'eps';
          else
            ptmp.lb_text   = num2str(ptmp.lb_number);
          end

          if(strcmp(ptmp.ptype, 'variance'))
            tmp_token = '<VP>';
          else
            tmp_token = '<P> ';
          end


          if(strcmp(ptmp.ptype, 'variance') | strcmp(cfg.estimation.set_name, 'default') )
            disp(sprintf('%s %-10s  %-10s %8s %10s   %-8s %5s %10s ', ...
                          tmp_token,     ...
                          ptmp.name,     ...
                          ptmp.value,    ...
                          ptmp.lb_text,  ...
                          ptmp.ub_text,  ...
                          ptmp.units,    ...
                          ptmp.editable, ...
                          ptmp.type));
          else
            disp(sprintf('<PSET:%s:%s> %-1s', ...
                         ptmp.set_name,   ...
                         ptmp.name, ...
                         ptmp.value));
          end
        end
      end
    catch
      % If we fail to calculate the solution statistics
      % a message is displayed and the parameter values are
      % also given back.
      vp(cfg, 'Solution statistics calculation failed');
      vp(cfg, 'This can happen when you have a parameter');
      vp(cfg, 'set that makes the system stiff, or if you');
      vp(cfg, 'do not have a license for the Statistics  ');
      vp(cfg, 'Toolbox. The final parameter estimates are:');
      for parameter_idx =1:length(parameters_est)
        disp(sprintf('#--> %s = %s', ...
               cfg.estimation.parameters.names{parameter_idx}, ...
               var2string(parameters_est(parameter_idx), 8, 5)));
      end
    end

  else
  % creating an emapty parameters_est 
  % to return 
  parameters_est = [];
  end


