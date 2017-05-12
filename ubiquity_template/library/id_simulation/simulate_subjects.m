function [p] = simulate_subjects(parameters, cfg)
%function [predictions] = simulate_subjects(parameters, cfg)
%
% Inputs:
%
% cfg - System configuration variable generated in the following manner:
%
% cfg  = select_set(cfg, 'default');
%
% parameters - Full parameter vector describing the parameters. The values for
% the current parameter set can be used:
%
% parameters = cfg.parameters.values;
%
% cfg.options.stochastic - data structure with the following fields:
%
%   ci - 
%      confidence interval (default 95)
%
%   nsub -
%      number of subjects to simulate  (default 100)
%
%   seed - 
%      seed for random number generator (default 8675309)
%
%   ponly -
%      generate only the parameters and do not perform the simulation
%      'yes', or 'no' (default)
%
% These values can then be modified as necessary.
%
% Output:
%
% The predictions data structure contains the following:
%
% predictions.subjects.all - 
%   Full parameter vector (one per column) for each subject
%
% predictions.subjects.name -
%   A field for each parameter with inter-individual variability. Each field
%   contains a vector of parameter values that were used.
%
% predictions.times
%   A field for every timescale containing the sample times from the
%   simulation.
%
% predictions.states and predictions.outputs -
%   There is a field for each state or output which contains a profile for
%   each subject (one per column) and each row corresponds to the sampling
%   times in predictions.times
% 
% predictions.states_stats and predictions.outputs_stats -
%   There is a field for each state or output which contains the following
%   fields:
%      lb_ci:  lower bound of the confidence interval for that named value
%      ub_ci:  upper bound of the confidence interval for that named value
%      mean:   mean of the prediction for that named value 
%      median: median of the prediction for that named value 
%   These are all vectors corresponding to the sampling times in
%   prediction.times
%
% predictions.times_patch, predictions.states_patch and predictions.outputs_patch -
% These contain vectors to be used with the patch command to generate shaded
% regions. For example if you had an output called Coverage, the following
% would shade in the region representing the confidence interval, specified by
% options.ci = 95 (the default):
%
%   mc = fetch_color_codes;
%
%   % This line is only necessary if you have some negative values and you want
%   % to put it on a log scale:
%   predictions.outputs_patch.Coverage.ci(predictions.outputs_patch.Coverage.ci < 0) = eps;
%   patch(predictions.times_patch.weeks, ...
%         predictions.outputs_patch.Cp_Total.ci, ...
%         mc.light_blue, 'edgecolor', 'none');
%
% To all of the data for the Coverage output, use the following: 
% plot(predictions.times.weeks, ...
%      predictions.outputs.Coverage, '.', 'color', 'r');
%
% This plots the mean:
%   plot(predictions.times.weeks, ...
%        predictions.outputs_stats.Cp_Total.mean, 'b-');
%
% This plots the upper and lower confidence intervals:
%   plot(predictions.times.weeks, ...
%        predictions.outputs_stats.Cp_Total.lb_ci, 'r--');
%   plot(predictions.times.weeks, ...
%        predictions.outputs_stats.Cp_Total.ub_ci, 'r--');
%   

% Parsing options
% First we define the defaults
nsub         = 100;
seed         = 8675309;
ci           = 95;
ponly        = false;
state_names  = fieldnames(cfg.options.mi.states);
output_names = fieldnames(cfg.options.mi.outputs);

% Then we overwrite them based on what's in cfg

if(isfield(cfg.options.stochastic, 'nsub'))
  nsub = cfg.options.stochastic.nsub;
end

if(isfield(cfg.options.stochastic, 'seed'))
  seed = cfg.options.stochastic.seed;
end

if(isfield(cfg.options.stochastic, 'ci'))
  ci   = cfg.options.stochastic.ci;
end

if(isfield(cfg.options.stochastic, 'ponly'))
  ponly   = cfg.options.stochastic.ponly;
end
if(isfield(cfg.options.stochastic, 'states'))
  state_names = cfg.options.stochastic.states;
end
if(isfield(cfg.options.stochastic, 'outputs'))
  output_names = cfg.options.stochastic.outputs;

  % By default all outputs will include those with and without residual error.
  % If the user specifies outputs manually, then we also add in the output
  % with error if it has been defined in the system file.
  output_names_specified = output_names;
  for output_idx = 1:length(output_names_specified)
    output_name = output_names_specified{output_idx};
    if(isfield(cfg.ve, output_name))
      output_names(end+1) = {sprintf('SIOE_%s', output_name)};
    end
  end

end




% initializing the output
p.subjects.all  = struct();
p.subjects.name = struct();

if(not(ponly))
  p.states        = struct();
  p.outputs       = struct();
  p.outputs_stats = struct();
  p.states_stats  = struct();
  
  p.times_patch   = struct();
  p.states_patch  = struct();
  p.outputs_patch = struct();
end


max_errors = 100;

isgood = 1;


vp(cfg, sprintf('Simulating multiple subjects (%d)', nsub));
vp(cfg, sprintf('Integrating with: %s', cfg.options.simulation_options.integrate_with));
strTP = textprogressbar('  ', []);


% checking to make sure we've specified the iiv terms
if(isfield(cfg.iiv, 'parameters'))
  % trying to generate a random sample from the variance covariance matrix:
  if(min(eig((cfg.iiv.values + cfg.iiv.values')./2)) <=0)
    vp(cfg, '---------------------------------------------- ');
    vp(cfg, ' simulate_subjects.m                          ');
    vp(cfg, ' Warning: The variance/covariance matrix is not ');
    vp(cfg, ' positive semi-definite. Testing only the diagonal ');
    vp(cfg, ' elements. I.e. no covariance/interaction terms');
    try
      cfg.iiv.values = diag(diag(cfg.iiv.values));
      mvnrnd(zeros(1, length(cfg.iiv.values(:,1))), ...
             cfg.iiv.values, 1);
      vp(cfg, ' Using only the diagional elements seems to   ');
      vp(cfg, ' have worked. Understand that the results do  ');
      vp(cfg, ' not include any interaction.                 ');
  
    catch
      vp(cfg, ' Failed using only diagonal/variance elements.');
      vp(cfg, ' Check the specified IIV elements in');
      vp(cfg, ' cfg.iiv.values');
      isgood = 0;
    end
    vp(cfg, '---------------------------------------------- ');
  end
  
  % Seeding the random number generator
  rng(seed);
  
  
  %
  % simulating for each subject
  %
  sub_idx = 1;
  error_count = 0;
  while (sub_idx <= nsub) & isgood
    
    % Generating a subject:
    subject = generate_subject(parameters,  cfg);
    parameters_subject = subject.parameters.all;

  
    % attempting to simulate that subject if this fails we try again :)
    try
      % printing out subject every 10 percent 
      strTP = textprogressbar(sub_idx/nsub.*100, strTP);
      if(not(ponly))
        som          = run_simulation_ubiquity(parameters_subject, cfg);
      
      % Some things are going to be the same for all simulations
      % so we pull that information out of the first simulation
        if(sub_idx == 1)
          p.times      = som.times;
        
          % for each set of time units
          % we create a vector to be used in
          % generating patches
          time_names = fieldnames(p.times);
          for time_idx = 1:length(time_names)
            time_name = time_names{time_idx};
            time_values = getfield(p.times, time_name);
            eval(sprintf('p.times_patch.%s = [time_values; flipud(time_values)];', time_name));
          end
        end
        
        % creating a summary table for each state and output
        for state_idx = 1:length(state_names)
          state_name = state_names{state_idx};
          if(isfield(p.states, state_names))
            eval(sprintf('p.states.%s(:,end+1) = som.states.%s;', state_name, state_name));
          else
            eval(sprintf('p.states.%s = som.states.%s;', state_name, state_name));
          end
        end
        for output_idx = 1:length(output_names)
          output_name = output_names{output_idx};
          if(isfield(p.outputs, output_names))
            eval(sprintf('p.outputs.%s(:,end+1) = som.outputs.%s;', output_name, output_name));
          else
            eval(sprintf('p.outputs.%s = som.outputs.%s;', output_name, output_name));
          end
        end
      end
  
      % since we're successful up to this point
      % we're saving the subject parameter data
      pnames = fieldnames(subject.parameters.name);
      for p_idx = 1:length(pnames)
        pname = pnames{p_idx};
        if(sub_idx == 1)
          eval(sprintf('p.subjects.name.%s = subject.parameters.name.%s;', pname, pname));
        else
          eval(sprintf('p.subjects.name.%s = [p.subjects.name.%s, subject.parameters.name.%s];', pname, pname, pname));
        end
      end
      if(sub_idx == 1)
        eval(sprintf('p.subjects.all     = parameters_subject;'));
      else
        eval(sprintf('p.subjects.all     = [p.subjects.all     parameters_subject];'));
      end
  
      sub_idx = sub_idx + 1;
    catch TRYERROR 


      % When we hit the first error we log it:
      if(error_count == 1)
        system_log_tryerror(cfg, TRYERROR);
      end
      %system_log_tryerror(cfg, 'further errors will be suppressed');

      simerror = true;
      if (error_count < max_errors)
        % we had an error so we just add that to the count
        error_count = error_count + 1;
      else
        % we've reached the maximum number of errors
        % we're going to issue a warning and exit the 
        % for loop
        disp(sprintf('Warning we''ve reached the maximum'));
        disp(sprintf('number of errors (%d)', max_errors));
        disp(sprintf('Make sure that individual simulations can run succesfully'));
        disp(sprintf('Exiting, simulation of subjects failed.'));
        sub_idx = nsub  + 2;
      end
    end
  end
  
  
  if(not(ponly))
    % Calculating the statistics for the output
    % (mean, median, CI etc. for each time point)
    for output_idx = 1:length(output_names)
      output_name = output_names{output_idx};
      output_data = getfield(p.outputs, output_name);
      [tcs, tcp] = timecourse_stats(output_data, ci);
      eval(sprintf('p.outputs_stats.%s = tcs;', output_name));
      eval(sprintf('p.outputs_patch.%s = tcp;', output_name));
    end
    
    % Calculating the statistics for the states
    % (mean, median, CI etc. for each time point)
    for state_idx = 1:length(state_names)
      state_name = state_names{state_idx};
      state_data = getfield(p.states, state_name);
      [tcs, tcp] = timecourse_stats(state_data, ci);
      eval(sprintf('p.states_stats.%s = tcs;', state_name));
      eval(sprintf('p.states_patch.%s = tcp;', state_name));
    end
  end
else
    disp('---------------------------------------------- ');
    disp('  simulate_subjects.m                          ');
    disp('> Error:Trying to simulate subjects with       ');
    disp('>    variability, but no variance/covariance   ');
    disp('>    information was specified.                ');
    disp('>                                              ');
    disp('>    Modify the system.txt file to add the     ');
    disp('>    IIV information using the following:      ');
    disp('>     <IIV:?>      ?                           ');
    disp('>     <IIV:?:?>    ?                           ');
    disp('>     <IIVCOR:?:?> ?                           ');
    disp('---------------------------------------------- ');

end

textprogressbar(' -> Done', strTP);


function [tcs, tcp] = timecourse_stats(d, ci)
%
% Calculates stats across a matrix of outputs where each column is a subject
% and each row represents a sample time
%


 myci = ci/100;
 dsorted = sort(d,2);
 nsubs   = length(dsorted(1,:));
 lb_idx  = nsubs*(1-myci)/2 + 1;
 ub_idx  = nsubs - nsubs*(1-myci)/2;

 lb_ci   = mean(dsorted(:,[floor(lb_idx), ceil(lb_idx)]),2);
 ub_ci   = mean(dsorted(:,[floor(ub_idx), ceil(ub_idx)]),2);

 tcs.lb_ci  = lb_ci;
 tcs.ub_ci  = ub_ci;
 tcs.mean   = mean(dsorted, 2);
 tcs.median = median(dsorted, 2);

 tcp.ci  = [ub_ci; flipud(lb_ci)];  


function [subjects] = generate_subject(parameters, cfg)
% function [subjects] = generate_subject(parameters, cfg)
%
% Generates subjects with variability specified using the <IIV:?> descriptor
% in the system.txt file
%
% Inputs:
%
% cfg - system configuration variable generated in the following manner:
%
% % selecting the default parameter set:
% cfg  = select_set(cfg, 'default');
%
% parameters - vector of typical parameter values. This can be obtained from
% the cfg variable:
%
% parameters = cfg.parameters.values;
%
% This can be modified before subject generation
%
% Output:
%
% The data structure 'subjects' will be generated with the following fields:
%
% subjects.parameters.all  - matrix where each column represents a sample
% subjects.parameters.name - each parameter with variability will have a field
% under 'name'. This is a vector of sampled values.
%



subjects.parameters.all   = [];
subjects.parameters.name  = [];

%
% Generating the subjects
%
iiv_parameter_names = fieldnames(cfg.iiv.parameters);
% creating a temporary vector containing the typical values of all of the
% parameters:
TMP_parameters_all = parameters;

% defining the mean of the IIVs and the covariance matirx
covmatrix = cfg.iiv.values;
muzero = zeros(1, length(covmatrix(:,1)));


iiv_sample = mvnrnd(muzero,covmatrix);

% now looping through each parameter with inter-individual variability
for(pidx=1:length(iiv_parameter_names))

  % use mvrnd to generate to account for off diagional effects
  % http://radio.feld.cvut.cz/matlab/toolbox/stats/mvnrnd.html;

  %sampling from the IIV parameter distribution

  % getting name and typical value of the current parameter with specified IIV:
  TMP_parameter_name  = iiv_parameter_names{pidx};
  TMP_parameter_value = parameters(getfield(cfg.options.mi.parameters, TMP_parameter_name));

  % getting the type of distribution for this parameter:
  % TMP_distribution   = getfield(getfield(cfg.iiv.parameters, TMP_parameter_name), 'distribution');
  TMP_equation       = getfield(getfield(cfg.iiv.parameters, TMP_parameter_name), 'equation');

  % getting the IIV parameter associated with this parameter and the value
  % from the variance/covariance matrix:
  TMP_iiv_name       = getfield(getfield(cfg.iiv.parameters, TMP_parameter_name), 'iiv_name');

  % IIV (ETA) sample       = TMP_iiv_value
  % Typical value (THETA)  = TMP_parameter_value
  eval(sprintf('TMP_iiv_value = iiv_sample(cfg.options.mi.iivs.%s);',  TMP_iiv_name));

  TMP_subject_parameter_value = generate_parameter(parameters, cfg, TMP_parameter_value, TMP_iiv_value, TMP_equation);

  % storing the samples in the temporary vector with all parameters:
  TMP_parameters_all(getfield(cfg.options.mi.parameters, TMP_parameter_name)) = TMP_subject_parameter_value;

  % storing the samples in the specific parameter vector:
  eval( sprintf('subjects.parameters.name.%s =  [TMP_subject_parameter_value];', ...
                 TMP_parameter_name, TMP_parameter_name));
end

% storing the full vector in the return value
subjects.parameters.all(:,end+1)   = TMP_parameters_all;


function [SIMINT_PARAMETER_SAMPLE] = generate_parameter(SIMINT_parameters, SIMINT_cfg, SIMINT_PARAMETER_TV, SIMINT_IIV_VALUE, SIMINT_equation)    
% Generates the parameter value based on the IIV and the equation specified
% for that particular distribution

% Defining the parameters locally
% this is for distributions that depend on the IIV and 
% system parameters like a Box-Cox distribution:
for SIMINT_parameter_idx = 1:length(SIMINT_parameters)
  % this defines the parameter
  eval(sprintf('%s = SIMINT_parameters(SIMINT_parameter_idx);',SIMINT_cfg.parameters.names{SIMINT_parameter_idx}));
end

% Evaluating the parameter sample
eval(sprintf('SIMINT_PARAMETER_SAMPLE = %s;', SIMINT_equation));

