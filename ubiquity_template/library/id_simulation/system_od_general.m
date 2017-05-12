function [od, odp]=system_od_general(pest, cfg)
% function [od, odp]=system_od_general(pest, cfg)
%
% Generates observation details for cohorts specified by the
% system_define_cohort() function. 
%
%
% inputs:
% pest - Vector of parameters being estimated
%
% outputs:
% od  - Observation details data structure used for calculating 
%       the estimation objective.
%
% odp - Data structure containing observation and prediction 
%       information for each output in each cohort. This is used to 
%       generate figures with smooth profiles.
%
%
% For a given cohort (CHNAME) and output (OPNAME) you can access the 
% observations and corresponding times using the following:
%
%   odp.cohorts.CHNAME.OPNAME.od.obs
%   odp.cohorts.CHNAME.OPNAME.od.time
%
% The smooth model predictions are stored here:
%
%   odp.cohorts.CHNAME.OPNAME.od.time_smooth
%   odp.cohorts.CHNAME.OPNAME.od.pred_smooth
%
% to add:
%  overwrite parameter set at the cohort level

% Format of cohort in cfg:
% cfg.cohorts.ch_name.cf.route                        = 1.0;
% cfg.cohorts.ch_name.cf.dose_level                   = 0.1;
% cfg.cohorts.ch_name.cp.PNAME                        = VALUE;
% cfg.cohorts.ch_name.dataset                         = 'ds_name';
% cfg.cohorts.ch_name.observation_times               = [];                % derived 
% cfg.cohorts.ch_name.outputs.pk.obs.time             = 'obs_time_col';
% cfg.cohorts.ch_name.outputs.pk.obs.value            = 'obs_obs_col';
% cfg.cohorts.ch_name.outputs.pk.obs.missing          = -1;
% cfg.cohorts.ch_name.outputs.pk.model.time           = 'TS  ';
% cfg.cohorts.ch_name.outputs.pk.model.value          = 'MODOUT';
% cfg.cohorts.ch_name.outputs.pk.model.variance       = 'PRED';  
% cfg.cohorts.ch_name.outputs.pk.data.time            = [];                % derived  
% cfg.cohorts.ch_name.outputs.pk.data.obs             = [];                % derived   
% cfg.cohorts.ch_name.outputs.pk.data.simtime         = [];                % derived   
% cfg.cohorts.ch_name.outputs.pk.options.marker_color = 'r';
% cfg.cohorts.ch_name.outputs.pk.options.marker_shape = 'o';
% cfg.cohorts.ch_name.outputs.pk.options.marker_line  = '-'; 
% cfg.cohorts.ch_name.inputs.bolus.Cp.TIME            = 0.0;
% cfg.cohorts.ch_name.inputs.bolus.Cp.AMT             = 0.1;


odp = struct();

% Getting all of the cohorts
cohorts = fieldnames(cfg.cohorts);

% storing the names of the outputs
opall = struct();


for(chidx=1:length(cohorts))

  % Making a local cohort-specific copy of cfg
  chcfg = cfg;

  % Pulling out current cohort
  cohort = getfield(cfg.cohorts, cohorts{chidx});

  % Smooth output times
  % by default the output times will be those for the simulation
  choutput_times = cfg.options.simulation_options.output_times;

  % If this cohort has a different set of output times then 
  % we overwrite the defaults
  if(isfield(cohort.options, 'output_times'))
    choutput_times = cohort.options.output_times;
  end

  % If the output times are specified as a row we turn it into a column
  if(isrow(choutput_times))
     choutput_times = choutput_times';
  end

  % Adding all of the observation times to the output times to make sure the
  % simulations evaluate at these times
  choutput_times = sort(unique([choutput_times; cohort.observation_simtimes]));

  % Setting times to give a smooth profile, this will include the cohort
  % output times as well 
  chcfg = system_set_option(chcfg, 'simulation', 'output_times', choutput_times);

  % Getting the full parameter vector
  chparameters = fetch_full_parameters(pest, chcfg);

  % Overwriting cohort specific parameters
  if(isfield(cohort, 'cp'))
    pnames = fieldnames(cohort.cp);
    for pidx = 1:length(pnames)
      [chparameters] = system_set_parameter(chcfg, chparameters, ...
                                            pnames{pidx},  ...
                                            getfield(cohort.cp, pnames{pidx}));
    end
  end

  %
  % Setting up the inputs  
  %
  % zeroing out all events
  chcfg=system_zero_inputs(chcfg);

  % Bolus inputs:
  if(isfield(cohort, 'inputs'))
    if(isfield(cohort.inputs, 'bolus'))
      states = fieldnames(cohort.inputs.bolus);
      for(stidx = 1:length(states))
        chcfg=system_set_bolus(chcfg, ...
            states{stidx}, ...
            getfield(getfield(cohort.inputs.bolus, states{stidx}), 'TIME'), ... 
            getfield(getfield(cohort.inputs.bolus, states{stidx}), 'AMT'));     
      end
    end
    
    % Infusion rates:
    if(isfield(cohort.inputs, 'infusion_rates'))
      irates = fieldnames(cohort.inputs.infusion_rates);
      for(iridx = 1:length(irates))
        chcfg=system_set_rate(chcfg, ...
            irates{iridx}, ...
            getfield(getfield(cohort.inputs.infusion_rates, irates{iridx}), 'TIME'), ... 
            getfield(getfield(cohort.inputs.infusion_rates, irates{iridx}), 'AMT'));     
      end
    end
    % Covariates:    
    if(isfield(cohort.inputs, 'covariates'))
      covars = fieldnames(cohort.inputs.covariates);
      for(cvidx = 1:length(covars))
        chcfg=system_set_covariate(chcfg, ...
            covars{cvidx}, ...
            getfield(getfield(cohort.inputs.covariates, covars{cvidx}), 'TIME'), ... 
            getfield(getfield(cohort.inputs.covariates, covars{cvidx}), 'AMT'));     
      end
    end
  end

  %
  % Setting up the outputs 
  %

  % Pulling out the dataset 
  % for the cohort (chds)
  chds  = getfield(chcfg.data, cohort.dataset);
  % selecting just the records for this cohort
  %chvalues     = nm_select_records(chds, chds.values, cohort.cf);
  outputs = fieldnames(cohort.outputs);

  % Simulating the cohort  
  som = run_simulation_ubiquity(chparameters, chcfg);

  % Sampling the different outputs for this cohort
  for(opidx = 1:length(outputs))
    opall         = setfield(opall, outputs{opidx}, '');
    output        = getfield(cohort.outputs, outputs{opidx});
    %odchunk       = getfield(odchunks, outputs{opidx});

    % Pulling out the data for this cohort/output combination
    odchunk       =  getfield(getfield(cohort.outputs, outputs{opidx}), 'data');
    odchunk.pred  = interp1(getfield(som.times,   output.model.time),  ...
                           getfield(som.outputs, output.model.value), ...
                           odchunk.time,  'linear');

    % Getting the time and predictions corresponding to the observations
    % specified in the dataset:
    odchunk.time_smooth = getfield(som.times,   output.model.time);
    odchunk.pred_smooth = getfield(som.outputs, output.model.value);

    % Calculating the variance
    odchunk.var = calculate_variance(chparameters, output.model.variance, odchunk, chcfg);

    % We're storing the observation details here:
    eval(sprintf('od.%s.%s = [odchunk.time, odchunk.obs, odchunk.pred, odchunk.var];', cohorts{chidx}, outputs{opidx}));

    % Storing the smooth profiles for plotting
    eval(sprintf('odp.cohorts.%s.%s.od   = odchunk;', cohorts{chidx}, outputs{opidx}));
    eval(sprintf('odp.cohorts.%s.%s.som  = som;',     cohorts{chidx}, outputs{opidx}));
  end
end

% stores a list of all outputs which is iterated through during modeling
odp.meta.outputs = fieldnames(opall);

function [SIMINT_var] = calculate_variance(SIMINT_parameters, SIMINT_varstr, SIMINT_odchunk, SIMINT_cfg)
%
% Calculating the variance for a given cohort/output combination 
%


% Defining the system and variance parameters locally
for(SIMINT_idx = 1:length(SIMINT_cfg.parameters.names))
    SIMINT_PNAME = SIMINT_cfg.parameters.names{SIMINT_idx};
    eval(sprintf('%s  = SIMINT_parameters(SIMINT_idx);', SIMINT_PNAME));
end

% Defining the pred, time and obs locally
SIMINT_pred = SIMINT_odchunk.pred;
SIMINT_time = SIMINT_odchunk.time;
SIMINT_obs  = SIMINT_odchunk.obs ;

% vectorize the equation in varstr
SIMINT_varstr = strrep(SIMINT_varstr, '^', '.^');
SIMINT_varstr = strrep(SIMINT_varstr, '/', './');
SIMINT_varstr = strrep(SIMINT_varstr, '*', '.*');

% substituing the locally defined pred, time and obs
SIMINT_varstr = strrep(SIMINT_varstr, 'PRED', 'SIMINT_pred') ;
SIMINT_varstr = strrep(SIMINT_varstr, 'TIME', 'SIMINT_time') ;
SIMINT_varstr = strrep(SIMINT_varstr, 'OBS' , 'SIMINT_obs' ) ;

% calculating the variance
eval(sprintf('SIMINT_var = (%s);', SIMINT_varstr));

%
% If the variance is constant then it will be evaluated 
% as a scalar and we will need to make sure it has the 
% the same dimensions the PRED vector
%
if(length(SIMINT_var) == 1 & length(SIMINT_pred) > 1)
  SIMINT_var = SIMINT_var*ones(size(SIMINT_pred));
end
