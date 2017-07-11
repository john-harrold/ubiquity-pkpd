function [cfg]=system_define_cohort(cfg, cohort)
%function [cfg]=system_define_cohort(cfg, cohort)
%
%  cohort - Data structure describing the cohort being added. This contains
%  information about the dataset to be used, how the dataset should be
%  filtered to return only data for that cohort, the outputs being used,
%  inputs (dosing, infusions, etc) for this cohort, etc.
%
%  Naming the cohort:
%  ------------------
%  cohort.name - This is a required field that contains a string that must
%  begin with a letter and can contain letters, numbers and underscores. This
%  should be short and descriptive. 
%
%  Dataset name
%  ------------
%  cohort.dataset - Name of the data set defined using system_load_dataset().
%
%  Overwriting parameters: 
%  -----------------------
%  cohort.cp - This is a data structure with fields that are parameter names
%  and values that these parameters should be set for this specific cohort.
%
%
%  Filtering out the cohort data:
%  ------------------------------
%  cohort.cf - This is a filter that is applied to the dataset to return the
%  subset that applies to the current cohort. Where the fields of cohort.cf
%  are the column names of the dataset are 'anded' together to extract this
%  subset. Say there are columns in the dataset representing the genotype
%  (gtype) and the dose level (dose). And we are looking at the cohort where
%  the genotype is 1 and the dose is 1000. Then the following would be used:
%
%     cohort.cf.gtype = 1;
%     cohort.cf.dose  = 1000;
%
%  Similarly, if we wanted to group all of the cohorts given a dose of 1000
%  where the genotype was 1 or 3. We would do the following:
%
%     cohort.cf.gtype = [1, 3];
%     cohort.cf.dose  = 1000;
%
%  Inputs 
%  ------
%  cohort.inputs - There is a field here for each type of input (bolus,
%  infusion_rates, and covariates). Each of these inputs has the same format.
%  A field for the input name and under that input name there is an TIME and
%  AMT fields. These are vectors of times and amounts with the units described
%  in the system.txt file. It is only necessary to specify the bolus and
%  infusion rates that are nonzero for this cohort and the covariate that is
%  different from the parameter set for this cohort.
%  
%  >Bolus<
%
%  cohort.inputs.bolus.STATE.TIME          = [];
%  cohort.inputs.bolus.STATE.AMT           = [];
%  
%  where STATE is the state or compartment receiving the bolus. 
%  
%  >Infusion Rates<
%
%  cohort.inputs.infusion_rates.RNAME.TIME = [];
%  cohort.inputs.infusion_rates.RNAME.AMT  = [];
%  
%  where RNAME is the name specified by the <R:RNAME> descriptor.
%  
%  >Covariates<
%  
%  cohort.inputs.covariates.CNAME.TIME     = [];
%  cohort.inputs.covariates.CNAME.AMT      = [];
%  
%  where CNAME is the name specified by the <CV:CNAME> descriptor.
%     
%  Outputs
%  -------
%  cohort.outputs - There is a field for each output here. This field, like
%  the name for the cohort, should begin with a letter and can contain only
%  letters, numbers and underscores. For each output there is an 'obs' field
%  that describes the columns in the dataset to be used that contain both the
%  independent variable (time) and the dependent variable (output). If there
%  is missing data, an optional field of 'missing' can be specified and the
%  value used to indicate a missing observation can be specified (e.g. -1).
%
%  There is also a field called 'model' which contains the timescale from the
%  model which matches the 'time' in the dataset. It also contains a 'value'
%  field that maps the model output to the observations in the dataset.
%
%  Options 
%  -------
%  cohort.options - This field contains optional information used to customize
%  how this cohort is simulated and what is returned. 
%
%  cohort.options.output_times - By default, the simulation output times set
%  by the following command: 
%
%      cfg = system_set_option(cfg, 'simulation', ...
%                                   'output_times', ...
%                                   linspace(0,100,1000)');
%
%  If it is necessary to have finer control over this, a different vector of
%  times can be specified. 
%
%  Example:
%  --------
%  First we would want to make sure we clear any previously defined cohorts:
%
%  cfg = system_clear_cohorts(cfg);
%
%  Say we have a cohort in the mouse_pk dataset that consists of female (sex=2
%  in dataset) mice, that we wish to name mouse_female, and these mice were in
%  the 5mpk dosing group (dose_level=5 in dataset).
%
%  cohort.name                              = 'mouse_female';                  
%  cohort.label                             = 'Female Mice';  % optional
%  cohort.cf.dose_level                     = 5;           
%  cohort.cf.sex                            = 2;           
%  cohort.dataset                           = 'mouse_pk'; 
%
%  Next we need to match the outputs in the model to the outputs in the
%  dataset. We'll call the serum pk output simply 'pk'. The times and
%  observations in the dataset are found in the 'time_days' column and the
%  'obs' column (missing data specified by -1). These are mapped to the model
%  outputs (which MUST have the same units) 'days' and 'Cp_ng_ml'. 
%
%  cohort.outputs.pk.obs.time               = 'time_days';% column in dataset
%  cohort.outputs.pk.obs.value              = 'obs';      % column in dataset
%  cohort.outputs.pk.obs.missing            = -1;         % optional 
%  cohort.outputs.pk.model.time             = 'days';     % model timescale <TS:?>
%  cohort.outputs.pk.model.value            = 'Cp_ng_ml'; % model output <O>
%
%  You must also specify the variance model to use for this cohort/output
%  combination. This is a string that can  be '1' for least squares
%  estimation. Or a mathematical formula using any of the system parameters
%  and PRED, TIME, and OBS for the model prediction, sample time and
%  observation respectively. Below are some examples:
%
%  cohort.outputs.pk.model.variance         = '1';             
%  cohort.outputs.pk.model.variance         = 'PRED^2';        
%  cohort.outputs.pk.model.variance         = 'SLOPE*PRED^2';  
%
%  For each output you can optionally specify a color, shape and line type for
%  markers to be used when plotting. 
%
%  cohort.outputs.pk.options.marker_color   = 'r';        % optional 
%  cohort.outputs.pk.options.marker_shape   = 'o';        % optional 
%  cohort.outputs.pk.options.marker_line    = '-';        % optional 
%
%  If this cohort has multiple outputs, it may be necessary to further filter
%  the data for this output. If the dataset has a column called cmt that has a
%  value of 1 for the output pk, then we would add the following:
%
%  cohort.outputs.pk.of.cmt                  = 1;         % column in dataset  
%
%  Now we need to define the inputs. If there are two doses of 5mpk applied
%  one week apart (and it is defined this way in the system.txt file). Then 
%  the following would define this for the cohort:
%
%  cohort.inputs.bolus.Cp.AMT               = [5 5];           
%  cohort.inputs.bolus.Cp.TIME              = [0 1];           
%
%  Similar inputs statements can be made for infusion rates and covariates.
%  Where the AMT and TIME have the same units as the definitions in the
%  system.txt:
%
%  %cohort.inputs.infusion_rates.RNAME.AMT   = [];           
%  %cohort.inputs.infusion_rates.RNAME.TIME  = [];           
%  
%  %cohort.inputs.covariates.CNAME.AMT       = [];           
%  %cohort.inputs.covariates.CNAME.TIME      = [];           
%  
%  Next we add this cohort to the cfg variable:
%  cfg = system_define_cohort(cfg, cohort);


% Making sure the options field exists
if(~(isfield(cohort, 'options')))
  cohort.options = struct();
end

% default output options for an output cohort:
defopts.marker_color   = 'k';           
defopts.marker_shape   = 'o';           
defopts.marker_line    = '--';           

validopts = {'marker_color'
             'marker_shape'
             'marker_line'};

isgood      = true;
datasetgood = true;

% checking the name
if(isfield(cohort, 'name'))
  if(isfield(cfg.cohorts, cohort.name))
    isgood = false;
    vp(cfg, sprintf('Error: cohort with name >%s< has already been defined', cohort.name));
  else
    cohort_name = cohort.name;
    name_check = ubiquity_name_check(cohort.name);
    if(~name_check.isgood)
      isgood = false;
      vp(cfg, sprintf('Error: cohort with name >%s< is invalid', cohort.name))
      vp(cfg, sprintf('Problems: %s', name_check.msg))
    end
  end 



else
  isgood = false;
  vp(cfg, 'Error: cohort name not specified');
  cohort_name = 'no name specified';
end

% checking to see if the dataset is loaded
if(isfield(cohort, 'dataset'))
  if(~isfield(cfg.data, cohort.dataset))
    isgood      = false;
    datasetgood = false;
    vp(cfg, sprintf('Error: dataset >%s< not found, please load first', cohort.dataset));
  else
    % pulling the dataset out to test for fields below
    tmpdataset = getfield(cfg.data, cohort.dataset);
  end 
else
  isgood      = false;
  datasetgood = false;
  vp(cfg, 'Error: dataset not specified for the cohort');
end

% If we have cohort specific parameters, we check to make sure they exist
if(isfield(cohort, 'cp'))
  pnames = fieldnames(cohort.cp);
  for(pidx =1:length(pnames))
    if(~isfield(cfg.options.mi.parameters, pnames{pidx}))
      isgood      = false;
      vp(cfg, sprintf('Error: The parameter >%s< ', pnames{pidx}));
      vp(cfg, sprintf('       is not defined. Check the spelling'));
      vp(cfg, sprintf('       or define the parameter using <P> '));
    else
      if(isfield(cfg.estimation.mi, pnames{pidx}))
        isgood      = false;
        vp(cfg, sprintf('Error: The parameter >%s< ', pnames{pidx}));
        vp(cfg, sprintf('       is selected for estimation. It is '));
        vp(cfg, sprintf('       not possible to fix a parameter   '));
        vp(cfg, sprintf('       that is being estiamted.          '));
      end
    end
  end
end

% checking the cohort filter columns against the columns in the dataest
if(datasetgood)
  if(isfield(cohort, 'cf'))
    cffields = fieldnames(cohort.cf);
    % looping through each field making sure the column exists in the dataset
    for(cfidx = 1:length(cffields))
      if(~isfield(tmpdataset.column.names, cffields{cfidx}))
        isgood = false;
        vp(cfg, sprintf('Error: The column >%s< in the cohort filter ', cffields{cfidx}));
        vp(cfg, sprintf('       was not found in the data set >%s< ', cohort.dataset));
      end
    end
  else
    % creating empty cohort filter
    cohort.cf = struct();
    vp(cfg, sprintf('Warning: No cohort filter was specified.'));
  end
end


%
% Checking the inputs
%
if(isfield(cohort, 'inputs'))

  %
  % Bolus input
  %
  if(isfield(cohort.inputs, 'bolus'))
    binputs  = fieldnames(cohort.inputs.bolus);
    if(isfield(cfg.options.inputs, 'bolus'))
      for(bidx =1:length(binputs))
        % Making sure the specified state has a bolus definition
        if( isfield(cfg.options.inputs.bolus.species, binputs{bidx}))
          binput = getfield(cohort.inputs.bolus, binputs{bidx});
          if(isfield(binput, 'AMT') & isfield(binput, 'TIME'))
            if(length(binput.AMT) ~= length(binput.TIME))
              isgood = false;
              vp(cfg, sprintf('Error: For the bolus input >%s< the length of ', binputs{bidx}));
              vp(cfg, sprintf('       the AMT and TIME fields need to be the same'));
            end
          else
            isgood = false;
            vp(cfg, sprintf('Error: The bolus input >%s< needs an ''AMT'' and a ''TIME'' field', binputs{bidx}));
            vp(cfg, sprintf('       cohort.inputs.bolus.%s.AMT  = []', binputs{bidx}));
            vp(cfg, sprintf('       cohort.inputs.bolus.%s.TIME = []', binputs{bidx}));
          end
        else
          isgood = false;
          vp(cfg, sprintf('Error: The bolus input >%s< has not been defined for this system', binputs{bidx}));
          vp(cfg, sprintf('       <B:times>;  %s  []; scale; units', pad_string('', length(binputs{bidx}))));
          vp(cfg, sprintf('       <B:events>; %s; []; scale; units', binputs{bidx}));
        end
      end
    else
      isgood = false;
      vp(cfg, sprintf('Error: An bolus input was specified for this cohort but'));
      vp(cfg, sprintf('       there are no bolus inputs defined in the system.txt file.'));
      vp(cfg, sprintf('       <B:times>;         []; scale; units'));
      vp(cfg, sprintf('       <B:events>; STATE; []; scale; units'));
    end
  end
  %
  % Infusion rates
  %
  if(isfield(cohort.inputs, 'infusion_rates'))
    rinputs  = fieldnames(cohort.inputs.infusion_rates);
    if(isfield(cfg.options.inputs, 'infusion_rates'))
      for(bidx =1:length(rinputs))
        % Making sure the specified infusion has a definition
        if( isfield(cfg.options.inputs.infusion_rates, rinputs{bidx}))
          rinput = getfield(cohort.inputs.infusion_rates, rinputs{bidx});
          if(isfield(rinput, 'AMT') & isfield(rinput, 'TIME'))
            if(length(rinput.AMT) ~= length(rinput.TIME))
              isgood = false;
              vp(cfg, sprintf('Error: For the infusion rate >%s< the length of ', rinputs{bidx}));
              vp(cfg, sprintf('       the AMT and TIME fields need to be the same'));
            end
          else
            isgood = false;
            vp(cfg, sprintf('Error: The infusion_rate >%s< needs an ''AMT'' and a ''TIME'' field', rinputs{bidx}));
            vp(cfg, sprintf('       cohort.inputs.infusion_rate.%s.AMT  = []', rinputs{bidx}));
            vp(cfg, sprintf('       cohort.inputs.infusion_rate.%s.TIME = []', rinputs{bidx}));
          end
        else
          isgood = false;
          vp(cfg, sprintf('Error: The infusion rate >%s< has not been defined for this system', rinputs{bidx}));
          vp(cfg, sprintf('       <R:%s>; times;   [];    scale; units ', rinputs{bidx}));
          vp(cfg, sprintf('       <R:%s>; levels;  [];    scale; units ', rinputs{bidx}));
        end
      end
    else
      isgood = false;
      vp(cfg, sprintf('Error: An infusion rate was specified for this cohort but'));
      vp(cfg, sprintf('       there are no infusion rates defined in the system.txt file.'));
      vp(cfg, sprintf('       <R:RNAME>; times;   [];    scale; units '));
      vp(cfg, sprintf('       <R:RNAME>; levels;  [];    scale; units '));
    end
  end

  %
  % Covariates
  %
  if(isfield(cohort.inputs, 'covariates'))
    cinputs  = fieldnames(cohort.inputs.covariates);
    if(isfield(cfg.options.inputs, 'covariates'    ))
      for(bidx =1:length(cinputs))
        % Making sure the specified covariate has a definition
        if( isfield(cfg.options.inputs.covariates, cinputs{bidx}))
          rinput = getfield(cohort.inputs.covariates, cinputs{bidx});
          if(isfield(rinput, 'AMT') & isfield(rinput, 'TIME'))
            if(length(rinput.AMT) ~= length(rinput.TIME))
              isgood = false;
              vp(cfg, sprintf('Error: For the covariate >%s< the length of ', cinputs{bidx}));
              vp(cfg, sprintf('       the AMT and TIME fields need to be the same'));
            end
          else
            isgood = false;
            vp(cfg, sprintf('Error: The covariate >%s< needs an ''AMT'' and a ''TIME'' field', cinputs{bidx}));
            vp(cfg, sprintf('       cohort.inputs.covariates.%s.AMT  = []', cinputs{bidx}));
            vp(cfg, sprintf('       cohort.inputs.covariates.%s.TIME = []', cinputs{bidx}));
          end
        else
          isgood = false;
          vp(cfg, sprintf('Error: The covariate >%s< has not been defined for this system', cinputs{bidx}));
          vp(cfg, sprintf('       <CV:%s>; times;  [];  units', cinputs{bidx}));
          vp(cfg, sprintf('       <CV:%s>; values; [];  units', cinputs{bidx}));
        end
      end
    else
      isgood = false;
      vp(cfg, sprintf('Error: An covariate was specified for this cohort but'));
      vp(cfg, sprintf('       there are no covariate defined in the system.txt file.'));
      vp(cfg, sprintf('       <CV:CNAME>; times;  [];  units'));
      vp(cfg, sprintf('       <CV:CNAME>; values; [];  units'));
    end
  end

end

%
% Checking the outputs
%
if(isfield(cohort, 'outputs'))
  outputs = fieldnames(cohort.outputs);
  for(opidx=1:length(outputs))

    % Pulling the current output out
    clear output;
    output = getfield(cohort.outputs, outputs{opidx});


    %
    % Checking the dataset specifications
    %
    % cohort.outputs.outputname.obs
    %
    if(isfield(output, 'obs'))
      if(isfield(output.obs, 'time'))
        if(datasetgood)
          if(~isfield(tmpdataset.column.names, output.obs.time))
            isgood = false;
            vp(cfg, sprintf('Error: For the output >%s< the specified observation time ', outputs{opidx}));
            vp(cfg, sprintf('       column >%s< was not found in the dataset', output.obs.time));
          end
        end
      else
         isgood = false;
         vp(cfg, sprintf('Error: For the output >%s<the column for the "time" must be specified', outputs{opidx}));
         vp(cfg, sprintf('       cohort.outputs.%s.obs.time  = ''name''; ', outputs{opidx}));
      end
      if(isfield(output.obs, 'value'))
        if(datasetgood)
          if(~isfield(tmpdataset.column.names, output.obs.value))
            isgood = false;
            vp(cfg, sprintf('Error: For the output >%s< the specified observation value ', outputs{opidx}));
            vp(cfg, sprintf('       column >%s< was not found in the dataset', output.obs.value));
          end
        end
      else
         isgood = false;
         vp(cfg, sprintf('Error: For the output >%s< the column for the "value" must be specified', outputs{opidx}));
         vp(cfg, sprintf('       cohort.outputs.%s.obs.value = ''name''; ', outputs{opidx}));
      end
      if(datasetgood)
        if(isfield(output, 'of'))
          offields = fieldnames(output.of);
          % looping through each field making sure the column exists in the dataset
          for(ofidx = 1:length(offields))
            if(~isfield(tmpdataset.column.names, offields{ofidx}))
              isgood = false;
              vp(cfg, sprintf('Error: The column >%s< in the output filter ', offields{ofidx}));
              vp(cfg, sprintf('       was not found in the data set >%s< ', cohort.dataset));
            end
          end
        end
      end

    else
      isgood = false;
      vp(cfg, sprintf('Error: For the output >%s< no observation information was specified', outputs{opidx}));
      vp(cfg, sprintf('       cohort.outputs.%s.obs.time  = ''name''; ', outputs{opidx}));
      vp(cfg, sprintf('       cohort.outputs.%s.obs.value = ''name''; ', outputs{opidx}));
    end

    %
    % Checking the model specifications
    %
    % cohort.outputs.outputname.model
    %
    if(isfield(output, 'model'))
      % 
      % Time component
      % 
      if(isfield(output.model, 'time'))
        % making sure there are timescales defined for this model
        if(isfield(cfg.options, 'time_scales'))
          % checking to see if the specified timescale exists.
          if(~isfield(cfg.options.time_scales, output.model.time))
            isgood = false;
            vp(cfg, sprintf('Error: For the output >%s< the specified model timescale >%s<', outputs{opidx}, output.model.time));
            vp(cfg, sprintf('       does not appear to have been defined in the system.txt file'));
            vp(cfg, sprintf('       <TS:%s> value ', output.model.time));
          end
        else
          isgood = false;
          vp(cfg, sprintf('Error: The model contains no timescale information >%s<', output.model.time));
          vp(cfg, sprintf('       The system file needs to be modified to include this information'));
          vp(cfg, sprintf('       <TS:%s> value ', output.model.time));
        end
      else
        isgood = false;
        vp(cfg, sprintf('Error: For the output >%s< the model timescale must be specified', outputs{opidx}));
        vp(cfg, sprintf('       cohort.outputs.%s.model.time  = ''name''; ', outputs{opidx}));
      end
      % 
      % Output component
      % 
      if(isfield(output.model, 'value'))
        if(~isfield(cfg.options.mi.outputs, output.model.value))
            isgood = false;
            vp(cfg, sprintf('Error: For the output >%s< the specified model output >%s<', outputs{opidx}, output.model.value));
            vp(cfg, sprintf('       does not appear to have been defined in the system.txt file'));
            vp(cfg, sprintf('       <O> %s = value ', output.model.value));
        end
      else
        isgood = false;
        vp(cfg, sprintf('Error: For the output >%s< the model output must be specified', outputs{opidx}));
        vp(cfg, sprintf('       cohort.outputs.%s.model.output  = ''name''; ', outputs{opidx}));
      end
      % 
      % variance component
      % 
      if(~isfield(output.model, 'variance'))
        isgood = false;
        vp(cfg, sprintf('Error: For the output >%s< the model variance must be specified', outputs{opidx}));
        vp(cfg, sprintf('       cohort.outputs.%s.model.variance = ''PRED^2''; ', outputs{opidx}));
      end
    else
      isgood = false;
      vp(cfg, sprintf('Error: For the output >%s< no model information was specified', outputs{opidx}));
      vp(cfg, sprintf('       cohort.outputs.%s.model.time  = ''name''; ', outputs{opidx}));
      vp(cfg, sprintf('       cohort.outputs.%s.model.value = ''name''; ', outputs{opidx}));
    end

    %
    % Checking the options. 
    %
    if(isfield(output, 'options'))
     defoptnames = fieldnames(defopts);
     opoptnames  = fieldnames(output.options);
     % First we check to see if all of the specified optiosn
     % are valid (e.g. they have default values). 
     for(optidx = 1:length(opoptnames))
       if(isempty(strmatch(opoptnames{optidx}, validopts)))
         vp(cfg, sprintf('Error: For output >%s< the specified option >%s< is invalid', outputs{opidx}, opoptnames{optidx}));
         vp(cfg, sprintf(' This option will be ignored'));
       end
     end

     % Now we loop through the default options and fill in those that were
     % not specified with the default values
     for(optidx = 1:length(defoptnames))
       if(~isfield(output.options, defoptnames{optidx}))
         eval(sprintf('output.options.%s = defopts.%s;', defoptnames{optidx}, defoptnames{optidx}));
       end
     end

    else
      % If they don't exist we assign the defaults
      output.options = defopts;
    end

    % Putting the current output back
    eval(sprintf('cohort.outputs.%s = output;', outputs{opidx}));
  end
else
  isgood = false;
  vp(cfg, 'Error: No outputs were specified');
end



% If everything checks out (dataset exists, columns specified for the
% outputs exists, etc.) If that's the case we extract the data from the datasets
if(isgood)
  % Pulling out the dataset 
  % for the cohort (chds)
  chds  = getfield(cfg.data, cohort.dataset);
  % selecting just the records for this cohort
  chvalues     = nm_select_records(chds, chds.values, cohort.cf);

  choutput_times = [];
  
  outputs = fieldnames(cohort.outputs);
  for(opidx = 1:length(outputs))
    % making sure the output is stored:
    output = getfield(cohort.outputs, outputs{opidx});
  
    % If the cohort needs to be further filtered for this output
    % we apply that here otherwise we just pass all of the cohort 
    % values through:
    if(isfield(output, 'of'))
      opvalues     = nm_select_records(chds, chvalues, output.of);
    else
      opvalues = chvalues;
    end
  
    if(isempty(opvalues))
     vp(cfg, sprintf('Unable to fetch observations:'));
     vp(cfg, sprintf('Cohort: %s',cohort.name));
     vp(cfg, sprintf('Output: %s',outputs{opidx}));
     vp(cfg, sprintf('Check the filters (cf, of), See:'));
     vp(cfg, sprintf('help system_define_cohort'));
     vp(cfg, sprintf('for more information'));
    end
  
    clear tmpop
    tmpop.time    = nm_fc(chds, opvalues, output.obs.time);
    tmpop.obs     = nm_fc(chds, opvalues, output.obs.value);
  
    % If a missing observation flag has been specified we strip those out here
    if(isfield(output.obs, 'missing'))
      [tmpop.time, tmpop.obs] = strip_missing(tmpop.time, tmpop.obs, output.obs.missing);
    end
  
    % converting the times in the dataset to simtimes
    tmpop.simtime = system_ts_to_simtime(cfg, tmpop.time, output.model.time);
  
    % Adding observation times to the smooth output times
    choutput_times = sort(unique([tmpop.simtime; choutput_times]));
    %eval(sprintf('odchunks.%s = tmpop;', outputs{opidx}));
  
    % Placing the data from the data file into the cohort
    eval(sprintf('cohort.outputs.%s.data.time    = tmpop.time;',   outputs{opidx}))
    eval(sprintf('cohort.outputs.%s.data.simtime = tmpop.simtime;',outputs{opidx}))
    eval(sprintf('cohort.outputs.%s.data.obs     = tmpop.obs ;',   outputs{opidx}))
  end

  % storing all of the observation times for the cohort
  eval(sprintf('cohort.observation_simtimes = choutput_times ;'))
end


if( isgood == true)
  % storing the name field
  cohort_name = cohort.name;
  % chopping off the name field
  cohort = rmfield(cohort, 'name');
  % storing the cohort
  eval(sprintf('cfg.cohorts.%s = cohort;',cohort_name));
else
  vp(cfg, 'system_define_cohort()');
  vp(cfg, sprintf('Cohort name: >%s<', cohort_name));
  vp(cfg, 'There was an error and the cohort information was not set.');
end
