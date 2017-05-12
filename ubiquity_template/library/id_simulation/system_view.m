function [] = system_view(cfg, field) 
% function [] = system_view(cfg, field) 
%
% Displays information about the system described by 'cfg'. The second
% argument 'field' is optional and indicates the aspect of the system to
% display. If the field is not specified, then all of the
% system information will be displayed. The following are equivalent:
%
%   system_view(cfg)
%   system_view(cfg, 'all')
%
% To display only information about the parameters of the system the
% following can be used:
%
%   system_view(cfg, 'parameters')
%
% The other values field can take are:
%
%   'all' (default)
%   'parameters'      - system parameters and current parameter set
%   'bolus'           - bolus dosing information
%   'infusion rates'  - infusion rate information
%   'covariates'      - covariates
%   'iiv'             - iiv details
%   'simulation'      - simulation options
%   'solver'          - colver options
%   'stochastic'      - stochastic simulation options
%   'datasets'        - loaded datasets
%
%   'estimation'      - estimation options
%   'cohorts'         - cohort descriptions
%


 if(not(exist('field', 'var')))
   field = 'all';
 end

  %
  % Parameters
  %
  if(strcmp(field, 'all') | strcmp(field, 'parameters'))
      vp(cfg, sprintf('====================='))
      vp(cfg, sprintf('Parameter Information'));
      vp(cfg, sprintf('====================='))
      vp(cfg, sprintf('Parameter set:'));
      vp(cfg, sprintf('Short name:   %s  ',  cfg.estimation.set_name));
      vp(cfg, sprintf('Description:  %s  ',  getfield(getfield(cfg.parameters.sets, cfg.estimation.set_name), 'name')))
      vp(cfg, sprintf(''));
      vp(cfg, sprintf('Default parameters for current set:  '));
      vp(cfg, sprintf('%s |  %s | %s ',     ...
                  pad_string('name',  20),  ...
                  pad_string('value', 10),  ...
                  pad_string('units', 10)));
     vp(cfg, repmat('-', 1,50));
     for pidx =1:length(cfg.parameters.names)
      vp(cfg, sprintf('%s |  %s | %s ', ...
                  pad_string(cfg.parameters.names{pidx},  20),  ...
                  var2string(cfg.parameters.values(pidx), 10),  ...
                  pad_string(cfg.parameters.units{pidx},  10)))
     end
     vp(cfg, repmat('-', 1,50));
     vp(cfg, sprintf(''));
  end

  %
  % Bolus dosing
  %
  if(strcmp(field, 'all') | strcmp(field, 'bolus'))
    vp(cfg, sprintf(''))
    vp(cfg, sprintf('===================='))
    vp(cfg, sprintf('Bolus dosing details'))
    vp(cfg, sprintf('===================='))
    if(isfield(cfg.options.inputs, 'bolus'))
      vp(cfg, sprintf('%s |  %s | %s | %s', ...
                   pad_string('field',   20), ...
                   pad_string('values',  20), ...
                   pad_string('scaling', 30), ...
                   pad_string('units',   10)))
      vp(cfg, repmat('-', 1,90));
      vp(cfg, sprintf('%s |  %s | %s | %s',  ...
                   pad_string('times', 20), ...
                   pad_string(num2str(cfg.options.inputs.bolus.times.values), 20), ...
                   pad_string(cfg.options.inputs.bolus.times.scale, 30), ...
                   pad_string(cfg.options.inputs.bolus.times.units, 10)))
      species = fieldnames(cfg.options.inputs.bolus.species);
      for(spidx = 1:length(species))
       vp(cfg, sprintf('%s |  %s | %s | %s',...
                    pad_string(species{spidx}, 20), ...
                    pad_string(num2str(getfield(getfield(cfg.options.inputs.bolus.species, species{spidx}), 'values')), 20), ...
                    pad_string(getfield(getfield(cfg.options.inputs.bolus.species, species{spidx}), 'scale'), 30), ...
                    pad_string(getfield(getfield(cfg.options.inputs.bolus.species, species{spidx}), 'units'), 10)))
      end
      vp(cfg, repmat('-', 1,90));
    else
      vp(cfg, sprintf('No bolus information found'));
    end
    vp(cfg, sprintf(''));
  end

  %
  % Infusion rates
  %
  if(strcmp(field, 'all') | strcmp(field, 'infusion rates'))
    vp(cfg, sprintf(''))
    vp(cfg, sprintf('====================='))
    vp(cfg, sprintf('Infusion rate details'))
    vp(cfg, sprintf('====================='))
    if(isfield(cfg.options.inputs, 'infusion_rates'))
      vp(cfg, sprintf(' %s | %s | %s | %s | %s', ...
                   pad_string('Rate ',   10), ...
                   pad_string('field',   10), ...
                   pad_string('values',  10), ...
                   pad_string('scaling', 10), ...
                   pad_string('units',   10)))
      vp(cfg, repmat('-', 1,65));
      for(ridx=1:length(cfg.options.inputs.infusion_rate_names))
        current_rate = getfield(cfg.options.inputs.infusion_rates, cfg.options.inputs.infusion_rate_names{ridx});
        vp(cfg, sprintf(' %s | %s | %s | %s | %s', ...
                     pad_string(cfg.options.inputs.infusion_rate_names{ridx}, 10), ...
                     pad_string('time', 10), ...
                     pad_string(num2str(current_rate.times.values), 10), ...
                     pad_string(        current_rate.times.scale,   10), ...
                     pad_string(        current_rate.times.units,   10)))
        vp(cfg, sprintf(' %s | %s | %s | %s | %s', ...
                     pad_string('', 10), ...
                     pad_string('levels', 10), ...
                     pad_string(num2str(current_rate.levels.values), 10), ...
                     pad_string(        current_rate.levels.scale,   10), ...
                     pad_string(        current_rate.levels.units,   10)))
      end
      vp(cfg, repmat('-', 1,65));
    else
      vp(cfg, sprintf('No infusion rate information found'));
    end
    vp(cfg, sprintf(''));
  end

  %
  % Covariates
  %
  if(strcmp(field, 'all') | strcmp(field, 'covariates'))
    vp(cfg, sprintf('================='))
    vp(cfg, sprintf('Covariate details'))
    vp(cfg, sprintf('================='))
    if(isfield(cfg.options.inputs, 'covariates'))
      vp(cfg, sprintf(' %s | %s | %s | %s', ...
                   pad_string('Rate ',   10), ...
                   pad_string('field',   10), ...
                   pad_string('values',  10), ...
                   pad_string('units',   10)))
      vp(cfg, repmat('-', 1,52));
      for(cidx=1:length(cfg.options.inputs.covariate_names))
        current_cov = getfield(cfg.options.inputs.covariates, cfg.options.inputs.covariate_names{cidx});
        vp(cfg, sprintf(' %s | %s | %s | %s', ...
                     pad_string(cfg.options.inputs.covariate_names{cidx}, 10), ...
                     pad_string('time', 10), ...
                     pad_string(num2str(current_cov.times.values), 10), ...
                     pad_string(        current_cov.times.units,   10)))
        vp(cfg, sprintf(' %s | %s | %s | %s', ...
                     pad_string('', 10), ...
                     pad_string('values', 10), ...
                     pad_string(num2str(current_cov.values.values), 10), ...
                     pad_string(        current_cov.values.units,   10)))
      end
      vp(cfg, repmat('-', 1,52));
    else
      vp(cfg, sprintf('No covariate information found'));
    end
    vp(cfg, sprintf(''));
  end

  %
  % IIV details
  %
  if(strcmp(field, 'all') | strcmp(field, 'iiv'))
    vp(cfg, sprintf('==========='))
    vp(cfg, sprintf('IIV details'))
    vp(cfg, sprintf('==========='))
    if(isfield(cfg.iiv, 'iivs'))
       
    vp(cfg, sprintf('IIV/Parameter set'))
    vp(cfg, sprintf('Short name:   %s', cfg.iiv.current_set))
    
    vp(cfg, sprintf('Variance/covariance matrix'))
      iivs = fieldnames(cfg.iiv.iivs);

      colwidth = max(cellfun('length', fieldnames(cfg.iiv.iivs))) +2;

      % creating the headers
      row = pad_string('', colwidth);
      for iidx =1:length(iivs)
        row = sprintf('%s%s', row,pad_string(iivs{iidx}, colwidth));
      end
      vp(cfg, row);
      for ridx =1:length(iivs)
        row = pad_string(iivs{ridx}, colwidth);
        for cidx =1:length(iivs)
          row = sprintf('%s%s', row,var2string(cfg.iiv.values(ridx,cidx), colwidth));
        end
        vp(cfg, row);
      end
        
    vp(cfg, sprintf(''))
      vp(cfg, sprintf('On parameters'))
      pnames = fieldnames(cfg.iiv.parameters);
      for pidx =1:length(pnames)
      
      vp(cfg, sprintf('%s, %s(%s)', ...
           pnames{pidx}, ...
           getfield(getfield(cfg.iiv.parameters, pnames{pidx}), 'iiv_name'), ...
           getfield(getfield(cfg.iiv.parameters, pnames{pidx}), 'distribution')));
      end
    else
      vp(cfg, sprintf('No IIV information found'));
    end
  end

  %
  % Simulation 
  %
  if(strcmp(field, 'all') | strcmp(field, 'simulation'))
    vp(cfg, sprintf(''))
    vp(cfg, sprintf('=================='))
    vp(cfg, sprintf('Simulation details'))
    vp(cfg, sprintf('=================='))

    % now we get the remaining simulation options
    if(isfield(cfg.options.simulation_options, 'integrate_with'))
      vp(cfg, sprintf('integrate_with                   %s', cfg.options.simulation_options.integrate_with));
    end
    if(isfield(cfg.options.simulation_options, 'include_important_output_times'))
      vp(cfg, sprintf('include_important_output_times   %s', cfg.options.simulation_options.include_important_output_times));
    end
    if(isfield(cfg.options.simulation_options, 'output_times'))
      vp(cfg, sprintf('output_times                     %s', var2string_gen(cfg.options.simulation_options.output_times)));
    end
  end

  %
  % Solver
  %

  if(strcmp(field, 'all') | strcmp(field, 'solver'))
    vp(cfg, sprintf(''))
    vp(cfg, sprintf('=============='))
    vp(cfg, sprintf('Solver details'))
    vp(cfg, sprintf('=============='))

    % pulling out the solver options
    if(isfield(cfg.options.simulation_options, 'default_simopts'))
      dsopts = fieldnames(cfg.options.simulation_options.default_simopts);
      if(~isempty(dsopts))
       vp(cfg, repmat('-', 1,50));
       for(dsoidx=1:length(dsopts))
         
         mylen = 30 - length(dsopts{dsoidx});
         vp(cfg, sprintf('%s%s   %s', ...
                 dsopts{dsoidx}, ...
                 repmat(' ', 1, mylen), ...
                 var2string_gen(getfield(cfg.options.simulation_options.default_simopts, dsopts{dsoidx}))))
       end
       vp(cfg, repmat('-', 1,50));
      end
    else
      vp(cfg, sprintf('No solver options found'));
    end
    vp(cfg, sprintf(''))
  end

  %
  % Stochastic 
  %

  if(strcmp(field, 'all') | strcmp(field, 'stochastic'))
    vp(cfg, sprintf('============================='))
    vp(cfg, sprintf('Stochastic simulation details'))
    vp(cfg, sprintf('============================='))
    if(isfield(cfg.options, 'stochastic'))
      sopts = fieldnames(cfg.options.stochastic);
      if(~isempty(sopts))
       for(soidx=1:length(sopts))
         mylen = 30 - length(sopts{soidx});
         vp(cfg, sprintf('%s%s   %s', ...
                 sopts{soidx}, ...
                 repmat(' ', 1, mylen), ...
                 var2string_gen(getfield(cfg.options.stochastic, sopts{soidx}))))
       end
      end

    else
      vp(cfg, sprintf('No stochastic options found'));
    end
    vp(cfg, sprintf(''))
  end


  %
  % dataset
  %
  if(strcmp(field, 'all') | strcmp(field, 'datasets'))
    vp(cfg, sprintf('==============='))
    vp(cfg, sprintf('Dataset details'))
    vp(cfg, sprintf('==============='))

    datasets_found = false;


    if(isfield(cfg, 'data'))
      dsets = fieldnames(cfg.data);
      for(dsidx=1:length(dsets))
        dset = getfield(cfg.data, dsets{dsidx});
        if(~isempty(fieldnames(dset)))
          % we've found a dataset so there's at least one :)

          vp(cfg, repmat('-', 1,15));
          vp(cfg, sprintf('Name:      %s', dsets{dsidx}));
          vp(cfg, sprintf('Data File: %s', dset.data_file.name));
          if(isfield(dset.data_file, 'sheet'))
            vp(cfg, sprintf('Sheet:     %s', dset.data_file.sheet));
          end
          vp(cfg, sprintf('Columns:   %s', [dset.header{1} sprintf(', %s', dset.header{2:end})]));
          vp(cfg, sprintf('Rows:      %d', length(dset.values(:,1))));


          datasets_found = true;
        end
      end
    end

    if(~datasets_found)
      vp(cfg, sprintf('No datasets found'));
    end
    vp(cfg, sprintf(''))
  end
  %{
  %}


  %
  % estimation
  %
  if(strcmp(field, 'all') | strcmp(field, 'estimation'))
    vp(cfg, sprintf('=================='))
    vp(cfg, sprintf('Estimation details'))
    vp(cfg, sprintf('=================='))
    vp(cfg, sprintf('Parameter set:        %s', cfg.estimation.set_name));
    vp(cfg, sprintf('Parameters estimated: %s', [cfg.estimation.parameters.names{1} sprintf(', %s', cfg.estimation.parameters.names{2:end})]));
    vp(cfg, sprintf('objective_type        %s', cfg.estimation.objective_type));
    vp(cfg, sprintf('effort                %d', cfg.estimation.effort));
    vp(cfg, sprintf('observation_function  %s', cfg.estimation.observation_function));

    vp(cfg, '');
    vp(cfg, sprintf('Monitor                 '));
    vp(cfg, repmat('-', 1,50));
    vp(cfg, sprintf('  status_function     %s', cfg.estimation.monitor.status_function));
    vp(cfg, sprintf('  exit_when_stable    %s', cfg.estimation.monitor.exit_when_stable));
    vp(cfg, sprintf('  iteration_history   %d', cfg.estimation.monitor.iteration_history));
    vp(cfg, sprintf('  slope_tolerance     %s', var2string(cfg.estimation.monitor.slope_tolerance,1)));
    vp(cfg, repmat('-', 1,50));
    eopts = fieldnames(cfg.estimation.options);
    if(~isempty(eopts))
      vp(cfg, '');
      vp(cfg, sprintf('Optimization Options    '));
      vp(cfg, repmat('-', 1,50));
      for(eidx=1:length(eopts))
        
        if(~isempty(getfield(cfg.estimation.options, eopts{eidx})))
        mylen = 17 - length(eopts{eidx});
        vp(cfg, sprintf('   %s%s   %s', ...
                eopts{eidx}, ...
                repmat(' ', 1, mylen), ...
                var2string_gen(getfield(cfg.estimation.options, eopts{eidx}))))
        end
      end
      vp(cfg, repmat('-', 1,50));
    end

    vp(cfg, sprintf(''))

  end

  %
  % cohorts  
  %
  if(strcmp(field, 'all') | strcmp(field, 'cohorts'))
    vp(cfg, sprintf('=============='))
    vp(cfg, sprintf('Cohort Details'));
    vp(cfg, sprintf('=============='))

    if(~isempty(fieldnames(cfg.cohorts)))
      cohorts = fieldnames(cfg.cohorts);

      for(chidx = 1:length(cohorts))
        cohort = getfield(cfg.cohorts, cohorts{chidx});
        vp(cfg, sprintf('Cohort:  %s', cohorts{chidx}));
        vp(cfg, repmat('-', 1,50));

        %
        % dataset
        %
        if(isfield(cohort, 'dataset'))
          vp(cfg, sprintf('dataset            %s', cohort.dataset));
          vp(cfg, sprintf(''))
        end

        vp(cfg, sprintf('Cohort options (options)'));
        if(length(fieldnames(cohort.options)) > 0)
          vp(cfg, sprintf(''))
        else
          vp(cfg, sprintf('  None set.'));
        end
        vp(cfg, sprintf(''))

        %
        % cohort filter
        %
        vp(cfg, sprintf('Cohort filter (cf)'));
        if(isfield(cohort, 'cf'))
          cfcs = fieldnames(cohort.cf);
          for(cfidx = 1:length(cfcs))
            vp(cfg, sprintf('%s %s', pad_string(cfcs{cfidx}, 10), ...
                    mat2str(getfield(cohort.cf, cfcs{cfidx}))));
          end
        else
          vp(cfg, sprintf('  None set.'));
        end
        vp(cfg, sprintf(''))

        %
        % parameter overwritten at the cohort level
        %
        vp(cfg, sprintf('Cohort-specific parameters (cp)'));
        if(isfield(cohort, 'cp'))
          cps = fieldnames(cohort.cp);
          for(cpidx = 1:length(cps))
            vp(cfg, sprintf('%s %s', pad_string(cps{cpidx}, 10), ...
               var2string_gen(getfield(cohort.cp, cps{cpidx}))));
          end
        else
          vp(cfg, sprintf('  None set.'));
        end
        vp(cfg, sprintf(''))


        %
        % Model Inputs
        %
        vp(cfg, sprintf('Inputs'));
        if(isfield(cohort, 'inputs'))
          % Bolus inputs
          if(isfield(cohort.inputs, 'bolus'))
            inpts = fieldnames(cohort.inputs.bolus);
            for(inptidx =1:length(inpts))
              vp(cfg, sprintf('  bolus: %s', inpts{inptidx}));
              input = getfield(cohort.inputs.bolus, inpts{inptidx});
              vp(cfg, sprintf('    TIME %s',  mat2str(input.TIME)));
              vp(cfg, sprintf('    AMT  %s',  mat2str(input.AMT)));
            end
            vp(cfg, sprintf(''))
          end
          % Infusion Rates
          if(isfield(cohort.inputs, 'infusion_rates'))
            inpts = fieldnames(cohort.inputs.infusion_rates);
            for(inptidx =1:length(inpts))
              vp(cfg, sprintf('  infusion_rates: %s', inpts{inptidx}));
              input = getfield(cohort.inputs.infusion_rates, inpts{inptidx});
              vp(cfg, sprintf('    TIME %s',  mat2str(input.TIME)));
              vp(cfg, sprintf('    AMT  %s',  mat2str(input.AMT)));
            end
            vp(cfg, sprintf(''))
          end
          % Covariates     
          if(isfield(cohort.inputs, 'covariates'))
            inpts = fieldnames(cohort.inputs.covariates);
            for(inptidx =1:length(inpts))
              vp(cfg, sprintf('  covariates: %s', inpts{inptidx}));
              input = getfield(cohort.inputs.covariates, inpts{inptidx});
              vp(cfg, sprintf('    TIME %s',  mat2str(input.TIME)));
              vp(cfg, sprintf('    AMT  %s',  mat2str(input.AMT)));
            end
            vp(cfg, sprintf(''))
          end
        else
          vp(cfg, sprintf('  No inputs have been defined.'));
          vp(cfg, sprintf(''))
        end

        %
        % Model Outputs 
        %
        vp(cfg, sprintf('Outputs'));
        if(isfield(cohort, 'outputs'))
          otpts = fieldnames(cohort.outputs);
          for(otptidx =1:length(otpts))
            otput = getfield(cohort.outputs,  otpts{otptidx});
            vp(cfg, sprintf('  -%s-', repmat('-', 1, length(otpts{otptidx}))));
            vp(cfg, sprintf('   %s ', otpts{otptidx}));
            vp(cfg, sprintf('  -%s-', repmat('-', 1, length(otpts{otptidx}))));
            vp(cfg, sprintf('               %s | %s ', pad_string('model', 12), pad_string('obs', 12)));
            vp(cfg, sprintf('               %s', repmat('-', 1,27)));
            vp(cfg, sprintf('     time      %s | %s ', pad_string(otput.model.time,  12), pad_string(otput.obs.time,  12)));
            vp(cfg, sprintf('     output    %s | %s ', pad_string(otput.model.value, 12), pad_string(otput.obs.value, 12)));
            vp(cfg, sprintf('               %s', repmat('-', 1,27)));
            vp(cfg, sprintf('     variance  %s ',otput.model.variance))
            if(isfield(otput.obs, 'missing'))
            vp(cfg, sprintf('     missing   %d ',otput.obs.missing))
            end
            vp(cfg, sprintf(''))

            %cohort.outputs.pk.of.cmt 
            
            vp(cfg, sprintf('     Output filter (of)'));
            if(isfield(otput, 'of'))
              ofcs = fieldnames(otput.of);
              for(ofidx = 1:length(ofcs))
                vp(cfg, sprintf('%s %s', pad_string(ofcs{ofidx}, 10), ...
                        mat2str(getfield(otput.of, ofcs{ofidx}))));
              end
            else
              vp(cfg, sprintf('       None set.'));
            end
          end
          vp(cfg, sprintf(''))
        else
          vp(cfg, sprintf('No outputs been defined'));
        end
        vp(cfg, repmat('-', 1,50));
        vp(cfg, sprintf(''))

      end


    else
      vp(cfg, sprintf('No cohorts have been defined'));
    end
      
  end


