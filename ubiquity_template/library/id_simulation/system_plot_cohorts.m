function [] = system_plot_cohorts(erp,  plot_opts, cfg);
% function [] = system_plot_cohorts(erp,  plot_opts, cfg);
%
% erp - Output from system_simulate_estimation_results.
% plot_opts - Data structure used to control plotting. 
%       The details are described below for an output 
%       named ONAME.
%
% Controlling the layout of the data and predictions
% --------------------------------------------------
% plot_opts.outputs.ONAME.panels   = 'no';
%
% With the value 'no', all of the cohorts for a given output will be plotted
% on the same axis. If this field is set to 'yes', then each cohort will be
% plotted on a separate axis. To control the number of rows and columns per
% figure set the following in conjunction with panels = 'yes':
%
% plot_opts.outputs.ONAME.rows     = 2;
% plot_opts.outputs.ONAME.cols     = 3;
%
% If there are more cohorts than panels on a single figure then the overflow
% will be pushed onto additional figures.
%
% Figure dimensions
% -----------------
% plot_opts.outputs.ONAME.Position = [];
%
% By default, Matlab will choose the position and dimensions of the figure. If
% you want to alter these dimensions simply resize the figure with the mouse
% until it appears the way you want it and use the following command to get
% the position:
%
% get(gcf, 'Position')
%
% Then take the four numbers that are returned and place them in the [] above.
% 
% Controlling the axes
% --------------------
% The following can be used 
%                                                       % Default
% plot_opts.outputs.ONAME.yscale   = 'log';             % 'linear'
% plot_opts.outputs.ONAME.ylim     = [1,    100];       %   auto
% plot_opts.outputs.ONAME.xlim     = [0,    100];       %   auto
% plot_opts.outputs.ONAME.ylabel   = 'OUTPUT (units)';  % model output
% plot_opts.outputs.ONAME.xlabel   = 'Time (units)';    % model TS
%
% Determining how data is presented
% ---------------------------------
% plot_opts.outputs.ONAME.datatype = 'all'  
%
% This field, when set to 'all' will simply plot all of the data for the given
% cohort/output. It the data are collected at the same times, it's possible to
% set this field to 'mean' and plot the mean and standard deviation instead. 
%
% When datatype is set to 'mean' it may be desirable to change how error is
% represented: 
% plot_opts.outputs.ONAME.spread   = 'stdev';  % standard deviation (default)
% plot_opts.outputs.ONAME.spread   = 'stderr'; % standard error 
%
%
% Writing figures
% ---------------
% plot_opts.save_figs = analysis_name;
% 
% If this field is specified it should be a string. That will write the
% figures generated to the output directory with this string as a prefix. If
% set to the analysis_name it will be grouped with the other information
% associated with the parameter estimation.
%
% The following can be created for each output. The defaults are listed here
% or as comments to the right
% plot_opts.outputs.ONAME.panels   = 'no';
% plot_opts.outputs.ONAME.rows     = 1;
% plot_opts.outputs.ONAME.cols     = 2;
% plot_opts.outputs.ONAME.datatype = 'mean'; % 'all' for data
% plot_opts.outputs.ONAME.Position = [];


def.panels   = 'no';
def.rows     = 1;
def.cols     = 2;
def.yscale   = 'linear';
def.datatype = 'all';
def.spread   = 'stderr';

% by default we don't save figures
save_figs  = 0;

figctr = 1;
spctr  = 1;

cohorts = fieldnames(erp.cohorts); 

for opidx = 1:length(erp.meta.outputs)
  figure(figctr); hold on;
  
  clear myylim;
  clear myxlim;
  
  % combining output options and default optoins
  op_plot_opts = fetch_output_plotopts(def, plot_opts, erp.meta.outputs{opidx});
  spr      = op_plot_opts.rows;
  spc      = op_plot_opts.cols;
  panels   = op_plot_opts.panels;
  yscale   = op_plot_opts.yscale;
  datatype = op_plot_opts.datatype;
  spread   = op_plot_opts.spread  ;
  if(isfield(op_plot_opts, 'ylim'))
     myylim   = op_plot_opts.ylim;
  end
  if(isfield(op_plot_opts, 'xlim'))
     myxlim   = op_plot_opts.xlim;
  end
  if(isfield(plot_opts, 'save_figs'))
    save_figs  = 1;
    fig_file   =  plot_opts.save_figs;
  end
  
  
  % Plotting timecourse
  spctr = 1;
  lh = [];
  ll = {};
  for chidx = 1:length(cohorts)
    cohort = getfield(erp.cohorts, cohorts{chidx});
  
    % Checking to see if the current cohort has the current output
    opexists = isfield(cohort, erp.meta.outputs{opidx});
  
    if(opexists)
      output = getfield(cohort, erp.meta.outputs{opidx});
      eval(sprintf('choptions = cfg.cohorts.%s.outputs.%s.options;', cohorts{chidx}, erp.meta.outputs{opidx}));
      eval(sprintf('chmodel   = cfg.cohorts.%s.outputs.%s.model  ;', cohorts{chidx}, erp.meta.outputs{opidx}));
    end
      
    
    % This will split out the cohorts on subplots,
    % otherwise they will be placed on a single plot
    if(strcmp(panels, 'yes'))
      spmax = spr*spc;
      subplot(spr, spc, spctr);
      hold on;
      lh = [];
      ll = {};
    end
    
    if(opexists)
      if(strcmp(datatype, 'all'))
         plot(output.od.time,        output.od.obs,         choptions.marker_shape, 'markerfacecolor', choptions.marker_color, 'color', choptions.marker_color);
      elseif(strcmp(datatype, 'mean'))
         [opsum] = find_data_mean(output.od.time, output.od.obs);


         opsum.ub = opsum.obsstd;
         opsum.lb = opsum.obsstd;
         %errorbar(opsum.times, opsum.obsmean, opsum.lb, opsum.ub, choptions.marker_shape, 'markerfacecolor', choptions.marker_color, 'color', choptions.marker_color);


         if(strcmp(spread, 'stdev'))
           ebh = errorbar(opsum.times, opsum.obsmean, opsum.obsstd, choptions.marker_shape, 'markerfacecolor', choptions.marker_color, 'color', choptions.marker_color);

         else
           ebh = errorbar(opsum.times, opsum.obsmean, opsum.obsse , choptions.marker_shape, 'markerfacecolor', choptions.marker_color, 'color', choptions.marker_color);
         end

         % A bit of a hack for the situation where the lb of 
         % the error bar is negative:
         if(strcmp(yscale, 'log'))
           % taken from:
           % https://www.mathworks.com/matlabcentral/answers/241261-how-do-i-get-working-y-axis-errorbars-using-a-log-scaleI% 
           ebh.LData = ebh.YData - max(10000*eps,ebh.YData-ebh.LData);
         end

      end
      tmplh = plot(output.od.time_smooth, output.od.pred_smooth, choptions.marker_line,  'markerfacecolor', choptions.marker_color, 'color', choptions.marker_color);
      lh(end+1) = tmplh;
      if(isfield(getfield(cfg.cohorts, cohorts{chidx}), 'label'))
        ll(end+1) = {getfield(getfield(cfg.cohorts, cohorts{chidx}), 'label')};
      else
        ll(end+1) = cohorts(chidx);
      end
    end
  
    if(exist('myylim', 'var'))
      ylim(myylim);
    end
    if(exist('myxlim', 'var'))
      xlim(myxlim);
    end
    set(gca, 'yscale',yscale);
    prepare_figure('present');
  
    if(isfield(op_plot_opts, 'xlabel'))
      xlabel(op_plot_opts.xlabel, 'interpreter', 'none');
    else
      xlabel(chmodel.time, 'interpreter', 'none');
    end
  
    if(isfield(op_plot_opts, 'ylabel'))
      ylabel(op_plot_opts.ylabel, 'interpreter', 'none');
    else
      ylabel(chmodel.value, 'interpreter', 'none');
    end
  
    if(isfield(op_plot_opts, 'Position'))
      set(gcf, 'Position', op_plot_opts.Position);
    end
    
      
    if(strcmp(panels, 'yes'))
      legend(lh, ll, 'location', 'best', 'interpreter', 'none');
  
      if((spctr < spmax) | (spmax == length(cohorts)));
        % The figure isn't full yet.
        spctr = spctr+1;
      else
        % The figure is full, move on to the next one
        spctr = 1;

        % if save_figs is set we write the current figure 
        if(save_figs)
          dump_figure(sprintf('output%s%s_timecourse_%s_%d', filesep, fig_file, erp.meta.outputs{opidx}, figctr))
        end

        figctr = figctr+1;
        figure(figctr);
      end
    end
  end
  
  if(strcmp(panels, 'no'))
    legend(lh, ll, 'location', 'best', 'interpreter', 'none');
  end
  if(save_figs)
    dump_figure(sprintf('output%s%s_timecourse_%s_%d', filesep, fig_file, erp.meta.outputs{opidx}, figctr))
  end
figctr = figctr+1;
end

flat_headers = erp.flat(1,:);
% getting the headers from flat
h_output_idx = strcmp(flat_headers,  'output');
h_cohort_idx = strcmp(flat_headers,  'cohort');
h_obs_idx    = strcmp(flat_headers,  'obs'   );
h_pred_idx   = strcmp(flat_headers,  'pred'  );
h_type_idx   = strcmp(flat_headers,  'type'  );

flat_data    = erp.flat(2:end,:);

% pulling out the records
flat_records = flat_data(strcmp(flat_data(:,h_type_idx), 'record'),:);


%
% Plotting obs vs pred
%
for opidx = 1:length(erp.meta.outputs)
  oname = erp.meta.outputs{opidx};
  % getting the rows where the output has data
  r_output_idx = strcmp(flat_records(:,h_output_idx),oname);

  op_plot_opts = fetch_output_plotopts(def, plot_opts, oname);


  % If there are outputs in this dataset ...
  if(sum(r_output_idx) > 0)
    % we pull all of the cohorts then iterate through them.
    cohorts = unique(flat_records(r_output_idx, h_cohort_idx));

    % all of the data for the output:
    output_data = flat_records(r_output_idx,:);
    figure(figctr); hold on;
    lh = [];
    ll = {};
    ad = [];
    for chidx = 1:length(cohorts)
      cname = cohorts{chidx};
      eval(sprintf('choptions = cfg.cohorts.%s.outputs.%s.options;', cname, oname));
      r_cohort_idx = strcmp(output_data(:,h_cohort_idx), cname);
      obs   = cell2mat(output_data(r_cohort_idx, h_obs_idx));
      pred  = cell2mat(output_data(r_cohort_idx, h_pred_idx));
      ad = [ad; pred];
      ad = [ad; obs];
      tmplh = plot(pred, obs, choptions.marker_shape, 'markerfacecolor', choptions.marker_color, 'color', choptions.marker_color);
      lh(end+1) = tmplh;
      if(isfield(getfield(cfg.cohorts, cohorts{chidx}), 'label'))
        ll(end+1) = {getfield(getfield(cfg.cohorts, cohorts{chidx}), 'label')};
      else
        ll(end+1) = cohorts(chidx);
      end
    end
    bounds = [min(ad), max(ad)];
    ylim(bounds);
    xlim(bounds);
    set(gca, 'yscale', op_plot_opts.yscale);
    set(gca, 'xscale', op_plot_opts.yscale);
    plot(bounds, bounds, '--k');
    prepare_figure('present');
    title(oname, 'interpreter', 'none')
    xlabel('Predicted')
    ylabel('Observed')
    legend(lh, ll, 'location', 'best', 'interpreter', 'none');
    if(save_figs)
      dump_figure(sprintf('output%s%s_obs_pred_%s_%d', filesep, fig_file, erp.meta.outputs{opidx}, figctr))
    end
    figctr = figctr+1;
  end
end

function [dsum] = find_data_mean(times, values)


all_times = unique(times);

dsum.times   = [];
dsum.obsmean = [];
dsum.obsstd  = [];
dsum.obsse   = [];

for tidx=1:length(all_times)
  dsum.times(end+1)   = all_times(tidx);
  dsum.obsmean(end+1) = mean(values(times == all_times(tidx)));
  dsum.obsstd(end+1)  =  std(values(times == all_times(tidx)));
  dsum.obsse (end+1)  =  std(values(times == all_times(tidx)))/sqrt(sum(times == all_times(tidx)));
end


function [opo] = fetch_output_plotopts(def, plot_opts, output)


% set the defaults
opo = def;

if(isfield(plot_opts, 'outputs'))
if(isfield(plot_opts.outputs, output))
 user_opo = getfield(plot_opts.outputs, output);
 opnames = fieldnames(user_opo);
 for(opidx = 1:length(opnames))
   eval(sprintf('opo.%s = user_opo.%s;', opnames{opidx}, opnames{opidx}));
 end
end
end 
