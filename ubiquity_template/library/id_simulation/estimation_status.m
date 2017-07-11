function [stop]=estimation_status(pest, state,cfg)
% function [stop]=estimation_status(pest, state, cfg)
% 
%  This function is called at each iteration of an optimization to plot the
%  progress of the optimization and to evaluate stopping criteria for the
%  optimization. It's possible to copy and modify this function. 
% 
% Inputs:
% pest  - parameters estimates at the current iteration
% state - data structure containing information about the current iteration
%        (e.g. objective function value, iteration number, etc.).
%
% cfg   - system variable sent to estimate_parameters.
%
% Outputs
%  stop - false if the optimization should continue and 'true' if it should
%  terminate


[pest, errorflag, objmult]=bound_parameters(pest , cfg.estimation.parameters.lower_bound, cfg.estimation.parameters.upper_bound);

stop=false;

% pulling out the initial guess for the 
% parameters and the bounds
pinit = cfg.estimation.parameters.guess;
if(isrow(pinit))
 pinit = pinit';
end
plb   = cfg.estimation.parameters.lower_bound;
pub   = cfg.estimation.parameters.upper_bound;

% calculating the percent shift needed to 
% reach the upper bounds
plb_percent  = ( plb./pinit - 1).*100;
pub_percent  = ( pub./pinit - 1).*100;
pest_percent = (pest./pinit - 1).*100;
ncol = 5;
nrow = ceil(length(pest)/ncol) + 1;


% Default termination options
iteration_history= 100;
slope_tolerance  = .001;
exit_when_stable = 'no';

% Overwriting stopping criteria with user defined values
if(isfield(cfg.estimation.monitor,                    'criteria'))
if(isfield(cfg.estimation.monitor.criteria,           'stability'))
if(isfield(cfg.estimation.monitor.criteria.stability, 'exit_when_stable'))
  exit_when_stable  = cfg.estimation.monitor.criteria.stability.exit_when_stable;
end
if(isfield(cfg.estimation.monitor.criteria.stability, 'slope_tolerance'))
  slope_tolerance   = cfg.estimation.monitor.criteria.stability.slope_tolerance;
end
if(isfield(cfg.estimation.monitor.criteria.stability, 'iteration_history'))
  iteration_history = cfg.estimation.monitor.criteria.stability.iteration_history;
end
end
end


if(state.iteration == 1)
   vp(cfg,'-----------------------------------');
   vp(cfg,'       estimation_status.m         ');
   vp(cfg,'       -------------------         ');
   vp(cfg,' Red series (parameters):          ');
   vp(cfg,'      percent change from initial  ');
   vp(cfg,'      guess.                       ');
   vp(cfg,'                                   ');
   vp(cfg,' Red series (objective):           ');
   vp(cfg,'      progress in minimizing the   ');
   vp(cfg,'      objective function           ');
   vp(cfg,'                                   ');
   vp(cfg,' Blue dashed lines: upper or lower ');
   vp(cfg,'      bound for parameters         ');
   vp(cfg,'-----------------------------------');

end


% If we haven't run this routine before we create a figure to hold all the
% goodness. This is determined by looking to see if a figure named 
% 'Estimation Status' exists. Note if you customize this function it's
% important that you retain the name 'Estimation Status' for the figure
% since it is that name that the estimate_parameters function will look for.
if(isempty(findobj('type', 'figure', 'name', 'Estimation Status')))
  % creating the new figure
  figure('name', 'Estimation Status');
  % Making it take up most of the screen. The user can make it smaller if they want:
  set(figure(findobj('type', 'figure', 'name', 'Estimation Status')), ...
      'Position',get(0,'ScreenSize') - [-30 -40 120 200]);
  %usuptitle('Percent parameter change from (initial guess) and OBJ progress ');

  % Making the subplots for the different parameters. Just setting up the
  % axes, customizing the tick labels, titling each axis, etc.
  for pest_idx = 1:length(pest)
    subplot(nrow, ncol, pest_idx); cla; hold on; 
    title(sprintf('%s (%s)', ...
          cfg.estimation.parameters.names{pest_idx}, ...
           num2str_pretty(pinit(pest_idx))), ...  
          'interpreter', 'none');
    grid on;
    prepare_figure('present')

    % Removing the xticklabels for all but the 
    % bottom parameter subplots
    if(pest_idx <= length(pest) - ncol)
      set(gca, 'xticklabel', []);
    end
  end
  % Creating the axis to track the progress of the objective function:
  subplot(nrow, ncol, [(nrow*ncol-ncol+1) nrow*ncol]); cla; hold on;
  grid on;
  prepare_figure('present')
  ylabel('Objective');
  xlabel('Estimation Iteration');
else
  % If that figure already exists, we make it the current figure:
  figure(findobj('type', 'figure', 'name', 'Estimation Status'));
end


% Making plots for the different parameters being 
% estimated
pest_data = [];

% default value for xlim overwritten dynamically below
myxlim = [0 1];
for pest_idx = 1:length(pest)
  sph = subplot(nrow, ncol, pest_idx); 

  % this plots the first point and subsequent points are then just added to
  % the series on the plot using the child handles 
  if((state.iteration == 0) &  (length(get(sph, 'children')) == 0))
    plot(state.iteration, pest_percent(pest_idx), 'r.-', 'DisplayName', 'pest');
  else


    child_pest  = findobj(get(sph, 'children'), 'DisplayName', 'pest');
    % updating the parameter series
    set(child_pest, 'XData', [state.iteration                         get(child_pest, 'XData')]) ;
    set(child_pest, 'YData', [pest_percent(pest_idx) get(child_pest, 'YData')]) ;
    pest_data = [pest_data; fliplr(get(child_pest, 'YData'))];

   %% updating the parameter bounds series
    if(not(pub_percent(pest_idx) == inf))
      if(length(findobj(get(sph, 'children'), 'DisplayName', 'pctub')) > 0)
        child_pctub = findobj(get(sph, 'children'), 'DisplayName', 'pctub');
        set(child_pctub, 'XData', [state.iteration                         get(child_pctub, 'XData')]) ;
        set(child_pctub, 'YData', [pub_percent(pest_idx)                   get(child_pctub, 'YData')]) ;
      else
        plot([xlim()], [pub_percent(pest_idx) pub_percent(pest_idx)], 'b--', 'DisplayName', 'pctub');
      end
    end
    if(not(plb_percent(pest_idx) == inf))
      if(length(findobj(get(sph, 'children'), 'DisplayName', 'pctlb')) > 0)
        child_pctlb = findobj(get(sph, 'children'), 'DisplayName', 'pctlb');
        set(child_pctlb, 'XData', [state.iteration                         get(child_pctlb, 'XData')]) ;
        set(child_pctlb, 'YData', [plb_percent(pest_idx)                   get(child_pctlb, 'YData')]) ;
      else
        plot([xlim()], [plb_percent(pest_idx) plb_percent(pest_idx)], 'b--', 'DisplayName', 'pctlb');
      end
    end

    if(pest_idx == 1)
      myxlim = [min(get(child_pest, 'XData')) max(get(child_pest, 'XData'))];
      if(myrange(myxlim) < 1)
        myxlim = [0 1];
      end
    end

%   if(not(plb_percent(pest_idx) == -inf))
%     child_pctlb = findobj(get(sph, 'children'), 'DisplayName', 'pctlb');
%     set(child_pctlb, 'XData', [state.iteration                         get(child_pctlb, 'XData')]) ;
%     set(child_pctlb, 'YData', [plb_percent(pest_idx)                   get(child_pctlb, 'YData')]) ;
%   end
    

  end
  grid on;
  xlim(myxlim);

  % Trying to set the limits on the y axis intelligently 
  % starting after the first iteration 
  if(state.iteration > 1)
    min_ydata = min(pest_data(pest_idx, :));
    max_ydata = max(pest_data(pest_idx, :));
    
    if(abs(min_ydata/plb_percent(pest_idx)) > 1.1)       % we're beyond the bound
      ymin = min_ydata;
    elseif(abs(min_ydata/plb_percent(pest_idx)) > 0.95)  % we're near the lower bound
      ymin = plb_percent(pest_idx);                      % so we make it the lower bound
    else                                                 % we're no where near the bound
      ymin = min_ydata;                                  % so we make it the lowest datum point
    end
    
    if(abs(max_ydata/pub_percent(pest_idx)) > 1.1)       % we're beyond the bound
      ymax = max_ydata;                                  
    elseif(abs(max_ydata/pub_percent(pest_idx)) > 0.95)  % we're near the upper bound
      ymax = pub_percent(pest_idx);                      % so we make it the upper bound
    else                                                 % we're no where near the bound
      ymax = max_ydata;                                  % so we make it the highest datum point
    end

    ylim([ymin-3, ymax+3]);
  end




end


% Pulling out the estimation data


% plotting the progress of the objective function.
sph = subplot(nrow, ncol, [(nrow*ncol-ncol+1) nrow*ncol]); hold on;
if((state.iteration == 0) & (length(get(sph, 'children')) == 0))
  plot(state.iteration, state.fval, 'r.-');
else
  children =get(sph, 'children');
  tmp = get(children(1), 'YData');
  abs(state.fval - tmp(end));
  set(children, 'XData', [state.iteration      get(children(1), 'XData')]) ;
  set(children, 'YData', [state.fval           get(children(1), 'YData')]) ;
  obj_data = fliplr(get(children, 'YData'));
  % if the objective function has changed a bunch
  % (spanning several logs) we switch to log scale
  if(log10(max(abs(obj_data))) - log10(min(abs(obj_data))) > 5)
    set(gca, 'yscale', 'log');
  end
end
grid on;

try
axis tight;
end


% forcing the figure to update after 
% each iteration
drawnow


if(strcmp(exit_when_stable, 'yes'))
   % we only check this after we've had at least iteration_history iterations
   if(state.iteration > iteration_history)
     slopes = [];
     % Calculating the slope of a trendline through the last
     % 'iteration_history' points and storing those in slopes. 
     % The value soln.m is the change in the parameter with 
     % respect to the iteration number
     for pest_idx = 1:length(pest)
       soln = linear_regression([0:iteration_history], pest_data(pest_idx,[(end-iteration_history):end]));
       slopes(end+1) = abs(soln.m);
     end

     % doing the same for the objective function
     soln = linear_regression([0:iteration_history], obj_data(1,[(end-iteration_history):end]));
     slopes(end+1) = abs(soln.m);

     % If the largest slope magnitude is less than 
     % slope tolerance then we stop the optimization
     if(max(slopes) < slope_tolerance)
       stop = true;
       vp(cfg,         '-----------------------------------');
       vp(cfg, sprintf('Function: %s.m ', cfg.estimation.monitor.status_function)) 
       vp(cfg,         'Estimation stopping criteria reached');
       vp(cfg, sprintf('Iteration history:      %d',iteration_history));
       vp(cfg, sprintf('slope_tolerance:        %s ', num2str_pretty(slope_tolerance)));
       vp(cfg, sprintf('Maximum slope of parameters and '));
       vp(cfg, sprintf('objective over history: %s ', num2str_pretty(max(slopes))))
       vp(cfg,         '-----------------------------------');
     end
   end
end
try
catch
  vp(cfg,         'estimation_status.m -- function failed');
end


function [val] = myrange(x)

  val = max(x) - min(x);
