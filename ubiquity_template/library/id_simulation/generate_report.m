function [files]=generate_report(parameters, solution_stats, cfg)
% function []=generate_report(parameters, solution_stats, cfg)
%
%   parameters     -- vector of parameter estimates from fminsearch
%   solution_stats -- output of solution_statistics.m
%   cfg            -- variable with parameter field defined
%
%   creates file 'report.txt' with summary information from analysis
%


report_file         = sprintf('output%sreport.txt', filesep) ;
parameters_all_file = sprintf('output%sparameters_all.csv', filesep) ;
parameters_est_file = sprintf('output%sparameters_est.csv', filesep) ;


header_1 = {'pname' 'guess'  'estimate' 'cvpct' 'cilb' 'ciub' 'units' 'notes'};

p_all(1,:) = header_1;
p_est(1,:) = header_1;

% parameter, original estimate, 95% CI, 95% CI, Units,
% name,      guess,           , Lower,   Upper,      ,



RFH = fopen(report_file, 'w');

% Start
% Dumping Estimates, confidence interval
%
fprintf(RFH, '         Estimate     95 %% Confidence Interval        Coeff. of Var  Notes \n');
fprintf(RFH, '                      Lower Bound    Upper Bound       (Percent)       \n');


% Dumping Estimates, confidence interval, etc.
% merging fixed and estimated parameters
for pidx=1:length(cfg.parameters.names)
  pname = cfg.parameters.names{pidx};
  if(isfield(cfg.estimation.mi, pname))
    % processing estimated parameters
    est_idx = getfield(cfg.estimation.mi, pname);
    notes = compare_estimate(cfg, parameters, pname);
    fprintf(RFH, '%s', pad_string(cfg.parameters.names{pidx}, 10));
    fprintf(RFH, '%s', var2string(parameters(est_idx), 12));
    fprintf(RFH, '%s', var2string(solution_stats.confidence_interval.lower_bound(est_idx),15));
    fprintf(RFH, '%s', var2string(solution_stats.confidence_interval.upper_bound(est_idx),15));
    fprintf(RFH, '%s', var2string(solution_stats.coefficient_of_variation(est_idx),15));
    fprintf(RFH, '%s', pad_string(notes, 5));
    fprintf(RFH, '\n');
    tmp_row = {pname                                                    ...  % name
               cfg.parameters.values(pidx)                              ...  % original
               parameters(est_idx)                                      ...  % estimate
               solution_stats.coefficient_of_variation(est_idx)         ...  % CV %
               solution_stats.confidence_interval.lower_bound(est_idx)  ...  % CI lower 
               solution_stats.confidence_interval.upper_bound(est_idx)  ...  % CI upper 
               cfg.parameters.units{pidx}                               ...  % units
               notes};                                                  ...  % notes
    p_est(end+1,:) = tmp_row;
  else
    notes = 'F';
    fprintf(RFH, '%s', pad_string(cfg.parameters.names{pidx}, 10));          % name     
    fprintf(RFH, '%s', var2string(cfg.parameters.values(pidx), 12));         % original 
    fprintf(RFH, '%s', pad_string('', 15));                                  % estimate (placeholder)
    fprintf(RFH, '%s', pad_string('', 15));                                  % CV %     (placeholder)
    fprintf(RFH, '%s', pad_string('', 15));                                  % CI lower (placeholder)
    fprintf(RFH, '%s', pad_string('', 15));                                  % CI upper (placeholder)
    fprintf(RFH, '%s', pad_string(notes, 5));                                % units    
    fprintf(RFH, '\n');                                                      % notes    

    tmp_row = {pname                                                    ...  % name
               cfg.parameters.values(pidx)                              ...  % original
               '---'                                                    ...  % estimate
               '---'                                                    ...  % CV %
               '---'                                                    ...  % CI lower 
               '---'                                                    ...  % CI upper 
               cfg.parameters.units{pidx}                               ...  % units
               notes} ;                                                 ...  % notes
    % processing fixed
  end
  p_all(end+1,:) = tmp_row;
end
notes_str = 'F=Fixed parameter, L=estimate at/near lower bound, U=estimate at/near upper bound'; 
p_all{end+1,1} = notes_str;
p_est{end+1,1} = notes_str;
fprintf(RFH, '---\n');
fprintf(RFH, notes_str);
fprintf(RFH, '\n');


fprintf(RFH, '\n\n\n');

% Start
% Dumping the covariance matrix
%
fprintf(RFH, '%s', pad_string('', 10));
fprintf(RFH, 'Variance -- Covariance Matrix \n');

fprintf(RFH, '%s', pad_string('', 10));
for i=1:length(parameters)
    fprintf(RFH, '%s', pad_string(sprintf('%s    ',char(cfg.estimation.parameters.names(i))), 12));
end
    fprintf(RFH, '\n');

for i=1:length(parameters)
    fprintf(RFH, '%s', pad_string(char(cfg.estimation.parameters.names(i)), 10));
for j=1:length(parameters)
    if(j<=i)
    fprintf(RFH, '%s', var2string(solution_stats.covariance(i,j),12));
    end
end
    fprintf(RFH, '\n');
end
%
% Dumping the covariance matrix
% End

fprintf(RFH, '\n\n\n');
fprintf(RFH, 'Misc Information \n');
fprintf(RFH, 'OBJ = %.6f\n',solution_stats.objective);
fprintf(RFH, 'AIC = %.6f\n',solution_stats.aic);
fprintf(RFH, 'BIC = %.6f\n',solution_stats.bic);





fclose(RFH);

% dumping parameter information to csv files
cell2csv(p_all, parameters_all_file);
cell2csv(p_est, parameters_est_file);

vp(cfg,         'Report generated and placed in: ');
vp(cfg, sprintf('   %s', report_file));
vp(cfg,         'Estimated parameter information ');
vp(cfg,         'summarized in CSV format: ');
vp(cfg, sprintf('   %s', parameters_est_file));
vp(cfg,         'All parameter information ');
vp(cfg,         'summarized in CSV format: ');
vp(cfg, sprintf('   %s', parameters_all_file));

files.report_file         = report_file         ;
files.parameters_all_file = parameters_all_file ;
files.parameters_est_file = parameters_est_file ;

function [notes]  = compare_estimate(cfg, parameters, pname)
%
% ceecking to see if the estimated parameter pname with the value in the
% parameters vector is close to the upper or lower bounds in
% cfg.parameters.upper_bound or cfg.parameters.lower_bound)
%

notes = '';

pidx        = getfield(cfg.options.mi.parameters, pname);
est_idx     = getfield(cfg.estimation.mi, pname);
pvalue      = parameters(est_idx);
lower_bound = cfg.parameters.lower_bound(pidx);
upper_bound = cfg.parameters.upper_bound(pidx);

lower_diff = abs(lower_bound - pvalue);
upper_diff = abs(upper_bound - pvalue);

if(lower_diff/lower_bound  <0.05)
notes = 'L';
elseif(upper_diff/upper_bound  <0.05)
notes = 'U';
end




