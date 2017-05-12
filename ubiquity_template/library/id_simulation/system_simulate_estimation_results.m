function [odp]=system_simulate_estimation_results(pest, cfg)
% function [odp]=system_simulate_estimation_results(pest, cfg)
%
% odp - Data structure containing observation and prediction 
%       information for each output in each cohort. This is used to 
%       generate figures with smooth profiles.
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
%
% If you prefer a flat data file the following field:
%
%   odp.flat
%
% can be saved as a csv file using the command:
%
% cell2csv(odp.flat, 'myfile.csv')
%
%

[od, odp]=system_od_general(pest, cfg);


% creating the flat format for exporting
flat = {'time' 'obs' 'pred' 'var' 'res' 'wres' 'cohort' 'output' 'type'};
% first we iterate through each cohort
cohorts = fieldnames(odp.cohorts);
for(chidx = 1:length(cohorts))
  cohort = getfield(odp.cohorts, cohorts{chidx});
  outputs = fieldnames(cohort);

  % next we go through each output
  for(opidx = 1:length(outputs))
    output = getfield(cohort, outputs{opidx});

    % lastly we go throgh each record for the cohort/output combination
    for(rcidx = 1:length(output.od.time))
      res  = output.od.obs(rcidx)  -  output.od.pred(rcidx);
      wres = (output.od.obs(rcidx) -  output.od.pred(rcidx))/output.od.pred(rcidx);
      %    
      %                time                         obs                   pred                         var                  res wres    cohort         output         type
      flat(end+1,:) = {output.od.time(rcidx)        output.od.obs(rcidx)  output.od.pred(rcidx)        output.od.var(rcidx) res wres   cohorts{chidx} outputs{opidx} 'record'};
    end
    for(rcidx = 1:length(output.od.time_smooth))
       flat(end+1,:) = {output.od.time_smooth(rcidx) -1                    output.od.pred_smooth(rcidx) -1                   -1  -1     cohorts{chidx} outputs{opidx} 'smooth'};
    end
  end
end

odp.flat = flat;
