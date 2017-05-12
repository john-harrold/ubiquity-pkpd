function [cfg]=system_set_bolus(cfg, state, times, values)
% function [cfg]=system_set_bolus(cfg, state, times, values)
%
% Overwrites the current bolus inputs for 'state' with 
% the supplied 'times' and % 'values'
%
%
% cfg    - System configuration variable 
% state  - State name to receive the bolus 
% times  - Vector of times the bolus's are administered
% values - Vector of bolus values corresponding to the times
%
% To a value of 10 at dose at times 6 and 13 use the following:
%
% cfg=system_set_bolus(cfg, 'Cp',    ...
%                           [ 6 13], ...
%                           [10 10]);
%



bolus_good = true;
% Checking inputs and system specifications 
if(isfield(cfg.options.inputs, 'bolus'))
  if(isfield(cfg.options.inputs.bolus, 'species'))
    
    if(~isfield(cfg.options.inputs.bolus.species, state))
      % Setting the bolus_good to indicate the cfg variable has bolus information
      vp(cfg, sprintf('The specified state (%s) was not found', state)) 
      bolus_good = false;
    end

    if(~(length(times) == length(values)))
      vp(cfg, sprintf('The number of time entries (%d) ', length(times)));
      vp(cfg, sprintf('is different from the number of  ' ));
      vp(cfg, sprintf('value entries (%d).', length(values)));
      bolus_good = false;
    end
    
  else
    bolus_good = false;
  end
else
  bolus_good = false;
end


% if the system and inputs check out then we update 
% the contents of the bolus fields
if(bolus_good)

  % storing the old values 
  bolus_old = cfg.options.inputs.bolus;

  % getting all of the species configured to receive a bolus
  species = fieldnames(cfg.options.inputs.bolus.species);


  % next all possible times are stored
  all_times = unique(sort([bolus_old.times.values, times]));
  % least one nonzero bolus being applied:
  all_times_keep = [];

  for time_index_all = 1:length(all_times)
    if(sum(all_times(time_index_all) == times)) 
      % if the current time is in the list 
      % of specified bolus times then we keep it
      all_times_keep(end+1) = all_times(time_index_all);
    elseif(sum(all_times(time_index_all) == bolus_old.times.values))
      % otherwise we loop through any other states 
      % and see if there is a nonzero dose 
      time_index = (all_times(time_index_all) == bolus_old.times.values);
      for(spidx=1:length(species))
        eval(sprintf('bolus_value = bolus_old.species.%s.values(time_index);', species{spidx}));
        if((bolus_value ~= 0) & (~strcmp(state, species{spidx})))
          all_times_keep(end+1) = all_times(time_index_all);
        end
      end
    end
  end

  % stripping out any duplicates
  all_times_keep = sort(unique(all_times_keep));


  % Now we create a new data structure for the bolus information
  % by copying the old one and zeroing out the values
  bolus_new = cfg.options.inputs.bolus;
  bolus_new.times.values = all_times_keep;
  for(spidx=1:length(species))
    eval(sprintf('bolus_new.species.%s.values = [];', species{spidx}));
  end


  % now we populate it with the nonzero values

  for(time_index_keep = 1:length(all_times_keep))
    for(spidx=1:length(species))
      % defaulting to 0
      species_bolus_value = 0;
      if(strcmp(species{spidx}, state))
        % if the currently kep time is present in the times vector 
        % we then overwrite the default value of zero with the value
        % that corresponds to that time
        time_index = (all_times_keep(time_index_keep) ==  times);
        if(sum(time_index))
          species_bolus_value = values(time_index);
        end
      else
        time_index = (all_times_keep(time_index_keep) ==  bolus_old.times.values);
        if(sum(time_index))
          eval(sprintf('species_bolus_value = bolus_old.species.%s.values(time_index);', species{spidx}));
        end
      end
      eval(sprintf('bolus_new.species.%s.values(end+1) = species_bolus_value;', species{spidx}));
    end
  end

  % replacing the bolus_new back in the cfg data structure
  cfg.options.inputs.bolus = bolus_new;

else
  vp(cfg, 'system_set_bolus()');
  vp(cfg, 'There was an error and the bolus information');
  vp(cfg, 'was not set. See above.');
end


