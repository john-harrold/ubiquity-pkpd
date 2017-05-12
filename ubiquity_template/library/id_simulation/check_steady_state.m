function []=check_steady_state(som)
% function []=check_steady_state(som)
% Checks to make sure the system is running at steady state
% where som is the mapped output from run_simulation_ubiquity 



% getting a list of all the states in the system:
states       = fieldnames(som.states);
offset_found = 'no';
for state_idx = 1:length(states)
    state = getfield(som.states, states{state_idx});

    state_max = max(abs(state));

    if state_max > 0
      if(range(state)/state_max > eps*100)
        if(strcmp(offset_found, 'no'))
          disp(sprintf(' Possible steady state offset'))
          disp(sprintf(' range       |             | state'))
          disp(sprintf(' (max-min)   | max(abs(s)) | name '))
          disp(sprintf('------------------------------------'))
          offset_found = 'yes';
        end
          disp(sprintf(' %.3e   | %.3e   | %s', range(state), state_max, states{state_idx}));
      end
    end
end

if(strcmp(offset_found, 'no'))
   disp(sprintf('|> No steady state offsets found'))
else
   disp(sprintf('#>------------------------------------'))
   disp(sprintf('#> Deviations from steady state found'))
   disp(sprintf('#> See above for details'))
   disp(sprintf('#> Machine Precision is %.3e', eps))
end
