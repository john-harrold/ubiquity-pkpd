function [mycsv]=sim2csv(som, cfg, csvfile)
% function [mycsv]=sim2csv(som, cfg, csvfile)
%
% Takes the contents of run_simulation_ubiquity in som and places the
% timescales, states, and outputs into csvfile. The contents of that file are
% also returned in mycsv
%
% Example
%   sim2csv(som, cfg, 'test.csv')
%





mycsv = {};

TS_names =  fieldnames(som.times);
for TS_idx = 1:length(TS_names)
  TS_name   = TS_names{TS_idx};
  TS_values = getfield(som.times, TS_name);

  mycsv(:,end+1) =  [{TS_name}
                     num2cell(TS_values)];


end

state_names =  fieldnames(som.states);
for state_idx = 1:length(state_names)
  state_name   = state_names{state_idx};
  state_values = getfield(som.states, state_name);

  mycsv(:,end+1) =  [{state_name}
                     num2cell(state_values)];

end

output_names =  fieldnames(som.outputs);
for output_idx = 1:length(output_names)
  output_name   = output_names{output_idx};
  output_values = getfield(som.outputs, output_name);

  mycsv(:,end+1) =  [{output_name}
                     num2cell(output_values)];

end


cell2csv(mycsv, csvfile);



