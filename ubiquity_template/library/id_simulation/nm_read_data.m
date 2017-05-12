function [data]=nm_read_data(data);
% function [data]=nm_read_data(data);
% This function can be used to read NONMEM data sets and output tab files. To
% read a dataset it can be in excel or csv format. If an excel file is used,
% it should be saved in an older format (95 or 98). 
%  
%
% For data files in excel format it's necessary to specify the 
% file name and the sheet it's located in:
%  data.data_file.name  = 'data.xls'
%  data.data_file.sheet = 'data_sheet';
%  
% For data files in csv format:
%  data.data_file.name  = 'data.csv'
%  
% Notes on data file format:
%   - First line should contain header labels
%  
% For nonmem tab output files
%  data.data_file.name  = 'output.tab'
%  
% Notes on table file format:
%   - First line will be ignored
%   - The second lines should contain header labels
%  
%  % The following will automatically read in the headers:
%  data.column.auto               = 'yes';
%
%  % Otherwise you can manually specify the names you want 
%  % to use to reference the headers
%  data.column.names.evid         = 'EVID';  
%  data.column.names.observations = 'DV';    
%
%  % in this example only the columns named EVID and DV in 
%  % the data file will be accessible and they will be referenced 
%  % using the names evid and observations
%
%  With these fields populated, the file can be read in:
%  [data]=nm_read_data(data);
%

% reading in the data sheet:
if(regexpi(data.data_file.name, 'xls$'));
  [data.raw]=fetch_excel_raw(data.data_file.name, data.data_file.sheet);
elseif(regexpi(data.data_file.name, 'csv$'));
  data.raw=fetch_csv_raw(data.data_file.name);
elseif(regexpi(data.data_file.name, 'tab$'));
  dsoptions.ignore    = [1];
  dsoptions.delimiter = ' ';
  data.raw=fetch_csv_raw(data.data_file.name, dsoptions);
end


% sometimes the headers are surrounded by "quotes"
% and I'm striping these off because they don't show up in excel
data.header = data.raw(1,:);
data.values = data.raw(2:end,:);


% If auto header naming has been selected
% then we set that up
if(isfield(data, 'column'))
if(isfield(data.column, 'auto'))
if(strcmp(data.column.auto, 'yes'))
  for row_idx = 1:length(data.header)
    % replacing '.' with '_'
    data.header{row_idx} = regexprep(data.header{row_idx}, '\.', '_');
    eval(sprintf('data.column.names.%s = ''%s'';', data.header{row_idx}, data.header{row_idx}));
  end
end
end
end


% sometimes the tab files will have headers printed ever so often.
% we're removing those
if(regexpi(data.data_file.name, 'tab$'));
  vtmp = {};
  for row_idx=1:length(data.values(:,1))
    if( not(isstr(data.values{row_idx,1})))
      vtmp(end+1,:) = data.values(row_idx,:);
    end
  end
  data.values = vtmp;
end




% placeholder for figures later


% Data format: first row is header information, rows 2-n are records
% defining the mapping between column names and their column index 
% These are required column names:
column_names = fieldnames(data.column.names);

for column_name_idx = 1:length(column_names)

% column_names{column_name_idx}
% try
% catch
% sprintf('tmp_col = find_element(data.header,  data.column.names.%s);', column_names{column_name_idx})
% keyboard
% sdfsdf
% end

eval(sprintf('tmp_col = find_element(data.header,  data.column.names.%s);', column_names{column_name_idx}));

  if(tmp_col == -1)
    disp(sprintf('> Unable to find column in dataset: %s', getfield(data.column.names, column_names{column_name_idx})));
  else
    eval(sprintf('%s_col = tmp_col;', column_names{column_name_idx}));
    % storing these mappings:
    eval(sprintf('data.column.numbers.%s = %s_col;', column_names{column_name_idx}, column_names{column_name_idx}));
  end

end
