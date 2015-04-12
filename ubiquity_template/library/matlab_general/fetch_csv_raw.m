function [raw] = fetch_csv_raw(filename,options)
%
%
% Read in a delimited file 'filename' with the following options
% 
%  options.delimiter = string with the delimiter 
%  default is ','
%
%  options.ignore = array of line numbers to ignore
%  default is none
%

delimiter = ',';
ignore    = [];

if(exist('options', 'var'))
  if(isfield(options, 'ignore'))
    ignore = options.ignore;
  end
  if(isfield(options, 'delimiter'))
    delimiter = options.delimiter;
  end
end

%----
% The first part was lifted from this web page:
%     http://stackoverflow.com/questions/4747834/import-csv-file-with-mixed-data-types
% Provided by user gnovice:
%     http://stackoverflow.com/users/52738/gnovice
% This is distributed according to the creative commons license:
%    http://creativecommons.org/licenses/by-sa/3.0/
%

%
% This reads in the csv file with each element a string
%
fid = fopen(filename,'r');   %# Open the file
lineArray = {};              %# initializing cell array (ideally slightly
                             %#   larger than is needed)
lineIndex = 1;               %# Index of cell to place the next line in
nextLine = fgetl(fid);       %# Read the first line from the file
while ~isequal(nextLine,-1)         %# Loop while not at the end of the file
  % only reading in lines that are not ignored
  if(not(sum(ismember(ignore, lineIndex))))
    lineArray(end+1,1) = {nextLine};  %# Add the line to the cell array
  end
  lineIndex = lineIndex+1;         %# Increment the line index
  nextLine = fgetl(fid);            %# Read the next line from the file
end
fclose(fid);                 %# Close the file

nrows = length(lineArray(:,1));

for iLine = 1:nrows
  if(strcmp(delimiter, ' '))
    lineArray{iLine} =   regexprep(lineArray{iLine}, '^\s+', '') ;
    lineArray{iLine} =   regexprep(lineArray{iLine}, '\s+$', '') ;
    lineArray{iLine} =   regexprep(lineArray{iLine}, '\s+', ' ') ;
  end
  lineData = textscan(lineArray{iLine},'%s',...  %# Read strings
                      'Delimiter',delimiter);
  lineData = lineData{1};              %# Remove cell encapsulation
  if strcmp(lineArray{iLine}(end),delimiter)  %# Account for when the line
    lineData{end+1} = '';                     %#   ends with a delimiter
  end
  lineArray(iLine,1:numel(lineData)) = lineData;  %# Overwrite line data
end
%----



% to make it consitent with the fetch_excel_raw function
% we need to convert the numeric values currently represented
% by strings into numbers
[nr, nc] = size(lineArray);

raw = cell(nr, nc);

for row_idx = 1:nr
for col_idx = 1:nc

  tmp_element = lineArray{row_idx, col_idx};

  if(not(isnan(str2double(tmp_element))))
    tmp_element =str2double(tmp_element);
  end

  raw{row_idx, col_idx} = tmp_element;
end
end




