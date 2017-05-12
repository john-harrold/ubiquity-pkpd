function [cfg]=system_load_dataset(cfg, data_name, data_file, data_sheet)
% function [cfg]=system_load_dataset(cfg, dsname, data_file, data_sheet)
%
% dsname = Short name given to the data set. This hsould begin with a letter and
%        can contain any combination of letters, numbers and _. No spaces  and
%        no dashes can be used.
%
% data_file = This is the name of the file and it can have the following
%        extensions: xls, csv, and tab. This file should be NONMEM-ish
%        with the column names in the first row.
%
% data_sheet = When the file type is xls, then the excel sheet can be specified as
%        an option.
%
d.column.auto           = 'yes';

% Checking to see if the data file exists:
if(exist(data_file, 'file') == 2)
  %loading the data file
  d.data_file.name        = data_file;
  if(regexpi(data_file, 'xls$'));
    d.data_file.sheet     = data_sheet;
  end
  [d]                     = nm_read_data(d);

  eval(sprintf('cfg.data.%s = d;', data_name));
else
  vp(cfg, sprintf('Error the data_file (%s) was not found',data_file))
end

