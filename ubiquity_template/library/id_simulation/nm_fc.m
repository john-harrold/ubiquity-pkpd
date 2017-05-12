function [col]=nm_fc(data, values, col_name)
% function [col]=nm_fc(data, values, col_name)
%
% fetches column from a NONMEM data set created by nm_read_data
%
%  data   - The value returned from nm_read_data
%
%  values - A cell array of data from the data. This can be either the entire
%  data set (data.values) or a subset that is returned from nm_select_records 
%
%  col_name - name of a column specified when the dataset was read (see the
%  help for nm_read_data for examples)
%

col = [];
if(isfield(data.column.numbers, col_name))
  cn = getfield(data.column.numbers, col_name);
  for ridx =1:length(values(:,1))
    if(isnumeric(values{ridx, cn}))
      col = [col; values{ridx, cn}];
    else
      error(sprintf('Non numeric value found in %s', col_name));
    end
  end

else
  error(sprintf('Error: unable to find column %s', col_name));
end

