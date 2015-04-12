function [data]=ed(ca, header)
%
% function [data]=ed(ca, header)
%
% ed - extract data
%
%  % This is a function used to extract columns of data from a cell array (ca).
%  % The header (string) specifies the column heading you wish to extract. The
%  % function returns 'data' a column vector of all the floating point numbers in
%  % the column designated by 'header' 
%   
%  % Consider the cell array mydata:
%  mydata = {'time'   'BW'  'Cp'
%            'day '   'kg'  'ng/ml'
%             1        50    25
%             10       55    26  
%             30       59    24  
%             40       62    27  
%             50       62    25};
%
%
% % To pull out the time and BW data you would use the following:
% myBW = ed(mydata, 'BW');
% mytime = ed(mydata, 'time');
%
% % These can then be plotted:
% plot(mytime, myBW, 'o'); 

data = [];

data_columns = strcmp((ca(1,:)),  header);

if(sum(data_columns) < 1)
  disp(sprintf('#->Error unable to find the specified header (%s)', header));
elseif(sum(data_columns) > 1)
  disp(sprintf('#->Error multiple headers found (%s)', header));
else
  % pulling out the numerical data
  dv = ca(:,data_columns);
  for idx =1:length(dv)
    if(isnumeric(dv{idx}))
      data = [data; dv{idx}];
    end
  end
end

