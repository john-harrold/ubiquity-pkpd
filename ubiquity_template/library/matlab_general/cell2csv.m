function []=cell2csv(cellarr, outputfile, delimiter)
%function []=cell2csv(cellarr, outputfile, delimiter)
%
% % example: 
%
% carr    = [{'a'} {'b,'} {'c'}
%              1   {'b'}     3
%            {'a'}   2   {'c'}
%            {'a'}   2     3 ];
% 
% ofile   = 'myfile.csv';
%  
%  
% cell2csv(carr, ofile)
%  
% Delimiter is optional. If it is specified, it will be used instead of a comma
%  
%  delim = '&';
%  
% cell2csv(carr, ofile, delim)
%  
%  

% defaulting to comma
use_del = ',';

if(exist('delimiter', 'var'))
    if(ischar(delimiter))
       use_del = delimiter;
    end
end
[num_rows, num_cols] = size(cellarr);

OFH = fopen(outputfile, 'wt');
for crow = 1:num_rows
  for ccol = 1:num_cols
    tmpstr = '';
    if iscellstr(cellarr(crow,ccol))
      %
      % current element is a string
      %
      tmpstr = char(cellarr(crow,ccol));
      % stripping out any delimiters from the string
      tmpstr = regexprep(tmpstr, use_del, '', 'ignorecase');
    else
      %
      % otherwise current element is a number
      %
      tmpstr = sprintf('%.4e', cell2mat(cellarr(crow,ccol)));
    end
    if(ccol == 1)
      fprintf(OFH, sprintf('%s', tmpstr));
    else
      fprintf(OFH, sprintf('%s%s', use_del, tmpstr));
    end
  end
  fprintf(OFH, sprintf('\n'));
end

fclose(OFH);


