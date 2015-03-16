function [elnumber]=find_element(cellarr, valueid, insensitive)
% function [elnumber]=find_element(cellarr, valueid, insensitive)
%
%  This function takes the one dimensional cell array (cellarr) and looks at
%  each element to see if it is the same as valueid and returns the element
%  number corresponding to that value. 
%
%  This is used to determine the row or column number which correpsonds to a
%  specific value (valueid). If valueid is a string, then a string comparison
%  will be used. If valueid is a number then it will be assuemd that cellarr 
%  is full of numeric values.
%
%  elnumber will have a value of -1 if no solution was found
%
%  Note: Any spaces at the beginning or end of the elements of cellarr are
%  ignored, and the search is case sensitive by default. If you define
%  'insensitive' as equal to 'yes', then it will perform a case insensitve
%  search.
%
%  Example: the following cell array contains a table of information about
%  three groups. Say you want clearance values for group 2. You can access
%  this value by looking at mydata(3,2). But it is undesirable to hard code
%  those values, a better way is to do the following:
%
%  mydata = [{'group'}  {'CL'} {'Vc'} {'Fb'}
%            {'  g1 '}   .21     75     .5
%            {' g2'}     .32     54     .7
%            {' g3 '}    .45     63     .2];
%
%
% % all of the rows in the first column
% tmp_row = find_element(mydata(:,1), 'g2');
%
%
% % all of the columns in the first row
% tmp_col = find_element(mydata(1,:), 'CL');
%
% mydata(tmp_row,tmp_col)
%


if not(exist('insensitive', 'var'))
  insensitive = 'no';
end

  elnumber    = -1;
  if(isstr(valueid))
    % stripping any extra whitespace at the beginning or end of the string
    valueid = regexprep(valueid, '^\s*', ''); % beginning white space
    valueid = regexprep(valueid, '\s*$', ''); % end       white space
  end

  for arr_idx = 1:length(cellarr)
    if(isstr(valueid) & ischar(cellarr{arr_idx}))
      tmpstr = char(cellarr(arr_idx));
      % stripping any extra whitespace at the beginning or end of the string
      tmpstr = regexprep(tmpstr, '^\s*', ''); % beginning white space
      tmpstr = regexprep(tmpstr, '\s*$', ''); % end       white space
      % now we check to see if the strings match
      if(strcmp(insensitive, 'yes'))
        if(strcmpi(tmpstr, valueid))
          elnumber = arr_idx ;
        end
      else
        if(strcmp(tmpstr, valueid))
          elnumber = arr_idx ;
        end
      end
    elseif(isnumeric(valueid) & isnumeric(cellarr{arr_idx}))
      tmpval = cell2mat(cellarr(arr_idx));
      if(valueid == tmpval)
        elnumber = arr_idx;
      end
    end

  end
