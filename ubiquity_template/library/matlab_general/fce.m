function [value] = fce(ca, row_header, col_header)
% function [value] = fce(ca, row_header, col_header)
%
% fce - fetch cell element
%
% % This function takes ca (cell array) and returns an element found on 
% % the row with the first entry 'row_header' (string) and the column with the first entry
% % 'col_header' (string).
%
% % For example consider the cell array mydata below:
%
%  mydata = [{'group'}  {'CL'} {'Vc'} {'Fb'}
%            {'  g1 '}   .21     75     .5
%            {' g2'}     .32     54     .7
%            {' g3 '}    .45     63     .2];
%
% % If you wanted to get the value on row g2 and column Vc you would execute
% % the following:
%
%  fce(mydata, 'g2', 'Vc')
%

% checking the input data
row_index = find_element(ca(:,1), row_header);
col_index = find_element(ca(1,:), col_header);


if(row_index > 0 & col_index > 0)
  value = ca{row_index, col_index};
else
  disp(sprintf('Unable to find row: %s or column: %s in the cell array', row_header, col_header));
end
