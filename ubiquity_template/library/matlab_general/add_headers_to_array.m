function [mycell] = add_headers_to_array(myarray, myheaders)
%
% function [mycell] = add_headers_to_array(myarray, myheaders)
% 
% Takes the cell array headers and the numerical array myarray and combines
% them into a cell array with the headers as the first row
%  
% Example:
%  
% myarray  = [1  2  3  4
%             5  6  7  8
%             9 10 11 12];
% 
% myheaders = {'col a' 'col b' 'col c' 'col d'};
%  
% 
%
% [mycell] = add_headers_to_array(myarray, myheaders)
%
%
%  mycell should look like this:
%
%   'col a'    'col b'    'col c'    'col d'
%   [    1]    [    2]    [    3]    [    4]
%   [    5]    [    6]    [    7]    [    8]
%   [    9]    [   10]    [   11]    [   12]



% the first row of mycell is myheaders
% this assignment a = b(:)' will make 
% sure it's a row vector
myheaders = myheaders(:)';
mycell = {};

% checking the number of columns in header 
% and array to make sure they are the same
if(length(myarray(1,:)) == length(myheaders))
  mycell = [myheaders; num2cell(myarray)];
else
  disp('Error the mylength of the header should be equal to the number of columns in myarray');
end

