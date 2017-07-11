function [chkres]=ubiquity_name_check(test_name)
%
% Error checking function to make sure the test_name 
% matches the following rules:
%
%  - starts with a letter
%  - only conatins letters, numbers, and _
%

 chkres.isgood = 1;
 chkres.msgs   = {};
 

 % Making sure it starts with a letter
 if(length(regexp(test_name, '^[a-z,A-z]')) ~= 1)
   chkres.msgs(end+1) = {'Does not begin with a letter'};
 end

 % now we remove all of the allowed characters and see what's left
 % there should be nothing left :)
 test_name_trim = regexprep(test_name, '[a-z,A-z,0-9,_]', '');

 if(length(test_name_trim)  > 0)
   chkres.msgs(end+1)  = {'Should only contain letters, numbers and _'};
 end

 % If there are any messages we flip the isgood to 
 % false and concatenate them together

 if( length(chkres.msgs) > 0)
   chkres.isgood = 0;
   chkres.msg    = strjoin(chkres.msgs, ', ');
   
 end

