function str = var2string(var,maxlength, mydp) 
%  function str = var2string(var, maxlength, mydp) 
%
%
%  Converts the numerical value 'var' to a padded string maxlength characters wide
%  with mydp (optional default=3) digits to the left of the decimal. 
%
%



if(~exist('mydp', 'var'))
  mydp = 3;
end


if( var == floor(var) )
   str = sprintf('%d', var);
elseif((var < .01 )| (var > 999))
     eval(sprintf('str = ''%s'';',eval(sprintf('sprintf(''%%.%de'', var)',mydp))))
else
     eval(sprintf('str = ''%s'';',eval(sprintf('sprintf(''%%.%df'', var)',mydp))))
end

str = pad_string(str, maxlength);


