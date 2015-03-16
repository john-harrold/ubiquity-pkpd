function str = var2string(var,maxlength) 
%  function str = var2string(var, 12) 
%  converts the numerical value 'var' to a padded string 12 characters wide


if((var < .01 )| (var > 999))
    str = sprintf('%.3e', var );
else
    str = sprintf('%.3f    ', var );
end

str = pad_string(str, maxlength);


