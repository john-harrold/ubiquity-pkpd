function [mystr] = var2string_gen(var)

mystr = 'Unknown variable type';

if(isnumeric(var))
  if(length(var) == 1)
    mystr = var2string(var, 0);
  else
    mystr = sprintf('min = %s; max = %s; length = %d', ...
            var2string(min(var), 0), ...
            var2string(max(var), 0), ...
            length(var));
  end
end
if(islogical(var))
  if(var)
    mystr = 'true';
  else
    mystr = 'false';
  end
end

if(isstr(var))
 mystr = var;
end
