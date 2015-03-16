function [str]=num2str_pretty(num)


if(num > 1000)
  str = sprintf('%.2e', num);
elseif(num < 1000 & num >= 100) 
  str = sprintf('%.0f', num);
elseif(num < 100 & num >= 1) 
  str = sprintf('%.1f', num);
elseif(num < 1 & num >= .1) 
  str = sprintf('%.2f', num);
else
  str = sprintf('%.2e', num);
end
