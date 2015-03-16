function [str] = make_number_pretty(number)
% function [str] = make_number_pretty(number)
%
%  takes a real number and turns it into a readable string.
%
%


if(number == 0.0)
    str = '0';
elseif((abs(number) < .099) | (abs(number) > 999 ))
    str = sprintf('%.2e', number);
else
    str = sprintf('%.3f', number);
end


