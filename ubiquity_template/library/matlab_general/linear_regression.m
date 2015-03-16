function [soln] = linear_regression(xdata, ydata)
%function [soln] = linear_regression(xdata, ydata)
%
% Takes the vectors xdata and ydata and returns a least squares fit 'soln'
% Which has the following form:
% 
%   soln.m         fitted slope
%   soln.b         fitted intercept
% 
%   soln.equation  string with the y=mx+b equation
%   soln.Rsq       R^2 coefficient of determination 
%   soln.xpred     x values (sorted in ascending order)
%   soln.ypred     corresponding y values using the linear prediction
%

if(isrow(xdata))
  xdata = xdata';
end
if(isrow(ydata))
  ydata = ydata';
end

data_sorted = sort([xdata, ydata], 1);

xdata = data_sorted(:,1);
ydata = data_sorted(:,2);

% finding the least squares parameters
params =  [xdata ones(size(xdata))]\ydata;

m = params(1);
b = params(2);

ymean = mean(ydata);
ypred = m*xdata + b;

SStot = sum((ypred-ymean).^2);
SSres = sum((ypred-ydata).^2);

Rsq = 1-SSres/SStot;


%giving the solution in y=mx+b form
if(b < 0)
 soln.equation = sprintf('y = %sx %s', var2string(m, 1), var2string(b, 1));
else
 soln.equation = sprintf('y = %sx + %s', var2string(m, 1), var2string(b, 1));
end
soln.equation = regexprep(soln.equation, ' ', '');

soln.Rsq   = Rsq;
soln.m     = m;
soln.b     = b;
soln.xpred = xdata;
soln.ypred = ypred;
