function [th]=calculate_halflife(times, values, tmin, tmax)
% function [th]=calculate_halflife(tsub, vsub)
%  tsub - timecourse of subject
%  vsub - observed values for subject
%  tmin - optional lower bound on values in tsub to consider
%  tmax - optional upper bound on values in tsub to consider



if(~exist('tmin', 'var'))
   tmin = min(times);
elseif(isempty(tmin))
   tmin = min(times);
end


if(~exist('tmax', 'var'))
   tmax = max(times);
elseif(isempty(tmax))
   tmax = max(times);
end


tsub =  times((times >= tmin) & (times <= tmax));
vsub = values((times >= tmin) & (times <= tmax));


lr       = linear_regression(tsub, log(vsub));
th.thalf = -log(2)/lr.m;
th.vsub  = vsub;
th.tsub  = tsub;
th.lr    = lr;

