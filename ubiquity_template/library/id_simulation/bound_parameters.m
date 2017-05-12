function [params, errorflag, objmult]=bound_parameters(params, lb, ub);
% function [params, errorflag, objmult]=bound_parameters(params, lb, ub);
%
%  params - vector of parameters
%  lb     - vector of corresponding lower bounds
%  ub     - vector of corresponding upper bounds
%
%  errorflag - if this value is 1 then one or more of the upper bounds are
%  less than the lower bounds.
%
%  objmult - multiplier for the objective function.
%

errorflag = 0;

objmult = 0.0;

if(sum(lb > ub))
errorflag = 1;
end
for i=1:length(params)
    if(params(i) > ub(i))
        objmult = objmult + 10*exp(abs(params(i) - ub(i)));
        params(i) = ub(i);
    elseif(params(i) < lb(i))
        objmult = objmult + 10*exp(abs(params(i) - lb(i)));
        params(i) = lb(i);
    end

end

