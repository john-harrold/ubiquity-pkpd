function [result] = AUC_trapezoid(data, extrap)
% function [result] = AUC_trapezoid(data, extrap)
% 
%  Calculating the
% 
%   with 'k' sample times  and 'r' subjects
%   time point the array 'data' is defined as:
%                  
%          % time  subj1  subj2... subjectr
%   data = [ t1    u11    u12  ... u1R
%            t2    u21    u22  ... u2R     
%            t3    u31    u32  ... u3R   
%            t4    u41    u42  ... u4R  
%            .     .      .        .   
%            .     .      .        .
%            .     .      .        .
%            tK    uK1    uK2  ... uKR ];
%                  
%  since it is an array, if there is less than r samples at each time point
%  a '-1' place holder should be used.
% 
%  If extrap is true we'll extrapolate final time to infinity 
% 




t   = data(:,1);
obs = data(:,2:end);
[R] = length(obs(1,:));   % number of subjects
[K] = length(t);

result.AUC_sub = zeros(R,1);
result.df      = R;

for sub_idx =1:R
  %getting rid of the missing values for this subject (-1's)
  [ttmp, subtmp] = strip_missing(t, obs(:,sub_idx));

  % adding the first and last pieces of the curve
  result.AUC_sub(sub_idx) = result.AUC_sub(sub_idx)...
                          + subtmp(1)        *(ttmp(2,  1)  - ttmp(1      ,1))/2;  % first part
                          + mean(subtmp(end))*(ttmp(end,1)  - ttmp(end-1  ,1))/2;  % last  part

  % now adding all the pieces in between
  for time_idx = 2:length(ttmp)-1
     result.AUC_sub(sub_idx) = result.AUC_sub(sub_idx) ...
                             + subtmp(time_idx)*(ttmp(time_idx+1)-ttmp(time_idx-1))/2;
  end

  if(exist('extrap', 'var'))
    if(extrap)
      if(subtmp(end) > subtmp(end-1))
        disp(' unable to extraploate to infinity');
      else
        %           ln(y2) - ln(y1)
        %  kel = -  ---------------
        %               t2 - t1
        kel = -(log(subtmp(end)) - log(subtmp(end-1)))/(ttmp(end) - ttmp(end-1));
        result.AUC_sub(sub_idx) = result.AUC_sub(sub_idx) + subtmp(end)/kel;
      end
    end
  end
end


result.AUC       = mean(result.AUC_sub);
result.var_AUC   =  var(result.AUC_sub);

result.confidence_interval    = [result.AUC+tinv(.975,result.df)*sqrt(result.var_AUC);
                                 result.AUC-tinv(.975,result.df)*sqrt(result.var_AUC)];
