function [result] = AUC_Bailors_method(data)
% function [result] = AUC_Bailors_method(data)
% 
% 
%  This is an implementation of Bailors method for calculating AUCs with
%  sparse sampling. It is taken from the following publication:
% 
%  Nedelman, J. R., Gibiansky, E., & Lau, D. T. (1995). Applying Bailer's
%  method for AUC confidence intervals to sparse sampling Pharmaceutical
%  Research, 12(1), 124-128.
% 
% 
%   with 'k' sample times and a maximum of 'r' possible observations at each
%   time point the array 'data' is defined as:
% 
%          % time  obs1 obs2 ... obsr
%   data = [ t1    u11  u12  ... u1R
%            t2    u21  u22  ... u2R     
%            t3    u31  u32  ... u3R   
%            t4    u41  u42  ... u4R  
%            .     .    .        .   
%            .     .    .        .
%            .     .    .        .
%            tK    uK1  uK2  ... uKR ];
% 
% 
%  since it is an array, if there is less than r samples at each time point
%  a '-1' place holder should be used.
% 
% 
%  result is a data structure with the following fields
% 
%   result.AUC        mean AUC 
%   result.sigma_AUC  variance from the mean
%   result.r          vector of sampels per time
%   result.df         degrees of freedum
%
%   % Upper and lower bounds on the confidence intervals
%   result.confidence_interval
% 
% 
% 


% calculating the weight vector
% 
%   w1 = (t(2)   -   t(1)  )/2 
%   wk = (t(k+1) -   t(k-1))/2    k = [2,K-1];
%   wk = (t(K)   -   t(K-1))/2    
% 


t   = data(:,1);
obs = data(:,2:end);

% Its' possible for a row of data to have no observations. I there are rows
% with only -1's then we remove those (times first then observations)
%t   =   t(boolean(sum((obs >0),2)),:);
%obs = obs(boolean(sum((obs >0),2)),:);
t   =   t((sum((obs >0),2)>0),:);
obs = obs((sum((obs >0),2)>0),:);
[K] = length(t);


w   = zeros(K,1);
r   = zeros(K,1);
u   = zeros(K,1);
ssq = zeros(K,1);

w(1) = (t(2,1)-t(1  ,1))/2;
w(K) = (t(K,1)-t(K-1,1))/2;
for k_idx=1:K
    % calculating weights
    clear nz_elements; % temporary variable to hold non zero elements
    if((k_idx > 1) & (k_idx<K))
      % calculating the middle weights
      w(k_idx) = (t(k_idx+1)-t(k_idx-1))/2;
    end  
    nz_elements     = obs(k_idx, obs(k_idx,:)>0);

    % samples per timepoint
    r(k_idx)        = length(nz_elements);

    % timepoint mean
    u(k_idx)        = mean(nz_elements);

    % timepoint variiacne
    ssq(k_idx) = var(nz_elements);

end


%keyboard

result.AUC       = sum(w.*u);
result.var_AUC   = sum(w.^2.*ssq./r);
result.r         = r;
% degrees of freedom calculated by 6a according to
% the reference above
df_num           = sum(w.^2.*ssq.^2./r).^2;
df_denom         = sum(w.^4.*ssq.^4./(r.^2.*(r-1)));
result.df        = df_num./df_denom;

%calculating confidence intervals
result.confidence_interval    = [result.AUC+tinv(.975,result.df)*sqrt(result.var_AUC);
                                 result.AUC-tinv(.975,result.df)*sqrt(result.var_AUC)];
