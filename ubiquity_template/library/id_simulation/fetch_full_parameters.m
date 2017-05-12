function [parameters_full] = fetch_full_parameters(parameters_subset, cfg);
%
%  function [parameters_full] = fetch_full_parameters(parameters_subset, cfg);
%
% This function is used to build a full parameter set from a subset, and is
% normally used during parameter estimation in the observation details
% function when the entire parameter vector is needed to simulate the system.
%
% The function select_set pulls out a parameter set and can optionally select
% only a subset for estimation:
%
%    cfg=system_select_set(cfg, 'default', {'Vp', 'CL'});
%
% The default values of this subset can be accessed in the following way:
%
%    parameters_subset = cfg.estimation.parameters.guess 
%
% The estimation routines will work with this reduced parameter set, but to
% run simulations the full set is needed. The full values can be retrieved 
% using the following: 
%
% parameters_full = fetch_full_parameters(parameters_subset, cfg);
%


parameters_full = cfg.parameters.values;

parameters_full(cfg.estimation.to_estimate) = parameters_subset;
