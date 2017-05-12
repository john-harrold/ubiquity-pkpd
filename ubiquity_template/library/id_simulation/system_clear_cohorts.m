function [cfg]=system_clear_chorts(cfg, bolus, infusion_rates)
% function [cfg]=system_zero_inputs(cfg, bolus, rates)
% 
% To clear all of the cohorts defined in the system:
%  cfg=system_clear_cohorts(cfg)

cfg.cohorts = struct();
