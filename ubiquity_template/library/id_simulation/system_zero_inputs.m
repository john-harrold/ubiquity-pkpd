function [cfg]=system_zero_inputs(cfg, bolus, infusion_rates)
% function [cfg]=system_zero_inputs(cfg, bolus, rates)
% 
% To set all of the system inputs to zero:
%  cfg=system_zero_inputs(cfg)
%
% To only bolus inputs to zero
%  cfg=system_zero_inputs(cfg, true, false)
%
% To only infusion rates to zero
%  cfg=system_zero_inputs(cfg, false, true)


if(isfield(cfg.options, 'inputs'))
if(~exist('bolus', 'var'))
 bolus = true;
end

if(~exist('infusion_rates', 'var'))
 infusion_rates = true;
end


if(bolus)
if(isfield(cfg.options.inputs, 'bolus'))
if(isfield(cfg.options.inputs.bolus, 'species'))
species = fieldnames(cfg.options.inputs.bolus.species);

cfg.options.inputs.bolus.times.values = 0;
for(spidx = 1:length(species))
 eval(sprintf('cfg.options.inputs.bolus.species.%s.values = 0;', species {spidx}))
end

end
end
end



if(infusion_rates)
if(isfield(cfg.options.inputs, 'infusion_rates'))
rates   = cfg.options.inputs.infusion_rate_names;

for(ridx = 1:length(rates))
  eval(sprintf('cfg.options.inputs.infusion_rates.Rname.times.values  = 0;', rates{ridx}));
  eval(sprintf('cfg.options.inputs.infusion_rates.Rname.levels.values = 0;', rates{ridx}));
end

end
end
end
