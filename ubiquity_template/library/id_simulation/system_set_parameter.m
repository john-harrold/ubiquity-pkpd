function [parameters] = system_set_parameter(cfg, parameters, pname, value)
% [parameters] = system_set_parameter(cfg, parameters, pname, value)
%
%  parameters = full parameter vector  obtained using the following:
%
%  pname = name of the parameter to set
%  value = value of the parameter
%
%
%  To set the parameter Vc to a value of 3, the following would be used:
%
%  parameters = system_fetch_parameters(cfg);
%  parameters = system_set_parameter(cfg, parameters, 'Vc', 3);


if(sum(strcmp(cfg.parameters.names, pname)))
  parameters(getfield(cfg.options.mi.parameters, pname)) = value;
else
  vp(cfg, 'system_set_parameter()');
  vp(cfg, sprintf('parameter name (%s) not found', pname));
end
