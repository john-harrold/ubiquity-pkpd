function [pval] = system_fetch_parameter(cfg, parameters, pname) 
% function [pval] = system_fetch_parameter(cfg, parameters, pname) 
%
% Retrieve the value of pname from the parameters vector.

pval  = -1;
isgood      = true;


if(length(parameters) ~= length(cfg.parameters.values))
  vp(cfg, sprintf('Error: Length of parameters >%d< dose not match system parameters length', length(parameters)));
  vp(cfg, sprintf('       system parameters length >%d<', length(cfg.parameters.values)));
  isgood      = false;
end

if(~isfield(cfg.options.mi.parameters, pname))
  vp(cfg, sprintf('Error: parameter name >%s< was not found. ', pname));
  isgood      = false;
end

if(isgood)
   %parameters = cfg.parameters.values;
   pidx = getfield(cfg.options.mi.parameters, pname);
   pval = parameters(pidx);
else
  vp(cfg, 'system_fetch_parameter()');
  vp(cfg, 'There was an error and the parameter was fetched.');
end
