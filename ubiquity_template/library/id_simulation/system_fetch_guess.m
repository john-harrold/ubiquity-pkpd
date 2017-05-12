function [guess]      = system_fetch_parameters(cfg, PNAME) 
guess  = cfg.estimation.parameters.guess;

if(exist('PNAME','var'))
  guess = guess(getfield(cfg.estimation.mi, PNAME))
end


