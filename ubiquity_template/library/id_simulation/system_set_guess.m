function [cfg]=system_set_guess(cfg, pname, value, lb, ub)
% [cfg] = system_set_guess(cfg,  pname, value, lb, ub)
%
%  pname = name of the parameter to set
%  value = value of the guess for the parameter
%  lb = lower bound (optional)
%  ub = upper bound (optional)
%
%  To set the initial guess for the parameter Vc to a value of 3, the
%  following would be used:
%
%  cfg = system_set_guess(cfg, 'Vc', 3);
%
%  To specify the guess and overwrite the upper bound on Vc and set it to 5
%
%  cfg = system_set_guess(cfg, 'Vc', 3, [], 5);
%

isgood = true;

if(sum(strcmp(cfg.parameters.names, pname)))
  if(isfield(cfg.estimation.mi, pname))
    cfg.estimation.parameters.guess(getfield(cfg.estimation.mi, pname)) = value;
  else
    vp(cfg, sprintf('parameter name (%s) was not selected for estimation', pname));
    vp(cfg, sprintf('see help for system_select_set '));

  end

  if(exist('lb', 'var'))
    if(~isempty(lb))
      cfg.estimation.parameters.lower_bound(getfield(cfg.estimation.mi, pname)) = lb;
    end
  end
  if(exist('ub', 'var'))
    if(~isempty(ub))
      cfg.estimation.parameters.upper_bound(getfield(cfg.estimation.mi, pname)) = ub;
    end
  end

else
  isgood = false;
  vp(cfg, sprintf('parameter name (%s) not found', pname));
end


if(isgood == false)
  vp(cfg, 'system_set_guess()');
end
