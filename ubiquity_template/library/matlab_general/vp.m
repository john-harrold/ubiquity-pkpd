function []=vp(cfg, str)
% function []=vp(cfg, str)
% vp -- verbose print
%
% Print out the message contained in 'str' if the verbose option in cfg is
% set. Example:
%
%  cfg.options.verbose = 'yes';
%
%  vp(cfg, 'Hellow world');
%


system_log_entry(cfg, str);

if(isfield(cfg, 'options'))
if(isfield(cfg.options, 'verbose'))
if(strcmp('yes', cfg.options.verbose))
  disp(sprintf('#> %s',str));
end
end
end
