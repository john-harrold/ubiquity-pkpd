function []=system_log_tryerror(cfg, TRYERROR)
% function []=system_log_tryerror(cfg, TRYERROR)
% 
%  Used to take the error from a catch/try/end and add it to the error log
%

vp(cfg, '---------------------------------');
vp(cfg, 'Error with catch/try/end');
vp(cfg, 'Message:');
vp(cfg, '--------');
vp(cfg, sprintf('   %s', TRYERROR.message));
vp(cfg, 'Stack:');
vp(cfg, '------');

stacknames = fieldnames(TRYERROR.stack);
for(sidx = 1:length(stacknames))

  sname  = stacknames{sidx};
  svalue = var2string_gen(getfield(TRYERROR.stack, stacknames{sidx}));

  vp(cfg, sprintf('   %s: %s', sname, svalue));

end
vp(cfg, '---------------------------------');
