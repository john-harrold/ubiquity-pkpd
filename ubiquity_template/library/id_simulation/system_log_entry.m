function [] = system_log_entry(cfg, entry)
%
% % Initialize the log file:
% % transient/ubiquity_log.txt
% cfg = system_log_init(cfg);
%
% Add a log entry
% system_log_entry(cfg, 'this is a log entry');
%

% If logging is disabled we don't do anything 
if(strcmp(cfg.options.logging.enabled, 'yes'))
  % If the log file doesn't exist we initialize it
  if(~exist(cfg.options.logging.file, 'file'))
   system_log_init(cfg);
  end
  
  % If the timestamp is enabled we prepend it to the
  % log message
  if(strcmp(cfg.options.logging.timestamp, 'yes'))
    entry = sprintf('%s %s', datestr(now, cfg.options.logging.ts_str ), entry);
  end

  % Now we dump it to the log file:
  FID = fopen(cfg.options.logging.file, 'a');
  fprintf(FID, '%s\n', entry);
  fclose(FID);
end
