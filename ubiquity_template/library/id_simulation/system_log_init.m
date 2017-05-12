function [cfg]=system_log_init(cfg)


try
  FID = fopen(cfg.options.logging.file, 'w');
  fclose(FID);
  cfg.options.logging.enabled = 'yes';
  system_log_entry(cfg, 'Ubiquity log init - Matlab');
catch
  disp(sprintf(' system_log_init()'));
  disp(sprintf(' Unable to create file >%s<', cfg.options.logging.file));
  disp(sprintf(' Logging disabled', cfg.options.logging.file));
  cfg.options.logging.enabled = 'no';
end
