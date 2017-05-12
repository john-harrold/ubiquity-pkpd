function []=archive_estimation(name, cfg)
%
%  archive_estimation(name, cfg)
%
%  Archives the estimation results by moving the output files to the same file
%  names with 'name' prepended to them. This prevents them from being
%  overwritten in a different analysis script the following files are
%  archived:
%
%   output/monitor_estimation_progress.jpg 
%   output/monitor_estimation_progress.pdf
%   output/monitor_estimation_progress.png
%   output/parameters_all.csv             
%   output/parameters_est.csv             
%   output/report.txt                     
%
%  Example:
%   archive_estimation('mysoln', cfg)
%
%   Would rename the files above 
%   output/mysoln-monitor_estimation_progress.jpg 
%   output/mysoln-monitor_estimation_progress.pdf
%   output/mysoln-monitor_estimation_progress.png
%   output/mysoln-parameters_all.csv             
%   output/mysoln-parameters_est.csv             
%   output/mysoln-report.txt                     
%



%        old filename                                                 new filename
files = {sprintf('output%smonitor_estimation_progress.jpg',filesep),  sprintf('output%s%s-monitor_estimation_progress.jpg',filesep, name)
         sprintf('output%smonitor_estimation_progress.pdf',filesep),  sprintf('output%s%s-monitor_estimation_progress.pdf',filesep, name)
         sprintf('output%smonitor_estimation_progress.png',filesep),  sprintf('output%s%s-monitor_estimation_progress.png',filesep, name)
         sprintf('output%sparameters_all.csv'             ,filesep),  sprintf('output%s%s-parameters_all.csv'             ,filesep, name)
         sprintf('output%sparameters_est.csv'             ,filesep),  sprintf('output%s%s-parameters_est.csv'             ,filesep, name)
         sprintf('output%sreport.txt'                     ,filesep),  sprintf('output%s%s-report.txt'                     ,filesep, name)};

for fidx = 1:length(files(:,1))
  oldfile = files{fidx,1};
  newfile = files{fidx,2};

  if(exist(oldfile, 'file'))
    movefile(oldfile, newfile);
  else
    if(strcmp(cfg.options.verbose, 'yes'))
      vp(cfg, sprintf('unable to archive file: %s', oldfile));
    end
  end
  end
end
