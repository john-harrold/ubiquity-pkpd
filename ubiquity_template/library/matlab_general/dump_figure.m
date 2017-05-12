function []=dump_figure(filename)
%
% Places the currently opened figure in 'filename'.
%
%    -> If you're on a mac it will create filename.pdf and png
%    -> If you're on windows it will create filename.tif, jpg, and png
%


if(ismac)
  eval(sprintf('export_fig %s -transparent -pdf',  filename));
  try
    eval(sprintf('export_fig %s -transparent -png -jpeg -m4  -painters', filename));
  end
elseif(ispc) 
  eval(sprintf('export_fig %s -transparent -tiff -jpeg -m4', filename));
end

% Saving the matlab fig file for editing later
 if(verLessThan('matlab', '8.2'))
   saveas(gcf, sprintf('%s.fig', filename))
 else
   savefig(gcf, sprintf('%s.fig', filename))
 end
