function []=dump_figure(filename)
%
% Places the currently opened figure in 'filename' as: pdf, png and fig
%

exportgraphics(gcf, sprintf("%s.pdf", filename))
exportgraphics(gcf, sprintf("%s.png", filename))
saveas(gcf, sprintf("%s.fig", filename))
% if(ismac)
%   eval(sprintf('export_fig %s -transparent -pdf',  filename));
%   try
%     eval(sprintf('export_fig %s -transparent -png -jpeg -m4  -painters', filename));
%   end
% elseif(ispc) 
%   eval(sprintf('export_fig %s -transparent -tiff -jpeg -m4', filename));
% end
% 
% % Saving the matlab fig file for editing later
%  if(verLessThan('matlab', '8.2'))
%    saveas(gcf, sprintf('%s.fig', filename))
%  else
%    savefig(gcf, sprintf('%s.fig', filename))
%  end
