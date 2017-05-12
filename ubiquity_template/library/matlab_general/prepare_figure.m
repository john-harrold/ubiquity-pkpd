function []=prepare_figure(purpose)
% function []=prepare_figure(purpose)
%
%  Prepares the current axis of a figure for either printing or presentation.
%  This is defined by the 'purpose'. 
%
%  # example: preparing a figure for printing
%  prepare_figure('print')
%
%  # example: preparing a figure for a presentation
%  prepare_figure('present')
%
%  # example: preparing a figure for a poster       
%  prepare_figure('poster')
%

if(strcmp(purpose,'present'))
   h=gca;
   set(h,'Fontsize',14);
   set(h,'linewidth',1.5);
   try
   set(get(gca,'children'),'linewidth',1.5);
   end
elseif(strcmp(purpose,'poster'))
   h=gca;
   set(h,'Fontsize',16);
   set(h,'linewidth',1.5);
   try
   set(get(gca,'children'),'linewidth',1.5);
   end
elseif(strcmp(purpose,'print'))
   h=gca;
   set(h,'Fontsize',12);
   set(h,'linewidth',1.1);
   try
   set(get(gca,'children'),'linewidth',1.);
   end
end
