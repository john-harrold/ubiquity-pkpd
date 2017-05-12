function [mc]=fetch_color_codes();
% [mc]=fetch_color_codes();
%
%
%  % mc is a data structure with elements that contain color codes in the
%  % pattern:
%
%  mc.weight_color
%
%  % Where weight can be:
%    -faint
%    -light
%    -medium
%    -dark
%
%  % And color can be:
%    -red
%    -orange
%    -yellow
%    -green
%    -bluegreen
%    -blue
%    -purple
%
%
%  % To get a light blue you can use mc.light_blue
%
%
%  % To see the different colors try:
%
%  close all;
%  provision_workspace;
%  mc = fetch_color_codes;
%  
%  colors = {'red' 'orange', 'yellow', 'green' 'bluegreen', 'blue', 'purple'  };
%  weights = {'faint', 'light', 'medium', 'dark'};
%  figure(1);
%  hold on;
%  for cidx= 1:length(colors)
%    color = colors{cidx};
%  
%    for widx = 1:length(weights)
%      weight  = weights{widx};
%  
%      pc      = getfield(mc, sprintf('%s_%s', weight, color));
%      
%      msize   = 20;
%      
%      plot(widx, cidx, 'o', 'color', pc,   'markerfacecolor', pc,  'markersize', msize);
%    end
%  end
%  axis tight
%  
%  
%  xlim([min(xlim())-.5, max(xlim())+.5]);
%  ylim([min(ylim())-.5, max(ylim())+.5]);
%  
%  set(gca, 'xtick',      [1:length(weights)]);
%  set(gca, 'xticklabel', weights);
%  set(gca, 'ytick',      [1:length(colors)]);
%  set(gca, 'yticklabel', colors); 
%



mc.faint_green      = [209 255 210]./255; 
mc.light_green      = [173 222 190]./255; 
mc.medium_green     = [ 84 209 128]./255; 
mc.dark_green       = [  0 135  47]./255;
mc.green            = mc.medium_green;

mc.faint_purple     = [241 231 255]./255;
mc.light_purple     = [241 186 255]./255;
mc.medium_purple    = [220  79 255]./255;
mc.dark_purple      = [100   0 150]./255;
mc.purple           = mc.medium_purple;


mc.faint_yellow     = [255 255 205]./255;
mc.light_yellow     = [255 255 102]./255;
mc.medium_yellow    = [255 255   0]./255;
mc.dark_yellow      = [255 240   0]./255;
mc.yellow           = mc.medium_yellow;



mc.faint_blue       = [231 237 253]./255;
mc.light_blue       = [171 226 254]./255;
mc.medium_blue      = [ 99 121 219]./255;
mc.dark_blue        = [  7  31 138]./255;
mc.blue             = mc.medium_blue;

mc.faint_red        = [253 231 231]./255;
mc.light_red        = [255 166 166]./255;
mc.medium_red       = [255   0   0]./255;
mc.dark_red         = [188  37  96]./255;
mc.red              = mc.medium_red;
                    
mc.faint_orange     = [253 249 231]./255;
mc.light_orange     = [255 201 132]./255;
mc.dark_orange      = [214 111   0]./255;
mc.medium_orange    = [255 174   0]./255;
mc.orange           = mc.medium_orange;

mc.faint_bluegreen  = [231 253 242]./255;
mc.light_bluegreen  = [150 177 232]./255;
mc.medium_bluegreen = [ 46 129 232]./255;
mc.dark_bluegreen   = [  0 161 153]./255;
mc.bluegreen        = mc.medium_bluegreen;


mc.black            = [0,0,0];
mc.white            = [1,1,1];
mc.light_grey       = [255 166 166]./255;
mc.light_grey       = [ .70  .75  .71];
