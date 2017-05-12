function [h] = mysubplot(rows, cols, idx)
% function [h] = mysubplot(rows, cols, idx)
%
% Has the same syntax as the subplot command but the gap between subplots is
% decreased. 
h = subplot(rows, cols, idx);
p = get(h, 'pos');

% making room for axis and tick labels
p(1) = p(1) - 0.035;
p(2) = p(2) - 0.035;
% shrinking the distance between the plots
p(3) = p(3) + 0.07;
p(4) = p(4) + 0.07;
set(h, 'pos', p);
