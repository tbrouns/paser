function psr_plot_stacking(h)

% PSR_PLOT_STACKING - Stack two subplots on top of each other
%
% Syntax:  psr_plot_stacking(h)
%
% References:
% [1] https://stackoverflow.com/a/11758610 (by Gunther Struyf)
% 
% Inputs:
%    h - Two-element vector of axis handles
%
% Example: 
%   ax1 = subplot(1,1,1); plot(x1,y1);
%   ax2 = subplot(2,1,2); plot(x2,y2);
%   psr_plot_stacking([ax1 ax2]);

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

% 
% h(1): top plot
% h(2): bottom plot

linkaxes([h(1),h(2)],'x')

set(h(1),'xticklabel',[]);        
pos    = get(h,'position');
bottom = pos{2}(2);
top    = pos{1}(2) + pos{1}(4);
plotspace = top - bottom;
pos{2}(4) = plotspace / 2;
pos{1}(4) = plotspace / 2;
pos{1}(2) = bottom + plotspace / 2;

set(h(1),'position',pos{1});
set(h(2),'position',pos{2});

% Fix alignment in case of colorbar
pos1 = get(h(1),'Position');
pos2 = get(h(2),'Position');
pos2(3) = pos1(3);
set(h(2),'Position',pos2)

% Move y-labels to the right
set(h(2),'YAxisLocation','right');
        
end