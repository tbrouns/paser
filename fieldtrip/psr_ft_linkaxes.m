function psr_ft_linkaxes(figHandle,axs,titles,climType)

% PSR_FT_LINKAXES - Links axes for FieldTrip figure
%
% Syntax:  psr_ft_linkaxes(figHandle,axs,titles,climType)
%
% Inputs:
%    figHandle - Handle to the figure we want to link the axes for
%    axs       - Array of axis handles for each subplot
%    titles    - Cell array of titles  for each subplot
%    climType  - Axis limits for the colorbar, can be set to:
%                'both'  : symmetric axis for colorbar
%                'other' : colorbar starts at zero
%
% Example: 
%    fig = figure;
%    [jpsth,jpsthShuff] = psr_ft_jpsth(psth,parameters);
%    axs(1,1) = subaxis(1,2,1); psr_ft_plot_jpsth(jpsth,     parameters,[1 2]); axs(1,2) = gca;
%    axs(2,1) = subaxis(1,2,2); psr_ft_plot_jpsth(jpsthShuff,parameters,[1 2]); axs(2,2) = gca;
%    psr_ft_linkaxes(fig,axs,{'Shift-predicted','Shift-predictor'},'both');
%
% See also: PSR_FT_PLOT_JPSTH,  PSR_FT_PLOT_ISIH

% PASER: Processing and Analysis Schemes for Extracellular Recordings
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept.
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

if (nargin < 4); climType = []; end

set(0,'CurrentFigure',figHandle);

nAxes = size(axs,1); % Number of subplots

% Colorbar limits
climits = zeros(nAxes,2);
for i = 1:nAxes; climits(i,:) = axs(i,1).CLim; end
climits = mean(climits);
climits = max(abs(climits)); 
if (strcmp(climType,'both')); climits = [-climits climits];
else,                         climits = [       0 climits];
end

for i = 1:nAxes
    set(figHandle,'currentaxes',axs(i,1));
    caxis(climits); 
    set(axs(i,2),'CLim',climits);
    if (length(titles) >= i); title(titles{i}); end
end

axesAll = findall(figHandle,'type','axes');
n = length(axesAll) / nAxes;
if (mod(n,1) == 0)
    for i = 1:n
        for j = 1:nAxes
            k = i + n * (j - 1);
            axesTemp(j) = axesAll(k); 
        end
        linkaxes(axesTemp);
    end
end
                
end