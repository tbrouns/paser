function [ax1,ax2] = psr_sst_plot_stability(spikes,clustID,parameters,metadata)
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% plot_stability - Plot amplitude and firing rate of a set of spikes over time
%
% Usage:
%      [ax1,ax2] = plot_stability( spikes, clus)
%
% Description:
%    Creates two axes in the same position to simultaneously plot two
% aspects of the stability of a given set of spikes.  In the upper half
% and using the right y-axis is a scatter plot of the peak-to-peak amplitude
% of waveforms taken from the full duration of the recording.  For efficiency
% a maximum number of randomly drawn data points is used as set by
% parameters.display.max_scatter.  In the lower half and using the left
% y-axis is a histogram of the firing rate as a function of time.  The
% unwrapped time is used (see unwrap_time.m).  The size of the bins used
% is set by parameters.display.stability_bin_size.
%
% Input:
%   spikes - a spike structure
%
% Optional input:
%   show          - array describing which events to show in plot
%if nargin < 1, spikes = [];    end

% Output:
%  ax1 - handle to axes the firing rate over time
%  ax2 - handle to axes that display the amplitude scatter over time
%

% Check input

if nargin < 3 || isempty_field(metadata,'metadata.trialonset'); blockOnsets = [];
else,                                                               blockOnsets = metadata.trialonset;
end

% Select waveforms

spikeIDs   = find(ismember(spikes.assigns, clustID));
waves      = spikes.waveforms(spikeIDs,:);
spiketimes = sort(spikes.spiketimes(spikeIDs));

% Grab artifacts

noisetimes  = [];
if (~isempty(blockOnsets))
    if (~isempty_field(spikes,'spikes.info.artifacts'))
        nBlocks = length(spikes.info.artifacts.lfp);
        for iBlock = 1:nBlocks
            artifactsBlock = spikes.info.artifacts.lfp{iBlock};
            fields = fieldnames(artifactsBlock);
            nFields = length(fields);
            for iField = 1:nFields
                artifacts  = artifactsBlock.(fields{iField});
                noisetimes = [noisetimes;artifacts + blockOnsets(iBlock)]; %#ok
            end
        end
        D = noisetimes(:,2) - noisetimes(:,1);
        [~,I] = sortrows(D,'descend');
        noisetimes = noisetimes(I,:);
        
        % Limit number of artifacts that we show
        nmax = parameters.display.max_artifacts;
        if (size(noisetimes,1) > nmax)
            noisetimes = noisetimes(1:nmax,:);
        end
    end
end

% Initialize axes
cla
ax1 = gca;

% We are ploting two axes on top of each other.
% The first axes will show the firing rate as a function of time

% Bin data
tlims    = [0 sum(spikes.info.dur)];
num_bins = round(diff(tlims) / parameters.display.stability_bin_size);
edges    = linspace(tlims(1),tlims(2),num_bins+1);
n        = histc(spiketimes,edges);
n(end)   = [];
vals     = n / mean(diff(edges));
edges    = edges(1:end-1) + mean(edges(1:2));

% Show histogram
hold on; bar(edges,vals,1.0,'FaceColor','k','EdgeColor','none','FaceAlpha',1.0); hold off
ylim_ax1 = max(get(gca,'YLim'));
set(gca, 'XLim', tlims,'YLim',[0 2 * ylim_ax1])
yticks = get(gca,'YTick'); set(gca,'YTick',yticks(yticks <= max(yticks)/2))

set(gca,'TickLabelInterpreter','Latex');
xlabel('\bf{Time \ (sec)}',         'Interpreter','Latex')
ylabel('\bf{Firing \ rate \ (Hz)}', 'Interpreter','Latex')

ylabh = get(gca,'ylabel');
set(ylabh,'Units','normalized');
set(ylabh,'position', get(ylabh,'position') - [0.01 0.25 0.00]);

% Second pass is scatter plot of amplitude versus time

% Get peak-to-peak amplitudes over time
amp = range(waves');

% Only plot a random subset
if isequal(parameters.display.max_scatter, 'all')
    ind = 1:length(amp);
else
    choice = randperm(length(amp));
    max_pos = min(length(amp), parameters.display.max_scatter);
    ind = choice(1:max_pos);
end

% Determine y-limit

med   = median(waves,1);
ymin  = min(med(:));
ymax  = max(med(:));
p2p   = ymax - ymin; 
if (p2p > 0); ymax = 2 * p2p;
else,         ymax = 1;
end

% Prettify axes
ax2 = axes('Parent',get(ax1,'Parent'),'Unit',get(ax1,'Unit'),'Position',get(ax1,'Position'));

% We also set Z coordinate to ensure that scatter plot is on the foreground
hold on; l = scatter3(spiketimes(ind),amp(ind),ones(size(spiketimes(ind)))); hold off; 
l1 = line(tlims, [0 0],'LineWidth',2);
set(l1,'Color',[0 0 0]);
hold off
ylim_ax2 = ymax;
set(l,'Marker','.','MarkerEdgeColor',[.3 .5 .3])
set(ax2,'Xlim',tlims,'YLim',ylim_ax2 * [-1 1],'XTickLabel',{},'YAxisLocation','right');
yticks = get(ax2,'YTick'); set(ax2,'YTick',yticks(yticks>=0))
set(ax2,'Color','none')

% Move y-label

set(gca,'TickLabelInterpreter','Latex');
ylabel('\bf{Peak-to-peak \ amplitude}','Interpreter','Latex')

ylabh = get(gca,'ylabel');
set(ylabh,'Units','normalized');
set(get(gca,'YLabel'),'Rotation',270)
set(ylabh,'position', get(ylabh,'position') + [0.03 0.25 0]);

% Draw trial boundaries

if (~isempty(blockOnsets))
    nBlocks = length(blockOnsets);
    hold on
    
    axes(ax1);
    for iBlock = 2:nBlocks
        t = blockOnsets(iBlock);
        line([t t],[0 1] * ylim_ax1,'Color','k','LineStyle','--','LineWidth',1.5);
    end
    
    axes(ax2);
    for iBlock = 2:nBlocks
        t = blockOnsets(iBlock);
        line([t t],[0 1] * ylim_ax2,'Color','k','LineStyle','--','LineWidth',1.5);
    end
    
    hold off
end

% Show LFP artifacts

hold on;
for i = 1:size(noisetimes,1)
    F = area(noisetimes(i,:),[1 1] * ylim_ax2,'EdgeColor','none','FaceColor','r');
    alpha(F, 0.2)
end
hold off

% Link properties of the 2 axes together

linkaxes([ax1 ax2],'x');
linkprop([ax1 ax2],'Position');
linkprop([ax1 ax2],'Unit');