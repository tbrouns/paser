function psr_sst_plot_waveforms(spikes, clustID, parameters, method)

% Based on 'plot_waveforms' from UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010

if (nargin < 4); method = 'heatmap'; end

% Parse inputs

spikeIDs  = ismember(spikes.assigns, clustID);
waveforms = spikes.waveforms(spikeIDs,:,:);
waveforms = psr_int16_to_single(waveforms,parameters);
Fs        = spikes.Fs;
ystep     = parameters.display.waveform_ystep;
nSpikes   = size(waveforms,1);

if (isempty(waveforms)); return; end

% Calculate average waveform maxima

med  = median(waveforms,1);
ymin = min(med(:));
ymax = max(med(:));
p2p  = ymax - ymin;
ymax = ymax + 0.5 * p2p;
ymin = ymin - 0.5 * p2p;
ymax = max(abs([ymin,ymax]));
if (ymin < ymax); ylims = [-ymax,ymax];
else,             ylims = [ -inf, inf]; % auto
end

clims = [0, parameters.display.waveform_clims];

% Plot waveforms

switch method
    
    case 'heatmap'
        
        cmap = parameters.display.cmap;
        [n,x,y] = histxt(waveforms(:,:),ystep);
        x = 1000 * (x - 1) / Fs; % in milliseconds
        imagesc(x,y,n,clims);
        set(gca,'Color', cmap(1,:)); % Set background color
        colormap(gca,cmap);
        caxis(clims);
        set(gca,'YDir','normal')
                
    case 'line'
        
        x = 1000 * ((1:size(waveforms(:,:),2)) - 1) / Fs;
        plot(x,waveforms(:,:),'Color',[0 0 0]);
        
end

% Make vertical lines to separate waveform into its channels
nSamples = size(waveforms,2);
nChans   = size(waveforms,3);
if nChans > 1
    l = zeros(nChans-1,1);
    for j = 1:nChans-1
        l(j) = line(1000 * (nSamples * j * [1 1] - 0.5) / Fs, 10 * ylims);
    end
    set(l,'Color',[0 0 0],'LineWidth',1.5) % electrode dividers
end

axis([x(1) x(end) ylims])

% Label axes

titleSubString = ['$\bf{Cluster \ \# ' num2str(clustID) ', \ N_{spikes} = ' num2str(nSpikes)];
if (parameters.display.metrics && ~psr_isempty_field(spikes,'spikes.clusters.metrics.rpv'))
    clustIDs = [spikes.clusters.metrics.id];
    I = find(clustIDs == clustID,1);
    fRPV = spikes.clusters.metrics(I).rpv;
    fRPV = round(100 * (100 * fRPV)) / 100;
    titleString = [titleSubString ' \ (' num2str(fRPV) ' \% \ RPVs)}$'];
else
    titleString = [titleSubString '}$'];
end

title(titleString,                 'Interpreter','Latex');
ylabel('$\bf{Voltage \ [\mu V]}$', 'Interpreter','Latex');
xlabel('$\bf{Time \ [ms]}$',       'Interpreter','Latex');
set(gca,'TickLabelInterpreter','Latex');

end

function [counts,t_inds,y_inds] = histxt(y,step)

% Based on 'histxt' from UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010

[nSpikes,nSamples] = size(y);

% Scale the data

yMin   = min(y(:));
yMax   = max(y(:));
yRange = yMax - yMin;
nSteps = round(yRange / step);

if (yRange == 0); y = repmat(nSteps / 2, size(y));
else,             y = (nSteps / yRange) .* (y - yMin);
end

y = round(y);

% Make bin centers/column indices

y_inds = linspace(yMin,yMax,nSteps);
t_inds = 1:nSamples;

% Histogram

if (size(y,1) == 1); y = repmat(y,2,1); nSpikes = 2; end % fix for "hist" function 

counts = hist(y,1:nSteps) / nSpikes;

end
