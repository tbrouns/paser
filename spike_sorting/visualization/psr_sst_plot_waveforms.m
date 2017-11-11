function psr_sst_plot_waveforms(spikes, show, parameters)

% Based on 'plot_waveforms' from UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010

cla reset;

clustID = show;

% Argument checking

if ~isfield(spikes,'waveforms'), error('No waveforms found in spikes object.'); end
if nargin < 2; show = 1:size(spikes.waveforms,1); end

% Which spikes are we showing?

show = get_spike_indices(spikes, show);

% Parse inputs

spiketimes = sort(spikes.spiketimes(show));
waveforms  = spikes.waveforms(show,:,:);
Fs         = spikes.Fs;
wstep      = spikes.params.display.waveform_step;

nRPV = sum(diff(spiketimes) <= 0.001 * parameters.spikes.ref_period);

% Calculate average waveform maxima 

med   = median(waveforms,1);
ymin  = min(med(:));
ymax  = max(med(:));
p2p   = ymax - ymin; 
ymax  = ymax + 0.5 * p2p;
ymin  = ymin - 0.5 * p2p;
ylims = [ymin, ymax];

barcolor = [0.5 0 0];

% Plot waveforms

cmap    = spikes.params.display.cmap;
[n,x,y] = histxt(waveforms(:,:),wstep);
x = 1000 * (x - 1) / Fs; % in milliseconds
imagesc(x,y,n);
colormap(cmap);
set(gca,'Color', cmap(1,:));
axis([x(1) x(end) ylims])
set(gca,'YDir','normal')

% Make vertical lines to separate waveform into its channels
num_channels = size(waveforms,3);
num_samples  = size(waveforms,2);
if num_channels > 1
    l = zeros(num_channels-1,1);
    for j = 1:num_channels-1
        l(j) = line(1 + num_samples * j * [1 1], ylims);
    end
    set(l,'Color',barcolor,'LineWidth',1.5) % electrode dividers
end

% Label axes

N = size(waveforms,1);
fRPV = round(100 * (100 * nRPV / N)) / 100;
ystr = ['$\bf{Cluster \ \# ' num2str(clustID) ', \ N = ' num2str(N) ' \ (' num2str(fRPV) ' \% \ RPVs)}$'];
title(ystr,                        'Interpreter','Latex');
ylabel('$\bf{Voltage \ [\mu V]}$', 'Interpreter','Latex');
xlabel('$\bf{Time \ [ms]}$',       'Interpreter','Latex');
set(gca,'TickLabelInterpreter','Latex');

end

function [counts,t_inds,x_inds] = histxt(x,step)

% Based on 'histxt' from UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010

[~,T] = size(x);

% Scale the data

oldmin   = min(x(:));  
oldmax   = max(x(:));
oldrange = oldmax - oldmin;
D        = round(oldrange / step);

if (oldrange == 0); x = repmat(D / 2, size(x));
else,               x = (D ./ oldrange) .* (x - oldmin);
end

x = round(x);

% Make bin centers/column indices 

x_inds = linspace(oldmin,oldmax,D);
t_inds = 1:T;

% Histogram

counts = hist(x,1:D);

end
