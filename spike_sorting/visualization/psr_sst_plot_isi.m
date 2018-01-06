function psr_sst_plot_isi(spikes, clustID, parameters)

% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% psr_sst_plot_isi - plot histogram of ISI distribution for a cluster
%
% Usage:
%       psr_sst_plot_isi(spikes, show, show_isi)
%
% Description:
%   Plots a histogram either representing the inter-spike interval (ISI)
% distribution for the selected spike events or the autocorrelation function
% in Hertz.  Also plotted is a gray area representing the "shadow"
% period during spike detection and a red area representing the user
% defined refractory period.
%
% The width of the histogram bins as well as the maximum time lag displayed
% are set by the following parameters:
%
%    parameters.display.isi_bin_size
%    parameters.display.max_isi_to_display
%    parameters.display.correlations_bin_size
%    parameters.display.max_autocorr_to_display
%
% The user can switch between displaying ISIs or autocorrelatoin by right-
% clicking the axes and selecting from a context menu.  The choice can
% also be imposed on all isi plots in the same figure.
%
% On the y-axis is listed the number of refractory period violationss (RPVs)
% along with an estimated contamination percentage and its 95% confidence
% interval under the assumption that contaminating spikes are independent
% events. See poisson_contamination.m for more details.
%
% Inputs:
%   spikes        - a spikes structure
%
% Optional inputs:
%   show          - array describing which events to show in plot
%                 - see get_spike_indices.m, (default = 'all')
%   show_isi      - 1 -> show ISI, 0 -> show autocorrelation
%                 - default is set by parameters.display.show_isi
%

% Set constants
show_isi        = parameters.display.show_isi; 
shadow          = 0.5 * parameters.spikes.window_size;
rp              = parameters.spikes.ref_period;
isi_maxlag      = parameters.display.max_isi_to_display;
autocorr_maxlag = parameters.display.max_autocorr_to_display;
corr_bin_size   = parameters.display.correlations_bin_size;
isi_bin_size    = parameters.display.isi_bin_size;

% Get the spiketimes
spikeIDs   = get_spike_indices(spikes, clustID);
spiketimes = sort(spikes.spiketimes(spikeIDs));

% ISI case
if show_isi
    maxlag = isi_maxlag;
    bins   = round(1000 * maxlag / isi_bin_size);
    
    % Make plot
    isis  = diff(spiketimes);
    isis  = isis(isis <= maxlag);
    [n,x] = hist(isis * 1000,linspace(0,1000 * maxlag,bins));
    ymax  = max(n) + 1;
    
    % Make patches to represent shadow and refractory period
    patch(  [0 shadow shadow 0], [0 0 ymax ymax], [0.5 0.5 0.5], 'EdgeColor', 'none');
    patch([shadow rp rp shadow], [0 0 ymax ymax], [1.0 0.0 0.0], 'EdgeColor', 'none');
    hold on; b2 = bar(x,n,1.0); hold off
    set(b2,'FaceColor',[0 0 0],'EdgeColor',[0 0 0])
    
    % Update axes
    set(gca,'YLim',[0 ymax],'XLim',[0 1000 * maxlag])
    xlabel('\bf{Interspike \ interval \ [ms]}', 'Interpreter','Latex');
    ylabel('\bf{No. \ of \ spikes}',            'Interpreter','Latex');
    set(gca,'TickLabelInterpreter','Latex');
else
    maxlag = autocorr_maxlag;
    
    % Calculate autocorrelation
    if length(spiketimes) > 1
        [cross,lags] = pxcorr(spiketimes, spiketimes, round(1000 / corr_bin_size), maxlag);
    else
        cross = 0;
        lags  = 0;
    end
    cross(lags == 0) = 0;
    
    % Place patches to represent shadow and refractory period
    ymax = max(cross) + 1;
    patch(  shadow * [-1 1 1 -1], [0 0 ymax ymax], [0.5 0.5 0.5], 'EdgeColor', 'none');
    patch( [shadow rp rp shadow], [0 0 ymax ymax], [1.0 0.0 0.0], 'EdgeColor', 'none');
    patch(-[shadow rp rp shadow], [0 0 ymax ymax], [1.0 0.0 0.0], 'EdgeColor', 'none');
    
    % Plot autocorrelation histogram
    hold on; bb = bar(1000 * lags,cross,1.0); hold off;
    set(bb,'FaceColor',[0 0 0],'EdgeColor',[0 0 0])
    
    % Set axes
    set(gca,'XLim', 1000 * maxlag * [-1 1]);
    set(gca,'YLim', [0 ymax]);
    xlabel('\bf{Time \ lag \ [ms]}', 'Interpreter','Latex');
    ylabel('\bf{Autocorrelation}',   'Interpreter','Latex');
    set(gca,'TickLabelInterpreter','Latex');
end

end