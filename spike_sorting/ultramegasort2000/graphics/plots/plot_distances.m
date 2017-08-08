function [x,y,z,dof] = plot_distances(spikes, show, method, display)
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% plot_distances - Histogram the distance of waveforms from cluster center.
%
% Usage:
%    [z,dof] = plot_distances( spikes, show, method );
%
% Description:
%    Plots histogram of Mahalanobis distance (z-value) from each waveform
% to the mean for a specified set of spike events.  Also plotted is the
% theoretical curve for a Gaussian distribution with the same number of
% degrees of freedom, as calculated from the chi2 distribution. This
% function is intended to be used as part of outlier removal. If the distribution
% of data points at large z-values is inconsistent with the Chi2 distribution,
% the user may want to remove these events as outliers.
%
% Input:
%  spikes - a spikes object
%
% Optional Inputs:
%   show          - array describing which events to show in plot
%                 - see get_spike_indices.m, (default = 'all')
%   method        - what covariance matrix to use when calculating distance
%                 - 1 => estimate covariance from cluster
%                 - 2 => estimate covariance from background noise
%                 - default value is in spikes.params.display.default_outlier_method
%
% Output:
%    z             - distance value for each specified event
%    dof           - degrees of freedom used in chi2 calculation
%

% parameter checking
if ~isfield(spikes,'waveforms'), error('No waveforms found in spikes object.'); end
if nargin < 2, show = 1:size(spikes.waveforms,1); end
if nargin < 3, method = spikes.params.display.default_outlier_method; end
if nargin < 4, display = 1; end

% which spikes are we showing?
show           = get_spike_indices(spikes, show);
data.waveforms = spikes.waveforms(show,:);
data.noise_cov = spikes.info.detect.cov;

% get z-scores
if method == 1
    data.data_cov  = cov(data.waveforms);
    [z,dof] = get_zvalues(data.waveforms, data.data_cov);
elseif method == 2
    [z,dof] = get_zvalues(data.waveforms, data.noise_cov);
end

% histogram binning
x     = 0:max(z);
[n,x] = histcounts(z,x);

[maxN,maxNi] = max(n);
maxNi        = maxNi + find(n(maxNi:end) > 0.01 * maxN, 1, 'last');
x            = 0:x(maxNi);

% calculate theoretical values

y = chi2pdf(x,dof);
y = y * length(z) * (x(2) - x(1));

% plotting

if (display)

    % initialize axes
    cla reset
    set(gca,'UserData', data.waveforms);
    
    histogram(z,x);
    xlim([0 200]);
    xlabel('Z-score'); 
    ylabel('Count');

    hndl = findobj(gca,'Type','patch');
    set(hndl,'FaceColor',[0 0 1])

    l = line(x,y);
    set(l,'Color',[0 1 0],'LineWidth',1.5)

end

end
