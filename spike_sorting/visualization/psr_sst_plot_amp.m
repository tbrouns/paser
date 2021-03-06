function p = psr_sst_plot_amp(spikes,clustID,parameters)

% based on plot_detection_criterion from UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% plot_detection_criterion - histogram of detection criterion applied to
%                            each waveform
%
% Usage:
%     [thresh_val,m,stdev,missing] = plot_detection_criterion(spikes,clus)
%
% Description:
%    Plots a histogram of the detection metric applied to each waveform. A
% normalization is performed so that the threshold has a value of +/- 1
% (vertical black dotted line). The histogram is overlaid with a fitted
% Gaussian (red line) that is used to estimate the number of undetected
% spikes. This percentage is shown in the axes title.
%
%  See psr_sst_amp_gaussfit.m for more information.
%
% Input:
%   spikes     - a spike structure
%   clus       - cluster ID or array describing which events to show in plot
%              - see get_spike_indices.m
%
% Output:
%  p            - estimate of probability that a spike is missing because it didn't reach threshhold
%

% Grab data
[p,mu,stdev] = psr_sst_amp_gaussfit(spikes,clustID,parameters);
[x,n]        = psr_sst_amp_hist    (spikes,clustID,parameters,true);

% Determine global extreme
global_ex = max(abs(x));

% Histogram
bar(x,n,1.0,'FaceColor','k','EdgeColor','none','FaceAlpha',1.0);

if (spikes.info.detected)
    if (parameters.display.show_gaussfit) % Draw Gaussian fit
        a = linspace(0,global_ex,200);
        b = normpdf(a,mu,stdev);
        b = b * (max(n) / max(b));
        l = line(a,b); set(l,'Color',[1 0 0],'LineWidth',1.5)
    end
end

% Threshold line
l = line([1 1], get(gca,'YLim'));
set(l,'LineStyle','--','Color',[0 0 1],'LineWidth',2)

% Prettify axes
axis tight

titleString = [];
if (parameters.display.metrics)
    titleString = ['\bf{Sub-threshold \ spikes: \ ' num2str(p*100,'%2.1f') '\% }']; 
end

title(titleString,                    'Interpreter','Latex');
xlabel('\bf{Normalized \ amplitude}', 'Interpreter','Latex');
ylabel('\bf{No. \ of \ spikes}',      'Interpreter','Latex');
set(gca,'TickLabelInterpreter','Latex');

% Determine x limit

N = cumsum(n);
N = N / sum(n);
I = find(N > 0.99,1,'first');
global_ex = x(I);
xlim([0 global_ex])

end

