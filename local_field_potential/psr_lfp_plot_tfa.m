function h = psr_lfp_plot_tfa(timefreq,parameters)

% PSR_LFP_PLOT_TFA - Plot spectrogram
%
% Syntax:  h = psr_lfp_plot_tfa(timefreq,parameters)
%
% Inputs:
%    timefreq   - Output from PSR_LFP_TFA function
%    parameters - See PSR_PARAMETERS_ANALYSIS
%
% Outputs:
%    h - Handle for spectrogram plot
%
% See also: PSR_LFP_TFA

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

if (isempty(parameters.analysis.tfa.base.window))
    powtype = lower(parameters.analysis.tfa.plot.powtype);
    switch powtype; case 'decibel'; timefreq.powspctrm = 10 * log10(timefreq.powspctrm); end
else
    powtype = lower(parameters.analysis.tfa.base.type);
end
yLabelStr = [upper(powtype(1)) powtype(2:end) ' \ power'];

% Check data

nDims = ndims(timefreq.powspctrm);
nTimes = length(timefreq.time);
if (nTimes ~= size(timefreq.powspctrm,nDims))
    if     (nDims == 3); timefreq.powspctrm = timefreq.powspctrm(:,:,  1:nTimes);
    elseif (nDims == 4); timefreq.powspctrm = timefreq.powspctrm(:,:,:,1:nTimes);
    end
end

% Spectrogram (Power vs. Frequency vs. Time)

% Remove artifacts field
timefreq = psr_remove_field(timefreq,'missing');
timefreq = psr_remove_field(timefreq,'artifacts');

if (all(isnan(timefreq.powspctrm(:)))); return; end

% Plot spectrogram
cfg             = [];
cfg.colormap    = parameters.analysis.tfa.plot.colormap;
cfg.interactive = 'no';
cfg.fontsize    = 11; % Default font-size
ft_singleplotTFR(cfg,timefreq); h = gca;
xlabel('$\bf{Time \ [s]}$',      'Interpreter','Latex');
ylabel('$\bf{Frequency \ [Hz]}$','Interpreter','Latex');

c = colorbar;
ylabel(c,['$\bf{' yLabelStr '}$'],'Interpreter','Latex');
set(gca,'TickLabelInterpreter','Latex');

tmin = min(timefreq.time);
tmax = max(timefreq.time);
xlim([tmin tmax]);

title('');

cmap = colormap;
set(gca,'Color',cmap(round(0.5 * size(cmap,1)),:));

if (~isempty_field(parameters,'parameters.analysis.tfa.plot.tlim'));  xlim(parameters.analysis.tfa.plot.tlim); end
if (~isempty_field(parameters,'parameters.analysis.tfa.plot.flim'));  ylim(parameters.analysis.tfa.plot.flim); end
if (~isempty_field(parameters,'parameters.analysis.tfa.plot.plim')); caxis(parameters.analysis.tfa.plot.plim); end

end