function spikes = psr_kst_convert2spikes(rez,data,parameters)

% PSR_KST_CONVERT2SPIKES - Convert KiloSort output to PASER data format.
%
% Syntax:  spikes = psr_kst_convert2spikes(rez,data,parameters)
%
% Inputs:
%    rez - Output of KiloSort spike sorting
%    data - Filtered time series of extracellular recording 
%    (number of probe channels x number of data points)
%    parameters - See README
%
% Outputs:
%    spikes - See README
%
% See also: PSR_SST_SORTING_KST

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2017

%------------- BEGIN CODE --------------

Fs = parameters.Fs;
window_samples = round(Fs * parameters.spikes.window_size / 1000);
samples_hwidth = round(0.5 * window_samples);
win = -samples_hwidth:samples_hwidth;

sLength = size(data,2);

% extract info from rez
spikeTimes     = rez.st3(:,1);
spikeTemplates = rez.st3(:,2);
if (size(rez.st3,2) >= 5)
    spikeClusters = 1+rez.st3(:,5);
else
    spikeClusters = spikeTemplates;
end

% get raw data around spiketimes

times = bsxfun(@plus,spikeTimes,win);
times(times < 1) = 1;
times(times > sLength) = sLength;
times = times';
times = times(:);
waves = data(:,times);
waves = permute(waves,[3 2 1]);
waves = reshape(waves,length(win),[],parameters.general.nelectrodes);
waves = permute(waves,[2 1 3]);

%% Save results

data = psr_single(data',parameters);

spikes             = [];
spikes.assigns     = int16(spikeClusters');
spikes.spiketimes  = single((spikeTimes - 1)' / Fs); % in secs
spikes.waveforms   = waves;
spikes.info.dur    = (sLength - 1) / Fs; % in secs
spikes.info.thresh = -parameters.spikes.thresh * (median(abs(data)) / 0.6745);
spikes.info.stds   = std(data);
rez = rmfield(rez,'st3'); % Remove now-redundant field
spikes.info.kst           = rez;

%------------- END OF CODE --------------