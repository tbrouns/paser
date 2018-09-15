function timefreq = psr_lfp_baseline(timefreq,baseline,parameters)

% PSR_LFP_BASELINE - Baseline correction for power spectrum
% 
% Syntax:  data = psr_lfp_baseline(timefreq,baseline,parameters)
%
% Based on FieldTrip's FT_FREQBASELINE
%  
% Inputs:
%    timefreq    - Output structure from PSR_LFP_TFA
%    baseline   - Input baseline, either with the same size as
%                 "timefreq.powspctrm" or can be set to empty to calculate
%                 the baseline from "timefreq.powspctrm".
%    parameters - See PSR_PARAMETERS_ANALYSIS
%
% Outputs:
%    timefreq - Baseline corrected power spectrum
%
% See also: FT_FREQBASELINE

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

if (ndims(baseline) ~= ndims(timefreq.powspctrm) || all(size(baseline) ~= size(timefreq.powspctrm)))
    
    nTrials = size(timefreq.powspctrm,1);
    baseTrials = 1:nTrials;
    baseCycles = [];
    baseWindow = [];
       
    if (strcmp(parameters.analysis.tfa.base.window,'pre'))
        baseWindow = [timefreq.time(1) 0];
    elseif (numel(parameters.analysis.tfa.base.window == 2))
        baseWindow = parameters.analysis.tfa.base.window;
    end
    
    if (isempty(baseWindow) || baseWindow(1) >= baseWindow(2))
        disp('Invalid baseline range'); return;
    end
    
    if (~isempty_field(parameters,'parameters.analysis.tfa.base.trials'))
        baseTrials = parameters.analysis.tfa.base.trials;
        baseTrials = baseTrials(baseTrials >= 1 & baseTrials <= nTrials);
    end

    if (~isempty_field(parameters,'parameters.analysis.tfa.base.ncycles'))
        baseCycles = parameters.analysis.tfa.base.ncycles; 
    end
        
    cfg = [];
    cfg.cycles = baseCycles;
    cfg.trials = baseTrials;
    cfg.win    = baseWindow;
    
    baseline = getBaseline(timefreq,cfg);
end

switch parameters.analysis.tfa.base.type
    case 'absolute';   timefreq.powspctrm =  timefreq.powspctrm -  baseline;
    case 'relative';   timefreq.powspctrm =  timefreq.powspctrm ./ baseline;
    case 'relchange';  timefreq.powspctrm = (timefreq.powspctrm -  baseline) ./ baseline;
    case 'normchange'; timefreq.powspctrm = (timefreq.powspctrm -  baseline) ./ (timefreq.powspctrm + baseline);
    case 'decibel';    timefreq.powspctrm = 10 * log10(timefreq.powspctrm ./ baseline);
end

end

function baseline = getBaseline(data,cfg)

nDims = ndims(data.powspctrm);
if (nDims < 4); error('error:FaultyDimensions',['The "powspctrm" field has ' num2str(nDims) ' dimensions, but should have 4 dimensions.']); end

nFreqs  =  size(data.powspctrm,3);
nTimes  =  size(data.powspctrm,4);

T1 = cfg.win(1);
T2 = cfg.win(2);

% Extract window based on number of cycles, i.e. frequency dependent
minFreq = min(data.freq); % Lowest frequency to determine maximum window range
tc = cfg.cycles / minFreq; % Maximum window range
dt = T2 - tc;
T1 = max([dt T1]);

iT1  = find(data.time >= T1,1,'first'); if (isempty(iT1)); iT1 = 1;      end
iT2  = find(data.time <= T2,1,'last');  if (isempty(iT2)); iT2 = nTimes; end
idxT = iT1:iT2;
nT   = length(idxT);
del  = false(nFreqs,nT);

for iFreq = 1:nFreqs
    tc = cfg.cycles / data.freq(iFreq); % Baseline window based on frequency
    dt = T2 - tc;
    t1 = max([dt T1]);
    i1 = find(data.time(idxT) >= t1,1,'first'); % Start of sub-range  
    if (i1 < 0); i1 = 0; end % Can't be outside window range
    del(iFreq,1:i1-1) = true; % Tag everything before start of sub-range
end

% Extrace the baseline window
baseline = data.powspctrm(cfg.trials,:,:,idxT); % extract baseline window
baseline(:,:,del) = NaN; % Ignore data points that are outside of cycle window
baseline = nanmean(baseline,4); % average over time bins
baseline =  repmat(baseline,1,1,1,nTimes);

end