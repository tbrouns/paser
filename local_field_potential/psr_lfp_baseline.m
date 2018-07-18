function data = psr_lfp_baseline(data,baseline,parameters)

% Based on FT_FREQBASELINE

if (ndims(baseline) ~= ndims(data.powspctrm) || all(size(baseline) ~= size(data.powspctrm)))
    
    nTrials = size(data.powspctrm,1);
    baseTrials = 1:nTrials;
    baseCycles = [];
    baseWindow = [];
       
    if (strcmp(parameters.analysis.tfa.base.window,'pre'))
        baseWindow = [data.time(1) 0];
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

    if (~isempty_field(parameters,'parameters.analysis.tfa.base.ncycles')); baseCycles = parameters.analysis.tfa.base.ncycles; end
        
    cfg = [];
    cfg.cycles = baseCycles;
    cfg.trials = baseTrials;
    cfg.win    = baseWindow;
    
    baseline = getBaseline(data,cfg);
end

switch parameters.analysis.tfa.base.type
    case 'absolute';   data.powspctrm =  data.powspctrm -  baseline;
    case 'relative';   data.powspctrm =  data.powspctrm ./ baseline;
    case 'relchange';  data.powspctrm = (data.powspctrm -  baseline) ./ baseline;
    case 'normchange'; data.powspctrm = (data.powspctrm -  baseline) ./ (data.powspctrm + baseline);
    case 'decibel';    data.powspctrm = 10 * log10(data.powspctrm ./ baseline);
end

end

function baseline = getBaseline(data,cfg)

nDims = ndims(data.powspctrm);
if (nDims < 4); error('error:FaultyDimensions',['The "powspctrm" field has ' num2str(nDims) ' dimensions, but should have 4 dimensions.']); end

nTrials =  size(data.powspctrm,1);
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
baseline = nanmean(baseline,1); % average over trials
baseline = nanmean(baseline,4); % average over time bins
baseline =  repmat(baseline,nTrials,1,1,nTimes);

end