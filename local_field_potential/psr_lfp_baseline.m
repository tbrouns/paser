function data = psr_lfp_baseline(data,cfg)

if (all(size(cfg.baseline) == size(data.powspctrm)))
    baseline = cfg.baseline;
else
    baseline = getBaseline(data.powspctrm,cfg.base_t0,cfg.base_t1);
end
    
switch cfg.powtype
    case 'absolute';   data.powspctrm =  data.powspctrm -  baseline;
    case 'relative';   data.powspctrm =  data.powspctrm ./ baseline;
    case 'relchange';  data.powspctrm = (data.powspctrm -  baseline) ./ baseline;
    case 'normchange'; data.powspctrm = (data.powspctrm -  baseline) ./ (data.powspctrm + baseline);
    case 'decibel';    data.powspctrm = 10 * log10(data.powspctrm ./ baseline);
end

end

function baseline = getBaseline(data,t1,t2)

if (nargin < 2); t1 = -Inf; end
if (nargin < 3); t2 =  Inf; end

baseline = [];
nTrials  =   size(data.powspctrm,1);
nDims    =  ndims(data.powspctrm);
nTimes   = length(data.time);
i1 = find(data.time <= t1,1,'last');
i2 = find(data.time >= t2,1,'first');
if (isempty(i1)); i1 = 1; end
if (isempty(i2)); i2 = nTimes; end

if (nDims == 4)
    baseline = data.powspctrm(:,:,:,i1:i2); % extract baseline window
    baseline = nanmean(baseline,1); % average over trials
    baseline = nanmean(baseline,4); % average over time bins
    baseline =  repmat(baseline,nTrials,1,1,nTimes);
elseif (nDims == 3)
    baseline = data.powspctrm(:,:,i1:i2);
    baseline = nanmean(baseline,1); % average over trials
    baseline = nanmean(baseline,3); % average over time bins
    baseline =  repmat(baseline,nTrials,1,nTimes);
end

end