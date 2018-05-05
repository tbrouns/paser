function data = psr_lfp_baseline(data,baseline,parameters)

% Based on FT_FREQBASELINE

if (ndims(baseline) ~= ndims(data.powspctrm) || all(size(baseline) ~= size(data.powspctrm)))
    base_twin = [];
    if (strcmp(parameters.analysis.tfa.baseline,'pre'))
        base_twin = [data.time(1) 0];
    elseif (numel(parameters.analysis.tfa.baseline == 2))
        base_twin = parameters.analysis.tfa.baseline; 
    end
    
    if (isempty(base_twin) || base_twin(1) >= base_twin(2))
        disp('Invalid baseline range'); return; 
    end

    baseline = getBaseline(data,base_twin);
end
    
switch parameters.analysis.tfa.baselineType
    case 'absolute';   data.powspctrm =  data.powspctrm -  baseline;
    case 'relative';   data.powspctrm =  data.powspctrm ./ baseline;
    case 'relchange';  data.powspctrm = (data.powspctrm -  baseline) ./ baseline;
    case 'normchange'; data.powspctrm = (data.powspctrm -  baseline) ./ (data.powspctrm + baseline);
    case 'decibel';    data.powspctrm = 10 * log10(data.powspctrm ./ baseline);
end

end

function baseline = getBaseline(data,base_twin)

t1 = base_twin(1);
t2 = base_twin(2);

baseline = [];
nTrials  =   size(data.powspctrm,1);
nDims    =  ndims(data.powspctrm);
nTimes   =   size(data.powspctrm,nDims);
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