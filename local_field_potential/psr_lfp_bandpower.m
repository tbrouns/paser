function bandpower = psr_lfp_bandpower(timefrq,parameters)

bp_ranges  = parameters.analysis.bpw.frange;
nRanges    = size(bp_ranges,1);
nReps      = size(timefrq.powspctrm,1);
nTimes     = size(timefrq.powspctrm,4);
bandpower  = NaN(nReps,nRanges,nTimes);
freqs      = timefrq.freq;

for iRep = 1:nReps
    for iRange = 1:nRanges
        f1 = find(freqs <= bp_ranges(iRange,1),1,'last');
        f2 = find(freqs >= bp_ranges(iRange,2),1,'first');
        powspctrm = timefrq.powspctrm(iRep,:,f1:f2,:);
        powspctrm = nanmean(powspctrm,2); % Mean over channels
        powspctrm = nanmean(powspctrm,3); % Mean over frequencies in band
        powspctrm = squeeze(powspctrm);
        bandpower(iRep,iRange,:) = powspctrm;
    end
end

end