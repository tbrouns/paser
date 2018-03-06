function bandpwrNew = psr_lfp_combine_bp(bandpwr,parameters)

bandpwrNew = [];
nTrials = size(bandpwr,1);
for iTrial = 1:nTrials
    missing = bandpwr{iTrial,:,2};
    bp      = bandpwr{iTrial,:,1};
    keep    = missing <= parameters.analysis.bpw.maxmiss;
    bandpwrNew = [bandpwrNew;bp(keep,:)];
end

end