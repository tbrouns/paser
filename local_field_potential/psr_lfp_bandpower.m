function bandpower = psr_lfp_bandpower(timefreq,parameters)

% PSR_LFP_BANDPOWER - Calculate bandpower from power spectrum array
% 
% Syntax:  bandpower = psr_lfp_bandpower(timefreq,parameters)
%
% Inputs:
%    timefreq    - Output structure from PSR_LFP_TFA
%    parameters - See PSR_PARAMETERS_ANALYSIS
%
% Outputs:
%    bandpower - Bandpower array, with shape:
%                [Num. of trials x Num. of bands x Num. of time points]
%   
% See also: PSR_LFP_TFA

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

bp_ranges  = parameters.analysis.bpw.frange;
nRanges    = size(bp_ranges,1);
nTrials    = size(timefreq.powspctrm,1);
nTimes     = size(timefreq.powspctrm,4);
bandpower  = NaN(nTrials,nRanges,nTimes);
freqs      = timefreq.freq;

for iTrial = 1:nTrials
    for iRange = 1:nRanges
        f1 = find(freqs <= bp_ranges(iRange,1),1,'last');
        f2 = find(freqs >= bp_ranges(iRange,2),1,'first');
        powspctrm = timefreq.powspctrm(iTrial,:,f1:f2,:);
        powspctrm = nanmean(powspctrm,2); % Mean over channels
        powspctrm = nanmean(powspctrm,3); % Mean over frequencies in band
        powspctrm = squeeze(powspctrm);
        bandpower(iTrial,iRange,:) = powspctrm;
    end
end

end