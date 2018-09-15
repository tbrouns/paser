function spikes = psr_freq2spikes(spikes,freq)

% PSR_FREQ2SPIKES - Adds LFP data to the spike data structure
%
% Syntax:  spikes = psr_freq2spikes(spikes,freq)
%
% Inputs:
%    spikes - See README
%    freq   - See README
%
% Outputs:
%    spikes - We add an LFP artifact field to the "spikes" structure
%
% See also: PSR_WRAPPER

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

if (isfield(freq,'artifacts'))
    nTrials = size(freq,2);
    for iTrial = 1:nTrials
        x = freq(iTrial).artifacts;
        fields = fieldnames(x);
        nFields = size(fields,1);
        artifacts = [];
        for iField = 1:nFields
            artifacts.(fields{iField}) = (freq(iTrial).artifacts.(fields{iField}) - 1) / freq(iTrial).fsample; 
        end
        spikes.info.artifacts.lfp{iTrial,1} = artifacts;
    end
end
end