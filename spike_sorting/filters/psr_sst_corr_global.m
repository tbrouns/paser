function psr_sst_corr_global(filesSpikes,filesData)

% PSR_SST_CORR_GLOBAL - Calculates correlation of spikes across all channels. 
% This function calculates the correlation for each spike waveform across
% every channel of every probe. Artifact spikes are expected to be highly
% correlated, since the source of the signal is not proximate to be probe,
% unlike action potentials.
%
% Syntax:  psr_sst_corr_global(filesSpikes,filesData)
%
% Inputs:
%    filesSpikes - Files that contain "spikes" structures from which spikes
%    need to be extracted. 
%    filesData - Files that contain filtered time series of extracellular
%    recording, given by "spikes.data".
%
% Outputs:
%    Adds "correlations" field to "spikes" structure in files given by
%    "filesSpikes".
%
% See also: PSR_WRAPPER

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2017

%------------- BEGIN CODE --------------

% Initialize variables

spikes        = [];
parameters    = [];
spikeTimesAll = [];
probeIDAll    = [];

filesSpikes = filesSpikes(:,1);
k = find(~cellfun(@isempty,filesSpikes));
filesSpikes = filesSpikes(k);
filesData   = filesData(k,:);

% Set parameters

iFile = 1;
while (isempty(spikes))
    load(filesSpikes{iFile}); % Load file to extract parameters
    iFile = iFile + 1;
end

if (isempty(spikes)); disp('No spikes found. Exiting function.'); return; end

nProbes   = length(filesSpikes);
nTrials   = size(filesData,2);
nChans    = size(spikes.waveforms,3);
nChansTot = nProbes * nChans;
Fs        = spikes.Fs;

tLength  = sum(spikes.info.dur);
sLength  = floor(Fs * tLength);
tSection = 60 * parameters.spikes.twin; % Cut data in sections
nSection = ceil(tLength / tSection); % Process data in sections
tSection = tLength / nSection; % Process data in sections
sWindow  = floor(0.5 * Fs * (parameters.spikes.window_size / 1000)); % in samples

% Extract all spikes

MSGID = 'MATLAB:load:variableNotFound';
warning('off', MSGID);
for iProbe = 1:nProbes
    spikes = [];
    load(filesSpikes{iProbe},'spikes','parameters');
    if (~isfield(spikes,'corr_global') && isfield(spikes,'spiketimes'))
        spikeTimesAll = [spikeTimesAll, spikes.spiketimes]; %#ok
        probeIDAll = [probeIDAll, iProbe * ones(size(spikes.spiketimes),'int16')]; %#ok
    end
end
warning('on', MSGID);

if (isempty(spikeTimesAll)); return; end

% Sort spikes in chronological order

[spikeTimesAll,I] = sort(spikeTimesAll);
probeIDAll = probeIDAll(I); clear I;

% Process data in sections

nSpikesAll = length(spikeTimesAll);
correlationsAll = zeros(1,nSpikesAll,'single');
tStart = 0;
iCorr  = 1; % array iterator
for iSection = 1:nSection
    
    disp(['Processed ' num2str(round(100 * (iCorr - 1) / nSpikesAll),'%02d') '% of spikes']);
    
    % Set section limits
    
    tEnd   = tStart + tSection; 
    iStart = floor(tStart * Fs) - sWindow;
    iEnd   =  ceil(tEnd   * Fs) + sWindow;
    if (iStart < 1);       iStart = 1;       end
    if (iEnd   > sLength); iEnd   = sLength; end
    
    % Extract spikes within section
    
    id = (spikeTimesAll >= tStart & spikeTimesAll < tEnd);
    spikeTimes = spikeTimesAll(id) - tStart;
    
    % Load raw data
    
    sDataSection = iEnd - iStart + 1;
    data = zeros(nChansTot,sDataSection,'int16');
    
    for iProbe = 1:nProbes
        iOffset = 0;
        itr     = 0; % iterator indicating end-position of last trial
        dataProbe = []; % Minimum data that we need to load
        for iTrial = 1:nTrials
            if (itr < iEnd) % Check if we still need to load
                load(filesData{iProbe,iTrial},'metadata'); % Only interested in trial length
                sLengthTrial = Fs * metadata.duration + 1;
                itr = itr + sLengthTrial;
                if (itr + sLengthTrial > iStart) % Check if we need to start loading
                    load(filesData{iProbe,iTrial},'ts_Spikes');
                    dataProbe = [dataProbe, ts_Spikes.data];
                else % Move offset to first trial within range
                    iOffset = itr;
                end
            else % Outside upper limit
                break;
            end
        end        
        iChanStart = nChans * (iProbe - 1) + 1;
        iChanEnd   = nChans * iProbe;
        data(iChanStart:iChanEnd,:) = dataProbe(:,iStart-iOffset:iEnd-iOffset);
    end
          
    % Align adjacent spikes
    
    [spikeTimes,spikeIDs] = psr_sst_align(parameters,spikeTimes);
        
    % Convert to sample number
    
    spikeTimes = round(Fs * spikeTimes) + 1;
    spikeTimes(spikeTimes < 1 + sWindow)            = 1 + sWindow;
    spikeTimes(spikeTimes > sDataSection - sWindow) = sDataSection - sWindow;
    
    % Calculate correlations
    
    nspikes = length(spikeTimes);    
    samples = bsxfun(@plus,spikeTimes,(-sWindow:sWindow));
    data = data(:,samples');
    data = mat2cell(data,nChansTot,size(samples,2)*ones(1,nspikes));
    correlations = cellfun(@corrMean,data);
    correlations = correlations(spikeIDs);
    
    % Save
    
    nspikes = length(spikeIDs);
    correlationsAll(iCorr:iCorr+nspikes-1) = single(correlations);
    iCorr = iCorr + nspikes;
    
    tStart = tEnd; % Next data section
    
end

% Save correlations for each probe
MSGID = 'MATLAB:load:variableNotFound';
warning('off', MSGID);
for iProbe = 1:nProbes
    spikes = [];
    load(filesSpikes{iProbe});
    if (~isempty(spikes))
        spikes.corr_global = correlationsAll(probeIDAll == iProbe);
        save(filesSpikes{iProbe},'spikes','-append');
    end
end
warning('on', MSGID);

end

function [spikeTimesAligned,spikeIDs] = psr_sst_align(parameters,spikeTimes)

% Combine and align adjacent spikes

spikeTimes        =   sort(spikeTimes);
nspikes           = length(spikeTimes);
spikeDur          = (parameters.spikes.max_desync * parameters.spikes.window_size) / 1000; % sec
spikeTimesAligned = -1 * ones(nspikes,1);
spikeIDs          = zeros(nspikes,1);
iSpike            = 1;
kSpike            = 1;
T                 = spikeTimes(1);

while (iSpike <= nspikes)
    spikeTimesTemp = [];
    jSpike = 1;
    while (iSpike <= nspikes) % end condition
        t = spikeTimes(iSpike);
        dt = t - T; T = t;
        if (dt > spikeDur); break; end
        spikeTimesTemp(jSpike) = t; %#ok
        spikeIDs(iSpike) = kSpike;
        iSpike = iSpike + 1;
        jSpike = jSpike + 1;
    end
    spikeTimesAligned(kSpike) = mean(spikeTimesTemp);
    kSpike = kSpike + 1;
end

spikeTimesAligned(spikeTimesAligned < 0) = [];

end

function correlation = corrMean(x)

x = single(x');
R = triu(corr(x),1);
R = R(triu(true(size(R)),1));
correlation = mean(R);
        
end

%------------- END OF CODE --------------