function psr_sst_spikes_info(filesSpikes,filesData)

filesSpikes = filesSpikes(:,1);
k = find(~cellfun(@isempty,filesSpikes));
filesSpikes = filesSpikes(k);
filesData   = filesData(k,:);

nProbes = length(filesSpikes);
nTrials = size(filesData,2);

for iProbe = 1:nProbes
    
    load(filesSpikes{iProbe},'spikes');   
    
    if (~min(isfield(spikes.info,{'stds','thresh','dur'}))) 

        for iTrial = 1:nTrials

            load(filesData{iProbe,iTrial},'ts_Spikes','parameters');
            ts_Spikes.data = psr_single(ts_Spikes.data,parameters);
            
            if (iTrial == 1) % Initialize
                nChans = size(ts_Spikes.data,1); 
                stdAll = zeros(nTrials,nChans);
                thrAll = zeros(nTrials,nChans);
                durAll = zeros(nTrials,1);
            end

            stdAll(iTrial,:) = std(ts_Spikes.data,[],2)';
            thrAll(iTrial,:) = -parameters.spikes.thresh * psr_mad(ts_Spikes.data)';
            durAll(iTrial)   = size(ts_Spikes.data,2) / ts_Spikes.Fs;       
        end

        spikes.info.stds   = stdAll;
        spikes.info.thresh = thrAll;
        spikes.info.dur    = durAll;
        save(filesSpikes{iProbe},'spikes','-append');
    end
end