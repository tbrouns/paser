function output = psr_ft_convert2fieldtrip(input,parameters)

output = []; % Field trip data format

if (isfield(input,'data') && isfield(input,'timestamps'))
    
    % Convert LFP data to FT_DATATYPE_RAW format
    
    % See FieldTrip documentation:
    % http://www.fieldtriptoolbox.org/reference/ft_datatype_raw
    % http://www.fieldtriptoolbox.org/faq/how_can_i_import_my_own_dataformat
    
    data       = double(input.data);
    timestamps = double(input.timestamps);
    
    % Function that converts dataformat to field trip format
    
    Fs      = parameters.Fs;
    nchans  = size(data,1);
    labels  = 1:nchans;
    labels  = strtrim(cellstr(num2str(labels'))');
    
    output.label   = labels;              % Cell-array containing strings, Nchan*1
    output.fsample = Fs;                  % Sampling frequency in Hz, single number
    output.trial   = {data};              % Cell-array containing a data matrix for each trial (1 X Ntrial); each data matrix is a [Nchan X Npoints] matrix
    output.time    = {timestamps};        % Cell-array containing a time axis   for each trial (1 X Ntrial); each time axis   is a [    1 X Npoints] vector
    output.sampleinfo = [1 size(data,2)]; % Array containing [startsample endsample] of data
    
elseif (isfield(input,'spikes'))
    
    %     Convert spike data to FT_DATATYPE_SPIKE format
    %
    %     See FieldTrip documentation:
    %     http://www.fieldtriptoolbox.org/reference/ft_datatype_spike
    
    % input.onset: stimulus onset (t = 0)
    
    spikes   = input.spikes;
    nTrials  = size(input.trialtime,1);
    clustIDs = unique(spikes.assigns);
    nClusts  = length(clustIDs);
    
    if (nClusts == 0); return; end
    
    if (isfield(input,'probeID')); probeID = input.probeID;
    else,                          probeID = [];
    end
    
    itr = 1;
    
    for iClust = 1:nClusts
        
        clustID  = clustIDs(iClust);
        spikeIDs = ismember(spikes.assigns, clustID);
        
        waveformsTemp = spikes.waveforms (spikeIDs,:,:);
        timestampTemp = spikes.spiketimes(spikeIDs);
        
        trials = spikes.trials(:,spikeIDs);
        
        % Deal with spikes contained within more than one trial
        
        timestamp = [];
        trialidxs = [];
        waveforms = [];
        for iTrial = 1:nTrials
            spikeIDs_trial = find(trials(iTrial,:));
            trialidxs = [trialidxs,iTrial * ones(size(spikeIDs_trial))];
            timestamp = [timestamp,     timestampTemp(spikeIDs_trial)];
            waveforms = [waveforms;     waveformsTemp(spikeIDs_trial,:,:)];
        end
        
        time = timestamp;
        for iTrial = unique(trialidxs)
            spikeIDs_trial = ismember(trialidxs,iTrial);
            time(spikeIDs_trial) = timestamp(spikeIDs_trial) - input.onsets(iTrial);
        end
        
        output.id       {itr} = clustID;
        output.label    {itr} = ['P' num2str(probeID) '-Unit-' num2str(clustID)];
        output.waveforms{itr} = permute(waveforms,[3 2 1]);
        output.timestamp{itr} = timestamp;
        output.trial    {itr} = trialidxs;
        output.time     {itr} = time;
                
        itr = itr + 1;
    end
    
    output.trialtime = input.trialtime;
end