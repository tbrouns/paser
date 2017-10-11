function ept_sst_artifact_subtraction(files,data)

% Subtract global median of other tetrodes from tetrode channels

nprobes     = length(files);
nchans      = size(data,1);
nlength     = size(data,2);
nelectrodes = nchans / nprobes;

if (floor(nelectrodes) ~= nelectrodes) % Integer check
    disp('Number of electrodes not a multiple of number of probes.'); 
    return; 
end  

dataFiltered = single(zeros(size(data)));

for iProbe = 1:nprobes
    chanIDs = false(nchans,1);
    iFile = (iProbe - 1) * nelectrodes + 1;
    jFile = iProbe * nelectrodes;
    chanIDs(iFile:jFile) = 1;
    globalMedian = median(data(~chanIDs,:));
    dataTetrode  = data(chanIDs,:);
    dataTetrode  = bsxfun(@minus,dataTetrode,globalMedian);
    dataFiltered(chanIDs,:) = dataTetrode;    
end

clear data;

% Set new waveforms

for iProbe = 1:nprobes
    
    load(files{iProbe});
    
    if (isfield(spikes,'spiketimes'))
        Fs = spikes.params.Fs;
        window_samples = round(Fs * parameters.spikes.window_size / 1000);
        samples_hwidth = round(0.5 * window_samples);
        spiketimes = round(Fs * spikes.spiketimes);
        spiketimes = bsxfun(@plus,spiketimes,-samples_hwidth:samples_hwidth);
        spiketimes(spiketimes < 1)       = 1;
        spiketimes(spiketimes > nlength) = nlength; 
        waveforms = data(:, spiketimes);
        spikes.waveforms = waveforms;
    elseif (isfield(spikes,'data'))
        iFile = (iProbe - 1) * nelectrodes + 1;
        jFile = iProbe * nelectrodes;
        spikes.data = dataFiltered(iFile:jFile,:);
    end
    
    save(files{iProbe},'spikes','metadata','freq','parameters');
    
end

end