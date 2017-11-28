function [filesSpikes,filesData] = psr_stability_check(filesSpikes,filesData)

nProbes = size(filesData,1);
nTrials = size(filesData,2);

for iProbe = 1:nProbes
    
    % Check if variable exists
    MSGID = 'MATLAB:load:variableNotFound';
    warning('off', MSGID);
    load(filesData{iProbe,1},'perturb'); % Load first temp file
    warning('on', MSGID);
    if ~exist('perturb','var'); perturb = 0; % No perturbation done yet
    else,                       perturb = perturb + 1;
    end
    
    load(filesSpikes{iProbe}); % Load spikes
    
    % Save new spikes file
    
    vars = {'spikes','metadata','parameters'};
    filesSpikes{iProbe} = saveFile(filesSpikes{iProbe},perturb,vars,spikes,[],metadata,parameters);
    
    % Perturb signal
    
    if (perturb < 2)
        
        for iTrial = 1:nTrials
            load(filesData{iProbe,iTrial},'ts_Spikes');
            
            psr_parameter_config; % TEMP
            
            signal = ts_Spikes.data;
            parameters.Fs = ts_Spikes.Fs;
            which = find(spikes.trials == iTrial);
            spikesTrl = psr_sst_spike_removal(spikes,which,'keep');
            spikesTrl.spiketimes = spikesTrl.spiketimes - metadata.trialonset(iTrial);
            
            if (perturb == 0); signal = psr_stability_blurring(spikesTrl,signal,parameters);
            else,              signal = psr_stability_reversal(spikesTrl,signal,parameters);
            end
        
            ts_Spikes.data = signal;
            vars = {'ts_Spikes','perturb'};
            filesData{iProbe,iTrial} = saveFile(filesData{iProbe,iTrial},perturb,vars,[],ts_Spikes,metadata,parameters);
        end
    end
    
    filesSpikes{iProbe,3} = [];
    clear perturb;
end

end

function filename = saveFile(file,perturb,vars,spikes,ts_Spikes,metadata,parameters)

[fpath,filename,~] = fileparts(file);

if     (perturb == 0); filename = [filename,'_UNP'];
elseif (perturb == 1); filename = [filename,'_BLR'];
else,                  filename = [filename,'_RVS'];
end

filename = [fpath '\' filename '.mat'];
save(filename,vars{:});
    
end