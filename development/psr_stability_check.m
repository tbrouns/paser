function filesSpikes = psr_stability_check(filesSpikes,filesData,sortMethod)

nProbes    = size(filesData,1);
nTrials    = size(filesData,2);
sortMethod = upper(sortMethod);

for iProbe = 1:nProbes
    
    % Check if variable exists
    MSGID = 'MATLAB:load:variableNotFound';
    warning('off', MSGID);
    load(filesData{iProbe,1},'perturb'); % Load first temp file
    warning('on', MSGID);
    if ~exist('perturb','var') || ...
            isempty(perturb); perturbTemp = 0; perturb = []; % No perturbation done yet
    else,                     perturbTemp = perturb + 1;
    end
    
    % Save old spikes file with new name
    
    load(filesSpikes{iProbe,1});
    filename = getFilename(filesSpikes{iProbe,1},perturbTemp,sortMethod);
    save(filename,'spikes','metadata','parameters');
    
    trialOnsets = metadata.trialonset;
    
    for iTrial = 1:nTrials
        
        % Perturb signal
        
        if (perturbTemp < 2)
            
            % Save old data
        
            load(filesData{iProbe,iTrial});
            filename = getFilename(filesData{iProbe,iTrial},perturbTemp,sortMethod);
            save(filename,'metadata','parameters','ts_LFP','ts_Spikes','perturb');
                    
            % Load unperturbed data
            
            filename = getFilename(filesData{iProbe,iTrial},0,sortMethod);
            load(filename);
            
            signal        = ts_Spikes.data;
            parameters.Fs = ts_Spikes.Fs;
            which = find(spikes.trials == iTrial);
            spikesTrl = psr_sst_spike_removal(spikes,which,'keep');
            spikesTrl.spiketimes = spikesTrl.spiketimes - trialOnsets(iTrial);
            
            if (perturbTemp == 0); signal = psr_stability_blurring(spikesTrl,signal,parameters);
            else,                  signal = psr_stability_reversal(spikesTrl,signal,parameters);
            end
            
            ts_Spikes.data = signal;
            
            % Save
            perturb = perturbTemp;
            save(filesData{iProbe,iTrial},'metadata','parameters','ts_LFP','ts_Spikes','perturb');
        else
            fileRAW = filesData{iProbe,iTrial};
            fileUNP = getFilename(fileRAW,0,sortMethod);
            fileBLR = getFilename(fileRAW,1,sortMethod);
            delete(fileRAW);
            delete(fileBLR);
            movefile(fileUNP,fileRAW);
            if (iTrial == 1)
                delete(filesSpikes{iProbe,1});
                filesSpikes{iProbe,1} = [];
            end
        end
    end
        
    filesSpikes{iProbe,3} = [];
    clear perturb;
end

end

function filename = getFilename(file,perturb,sortMethod)

[fpath,filename,~] = fileparts(file);

if     (perturb == 0); filename = [filename '_' sortMethod '_UNP'];
elseif (perturb == 1); filename = [filename '_' sortMethod '_BLR'];
else,                  filename = [filename '_' sortMethod '_RVS'];
end

filename = [fpath '\' filename '.mat'];

end