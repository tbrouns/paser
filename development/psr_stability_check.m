function filesSpikes = psr_stability_check(filesSpikes,filesData,sortMethod)

nProbes    = size(filesData,1);
nBlocks    = size(filesData,2);
sortMethod = upper(sortMethod);

for iProbe = 1:nProbes
    
    % Check if variable exists
    MSGID = 'MATLAB:load:variableNotFound';
    warning('off', MSGID);
    load(filesData{iProbe,1},'perturb'); % Load first temp file
    warning('on', MSGID);
    if ~exist('perturb','var') || isempty(perturb); perturbTemp = 0; perturb = []; % No perturbation done yet
    else,                                           perturbTemp = perturb + 1;
    end
    
    % Save old spikes file with new name
    
    load(filesSpikes{iProbe,1});
    filename = getFilename(filesSpikes{iProbe,1},perturbTemp,sortMethod);
    save(filename,'spikes','metadata','parameters');
    
    blockOnsets = metadata.blockonset;
    
    for iBlock = 1:nBlocks
        
        if (perturbTemp < 2) % Perturb signal
            
            % Save old data
        
            load(filesData{iProbe,iBlock});
            filename = getFilename(filesData{iProbe,iBlock},perturbTemp,sortMethod);
            save(filename,'metadata','parameters','ts_LFP','ts_Spikes','perturb');
                    
            % Load unperturbed data
            
            filename = getFilename(filesData{iProbe,iBlock},0,sortMethod);
            load(filename);
            
            signal        = ts_Spikes.data;
            parameters.Fs = ts_Spikes.Fs;
            
            % Grab spikes from trial
            keep = find(spikes.blocks == iBlock);
            spikesBlock = psr_sst_remove_spikes(spikes,keep,'keep'); 
            spikesBlock.spiketimes = spikesBlock.spiketimes - blockOnsets(iBlock);
            
            if (perturbTemp == 0); signal = psr_stability_blurring(spikesBlock,signal,parameters);
            else,                  signal = psr_stability_reversal(spikesBlock,signal,parameters);
            end
            
            ts_Spikes.data = signal;
            
            % Save
            perturb = perturbTemp;
            save(filesData{iProbe,iBlock},'metadata','parameters','ts_LFP','ts_Spikes','perturb');
        else % Done stability check
            fileRAW = filesData{iProbe,iBlock};
            fileUNP = getFilename(fileRAW,0,sortMethod);
            fileBLR = getFilename(fileRAW,1,sortMethod);
            delete(fileRAW);
            delete(fileBLR);
            movefile(fileUNP,fileRAW);
            if (iBlock == 1)
                delete(filesSpikes{iProbe,1});
                filesSpikes{iProbe,1} = [];
            end
        end
    end
        
    filesSpikes{iProbe,4} = [];
    filesSpikes{iProbe,5} = [];
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