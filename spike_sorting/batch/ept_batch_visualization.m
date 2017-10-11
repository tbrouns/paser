function ept_batch_visualization(subject,loadPath_root,savePath_root,expType,plot_analysis,plot_spiking)

%%% TEMP

GPSTH_all = [];
PSTH_all  = [];

%%%

if (nargin < 2); loadPath_root = [];    end
if (nargin < 3); savePath_root = [];    end
if (nargin < 4); expType = 'all';       end
if (nargin < 5); plot_analysis = false; end
if (nargin < 6); plot_spiking  = false; end

ptrn = 'passive';

folder_names = dir([loadPath_root 'R*']);
folder_names = char(folder_names.name);

%% sort files

numfolders = length(folder_names(:,1));

folders_passive = cell(1,1);
folders_active  = cell(1,1);
for iFolder = 1:numfolders
    foldername = lower(folder_names(iFolder,:));
    foldername = strtrim(foldername);
    k          = strfind(foldername,ptrn);
    if isempty(k); folders_active {iFolder,1} = folder_names(iFolder,:);
    else,          folders_passive{iFolder,1} = folder_names(iFolder,:);
    end
end

folders_active  = folders_active (~cellfun('isempty',folders_active )); % remove empty cells
folders_passive = folders_passive(~cellfun('isempty',folders_passive)); % remove empty cells

if     (strcmp(expType,'all'));     folders = [folders_active;folders_passive];
elseif (strcmp(expType,'active'));  folders = folders_active;
elseif (strcmp(expType,'passive')); folders = folders_passive;
end

numfolders = length(folders);

for iFolder = 14:numfolders
    foldername = strtrim(folders{iFolder});
    
    loadPath = [loadPath_root, foldername];
    savePath = [savePath_root, foldername];
    
    if (loadPath(end) ~= '\'); loadPath = [loadPath, '\']; end %#ok
    if (savePath(end) ~= '\'); savePath = [savePath, '\']; end %#ok
    
%     %%%% TEMP
%     
%     MATfiles  = dir([loadPath '\data.mat']);
%     filename  = char(MATfiles.name);
%     
%     if (isempty(filename)); continue; end
%     
%     load([loadPath filename]);
%     
%     GPSTH_all = [GPSTH_all,GPSTH];
%     PSTH_all  = [ PSTH_all, PSTH];
%     
%     continue;
%     
%     %%%%
        
    baseline  = [];
    MATfiles  = dir([loadPath '\Spikes_' subject '*.mat']);
    numfiles  = size(MATfiles,1);
    filenames = char(MATfiles.name);
    if (isempty(filenames)); continue; end
    
    if (plot_analysis)
        [GPSTH,PSTH] = ept_batch_analysis(filenames,loadPath,savePath,plot_analysis); %#ok
        save([loadPath 'data.mat'],'GPSTH','PSTH');
    end
    
    if (plot_spiking)
        for iFile = 1:numfiles
            close all
            filename = filenames(iFile,:);
            load([loadPath filename]);
            [~,filename,~] = fileparts(filename);
            parameters = ept_parameter_config(); % TEMP
            if (isfield(spikes,'spiketimes'))
                
                ept_sst_plotting(spikes,parameters,freq,[],savePath,filename); 
                suptitle(['Tetrode: ' num2str(metadata.tetrode)]);
    
                % PRE-MERGE 
                spikes.assigns = spikes.assigns_prior;
                spikes = ept_sst_clusterfeatures(spikes,freq,parameters);
                ept_sst_plotting(spikes,parameters,freq,[],savePath,filename,'_preMerge'); 
                suptitle(['Tetrode: ' num2str(metadata.tetrode)]);
                
                %             ept_plotting_cls(spikes,clusters,      savePath,filename);
            end
            
            %%% FIX LFP PLOTTING
%             if (~isempty(spikes.stimtimes{1}))
%                 if (metadata.stimulus == 0); ept_plotting_lfp(freq,parameters,savePath,filename); baseline = getBaseline(freq); % must be loaded first
%                 else,                        ept_plotting_lfp(freq,parameters,savePath,filename,baseline);
%                 end
%             end

        end
    end
end

end

function base = getBaseline(freq) % TEMP
    base = freq.powspctrm;
    base = nanmean(base,3);
    base = repmat(base,1,1,size(freq.time,2));
end