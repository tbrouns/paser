function output = psr_example_analysis(cfg)

output = []; % Dummy output. You can set your own output anywhere in the function

% Parse input
if (~isfield(cfg,'configPath')); cfg.configPath = []; end

% Redefine input
loadPath  = cfg.loadpath;
savePath  = cfg.savepath;

% Find the processed data files
filenames = dir([loadPath 'PSR_*.mat']);
if (size(filenames,1) == 0); return; end
filenames = char(filenames.name);

% Initialize
data_SPK   = cell(0,0);
data_LFP   = cell(0,0);

nProbes = size(filenames,1);
for iProbe = nProbes:-1:1 % Convert data for each probe
    
    % Load the probe data
    
    filename = filenames(iProbe,:);
    filepath = [loadPath filename];
    load(filepath);
    
    % Add FieldTrip to path
    [~,~] = psr_ft_path(parameters,'add');
    
    % Load parameters. We have to run this function AFTER loading the
    % data, because we have to add to the loaded parameters
    
    parameters.analysis.configPath = cfg.configPath; % Path to your custom parameters script
    parameters = psr_load_parameters(parameters,'analysis'); % Loads your custom parameters script. If no path is provided, we load "psr_parameters_analysis" instead
    
    nBlocks = size(metadata.stimtimes,1); % The number of experimental blocks
    
    for iBlock = nBlocks:-1:1 % Process each block
        
        % Grab the data for the current block
        spikeIDs = find(spikes.blocks == iBlock);
        block_spikes = psr_sst_remove_spikes(spikes,spikeIDs,'keep'); % Only keep the spikes in the current block
        block_freq   = freq(iBlock);
        
        % Stimulus timestamps for the block
        stimTimes = metadata.stimtimes(iBlock,:);
        
        % Offset [t = 0] to start of block
        block_spikes.spiketimes = block_spikes.spiketimes - metadata.blockonset(iBlock);
        
        % Only accept high quality isolated single units
        block_spikes = qualityThresholding(block_spikes,1,3);
        
        % Extract a window around each stimulus
        [SPK,LFP] = psr_stimulus_window(block_spikes,block_freq,stimTimes,parameters);
        
        % Convert spike data to FieldTrip format
        input         = [];
        input.spikes  = SPK;
        input.probeID = iProbe;
        input.popUnit = true; % Save a population reponse as well, which is going to be the first unit
        SPK           = psr_ft_convert2fieldtrip(input);
        
        % Save output
        data_SPK{iProbe,iBlock} = SPK;
        data_LFP{iProbe,iBlock} = LFP;
        
    end
end

% Do the analysis and plot some figures

close all;
fig_1 = figure; set(gcf,'position',get(0,'screensize'));
fig_2 = figure; set(gcf,'position',get(0,'screensize'));

id_subject = metadata.subject;    % Subject name 
id_session = metadata.session{1}; % Session name 
id_blocks  = metadata.stimulus;   % Stimulus ID

% Create a new folder to store the generated figures
savePathFigs = strsplit(savePath,'\');
savePathFigs = cell2mat(join(savePathFigs(1:end-2),'\'));
savePathFigs = [savePathFigs,'\analysis_figures\'];
[~,~,~] = mkdir(savePathFigs);

for iProbe = 1:nProbes
    for iBlock = 1:nBlocks
        
        SPK = data_SPK{iProbe,iBlock};
        LFP = data_LFP{iProbe,iBlock};
        
        %% Spikes

        % Analyze
        psth  = psr_ft_psth (SPK,parameters); % Peri-stimulus time histogram
        isih  = psr_ft_isi  (SPK,parameters); % Inter-spike interval histogram
        XCorr = psr_ft_xcorr(SPK,parameters); % Correlogram
        [jpsth,jpsthShuff] = psr_ft_jpsth(psth,parameters); % Joint peri-stimulus time histogram

        % Visualize

        unitIDs = SPK.label; % Label of each unit
        nUnits  = length(unitIDs); % Number of units for this probe
        startID = 2; % We ignore the population response for this example

        for iUnit = startID:nUnits

            %% ISIH and PSTH

            figure(fig_1); clf;
            subplot(1,2,1); psr_ft_plot_isih  (    isih,parameters,iUnit); % Plot inter-spike interval histogram
            subplot(1,2,2); psr_ft_plot_raster(SPK,psth,parameters,iUnit); % Plot peri-stimulus time histogram + raster plot

            % Super label
            suplabel([                             ...
                'ISIH and PSTH for: '              ...
                'Subject: ' id_subject        ', ' ...
                'Session: ' id_session        ', ' ...
                'Block : '  id_blocks{iBlock} ', ' ...
                'Probe #: ' num2str(iProbe)   ', ' ...
                'Unit: '    unitIDs{iUnit}],[],[],'none');

            % Save the figure
            savestr = [savePathFigs ...
                'ISIH-PSTH'                '_' ...
                id_subject                 '_' ...
                id_session                 '_' ...
                id_blocks{iBlock}          '_' ...
                'P' num2str(iProbe,'%02d') '_' ...
                unitIDs{iUnit}];

            export_fig(savestr);
            savefig([savestr '.fig']);

            %% Correlogram and JPSTH

            for jUnit = iUnit:nUnits

                figure(fig_2); clf;
                subplot(1,2,1); psr_ft_plot_xcorr(XCorr,parameters,[iUnit jUnit]); % Plot correlogram
                subplot(1,2,2); psr_ft_plot_jpsth(jpsth,parameters,[iUnit jUnit]); % Plot joint peri-stimulus time histogram

                % Super label
                suplabel([                                ...
                    'Correlogram and JPSTH for: '         ...
                    'Subject: ' id_subject        ', '    ...
                    'Session: ' id_session        ', '    ...
                    'Block: '   id_blocks{iBlock} ', '    ...
                    'Probe #: ' num2str(iProbe)   ', '    ...
                    'Unit: '    unitIDs{iUnit}    ', vs ' ...
                    'Unit: '    unitIDs{jUnit}],[],[],'none');

                % Save the figure
                savestr = [savePathFigs ...
                    'CORR-JPSTH'               '_' ...
                    id_subject                 '_' ...
                    id_session                 '_' ...
                    id_blocks{iBlock}          '_' ...
                    'P' num2str(iProbe,'%02d') '_' ...
                    unitIDs{iUnit}             '-' ...
                    unitIDs{jUnit}];

                export_fig(savestr);
                savefig([savestr '.fig']);
            end
        end
        
        %% Local field potential
        
        % Set times of interest for sliding window
        % Normally you would set these parameters in your custom script for analysis parameters
        t = LFP.time{1};
        parameters.analysis.tfa.toi = t(1):0.01:t(end);
        parameters.analysis.tfa.pad = 10;
        parameters.analysis.tfa.base.type    = 'decibel';
        parameters.analysis.tfa.base.ncycles = 5;
        parameters.analysis.tfa.base.window  = [-1.00 0];
        
        % Time-frequency analysis
        timefreq = psr_lfp_tfa(LFP,parameters); % Time-frequency analysis
        timefreq = psr_lfp_baseline(timefreq,[],parameters); % Subtract pre-stimulus baseline
        
        % Take median over trials, ignoring NaNs. We have to take the median or the mean
        % before we plot the spectrogram to avoid errors.
        timefreq.powspctrm = nanmedian(timefreq.powspctrm,1);
        
        % Visualize
        figure(fig_1); clf;
        parameters.analysis.tfa.plot.tlim = [-0.25 0.50];
        parameters.analysis.tfa.plot.plim = [-5 5];
        psr_lfp_plot_tfa(timefreq,parameters);
        
        % Super label
        
        suplabel([                                ...
            'Spectrogram for: '                   ...
            'Subject: ' id_subject        ', '    ...
            'Session: ' id_session        ', '    ...
            'Block: '   id_blocks{iBlock} ', '    ...
            'Probe #: ' num2str(iProbe)],[],[],'none');
        
        % Save the figure
        
        savestr = [savePathFigs ...
            'TFA'                      '_' ...
            id_subject                 '_' ...
            id_session                 '_' ...
            id_blocks{iBlock}          '_' ...
            'P' num2str(iProbe,'%02d')];
        
        export_fig(savestr);
        savefig([savestr '.fig']);
        
    end
end

end

function spikes = qualityThresholding(spikes,threshLow,threshHigh)

% threshLow:  minimum quality of units we accept for the population response
% threshHigh: minimum quality of units we accept for single unit reponse

% Quality 3: isolated single unit
% Quality 2: partial  single unit
% Quality 1: multi-unit
% Quality 0: noise cluster

unitIDs = [spikes.clusters.metrics.id];
unitIDs_high = unitIDs([spikes.clusters.metrics.quality] >= threshHigh); % Sinlge unit response
unitIDs_low  = unitIDs([spikes.clusters.metrics.quality] >= threshLow);  % Population  response

spikesHigh = ismember(spikes.assigns,unitIDs_high);
spikesLow  = ismember(spikes.assigns,unitIDs_low);

spikes.assigns(spikesLow & ~spikesHigh) = 0;
spikes = psr_sst_remove_spikes(spikes,spikesLow,'keep');

end