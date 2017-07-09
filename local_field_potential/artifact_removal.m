close all
clear all

files       = dir('raw_FT_*');
files       = char(files.name);
numfiles    = length(files(:,1));

Fresample  = 1000;

acutoff = 2000; % microvolt
zcutoff = 2.0;  % z-value threshold for artifact detection
artipad = 0.05; % artifact padding
Z_SCORE = 0;    % do z-score thresholding

fNoise = 50;   % Frequency [Hz]
aNoise = 100;  % Amplitude

nsections = 10;
nbuffer   =  1; % sec

show_plot = 'no';

bp_lo = 0.5; % hz
bp_hi =  50;

for iTetrode = 1:numfiles
    
    data_input = importdata(files(iTetrode,:));
    
    % Re-sampling
    
    data        = data_input.trial{1};
    timestamps  = data_input.time{1};
    
    Fsample     = data_input.fsample;
    
    if (Fresample < data_input.fsample) % only resample if new frequency is lower
        
        nsamples = data_input.sampleinfo(2);
        
        data_resampled        = [];
        timestamps_resampled  = [];
        
        itr = 1;
        
        for iSection = 1:nsections
            
            nsamples_section = nsamples / nsections;
            
            iStart = itr - nbuffer * Fsample;
            iEnd   = itr + nbuffer * Fsample + nsamples_section;
            
            if (iStart < 1);
                iStart = 1; % turn this into extrapolated buffer
                jStart = 1;
            else
                jStart = nbuffer * Fresample;
            end
            
            if (iEnd   > nsamples);
                iEnd = nsamples; % turn this into extrapolated buffer
                jEnd = 0;
            else
                jEnd = nbuffer * Fresample;
            end
            
            data_section        =     data(:,iStart:iEnd);
            timestamps_section  = timestamps(iStart:iEnd);
            
            [data_section,timestamps_section] = resample(data_section',timestamps_section', Fresample); % resample
            
            data_section        =       data_section(jStart:end-jEnd,:);
            timestamps_section  = timestamps_section(jStart:end-jEnd);
            
            data_resampled        = [      data_resampled;      data_section]; %#ok
            timestamps_resampled  = [timestamps_resampled;timestamps_section]; %#ok
            
            itr = itr + nsamples_section;
            
        end
        
        data_input.trial        = {data_resampled'};
        data_input.time         = {timestamps_resampled'};
        data_input.fsample      = Fresample;
        data_input.sampleinfo   = [1 length(timestamps_resampled)];
        
    end
    
    clear data ...
        data_resampled ...
        data_section ...
        timestamps ...
        timestamps_resampled ...
        timestamps_section
    
    if (~Z_SCORE)
        
        %% Amplitude threshold
        
        asamples  = round(artipad*Fresample);

        signal  = data_input.trial{1};
        artifacts = [];
        for iChan = 1:size(signal,1)
            zscore          = (signal(iChan,:) - mean(signal(iChan,:))) / std(signal(iChan,:));
            stdev           = median(abs(zscore)) / 0.6745;
            acutoff         = 4 * stdev;
            id              = find(abs(zscore) > acutoff);
            artifacts       = [artifacts id id-asamples id+asamples]; %#ok
        end
        
        artifacts = artifacts(artifacts > 0 & artifacts < size(data_input.trial{1},2));
        artifacts = sort(artifacts);
        
        interArtifactInterval = diff(artifacts);
        
        artend   = find(interArtifactInterval > asamples);
        artstart = artend + 1;
        
        artend   = [artend length(interArtifactInterval)+1];    %#ok
        artstart = [1 artstart];                                %#ok
        
        artifact = [artifacts(artstart)', artifacts(artend)'];
        
        artifact(:,1) = sort(artifact(:,1));
        artifact(:,2) = sort(artifact(:,2));
        
        % Remove artifact
        
        cfg = [];
        cfg.artfctdef.reject            = 'nan'; % partial artifact rejection
        cfg.artfctdef.zvalue.artifact   = artifact;
        data_no_artifacts               = ft_rejectartifact(cfg,data_input);
        
        noise  = aNoise*sin(2*pi.*data_no_artifacts.time{1}.*fNoise);
        noise  = repmat(noise,4,1);
        
        [~, col] = find(isnan(data_no_artifacts.trial{1}(1,:)));
        
        x           = data_no_artifacts.trial{1};
        x(:,col)    = (normrnd(mean(nanmean(x,2)),max(nanstd(x,0,2)),size(x,1),length(col)));
        t           = data_no_artifacts.time{1};
        % plot result   
        
        ymax = 5e3;
        figure;
        subplot(1,2,1);
        plot(data_input.time{1},data_input.trial{1});
        ylim([-ymax ymax]);
        subplot(1,2,2);
        plot(t,x);
        ylim([-ymax ymax]);
        pause(0.5);
        
        data_new            = data_no_artifacts;
        data_new.trial{1}   = x;
        data_new.time{1}    = t;
        data_new.sampleinfo = [1 length(t)];
        save(['data_AR_RS_' sprintf('%04d',iTetrode)],'data_new');
        
    else
        
        %% Z-score threshold
        
        % jump
        
        cfg                     = [];
        cfg.trl                 = [];
        cfg.continuous          = 'yes';
        
        % channel selection, cutoff and padding
        cfg.artfctdef.zvalue.channel    = 'all';
        cfg.artfctdef.zvalue.cutoff     = zcutoff;
        cfg.artfctdef.zvalue.trlpadding = 0;
        cfg.artfctdef.zvalue.artpadding = artipad;
        cfg.artfctdef.zvalue.fltpadding = 0;
        
        % algorithmic parameters
        cfg.artfctdef.zvalue.cumulative    = 'no';
        cfg.artfctdef.zvalue.medianfilter  = 'no';
        cfg.artfctdef.zvalue.medianfiltord = 9;
        cfg.artfctdef.zvalue.absdiff       = 'no';
        
        cfg.artfctdef.zvalue.bpfilter   = 'yes';
        cfg.artfctdef.zvalue.bpfilttype = 'but';
        cfg.artfctdef.zvalue.bpfreq     = [bp_lo bp_hi]; % band-pass filtering (Hz)
        cfg.artfctdef.zvalue.bpfiltord  = 4;
        cfg.artfctdef.zvalue.hilbert    = 'yes';
        
        % make the process interactive
        cfg.artfctdef.zvalue.interactive = show_plot;
        
        data_new            = data_input;
        data_new.trial{1}   = abs(data_new.trial{1});
        [cfg, artifact]     = ft_artifact_zvalue(cfg, data_new);
        
        % Remove artifact
        
        cfg = [];
        cfg.artfctdef.reject            = 'nan'; % partial artifact rejection
        cfg.artfctdef.zvalue.artifact   = artifact;
        data_no_artifacts               = ft_rejectartifact(cfg,data_input);
        
        % Removed artifacts
                
        [~, col] = find(isnan(data_no_artifacts.trial{1}(1,:)));
        
        x           = data_no_artifacts.trial{1};
        t           = data_no_artifacts.time{1};
        x(:,col)    = [];
        t(col)      = [];
        
        data_new            = data_no_artifacts;
        data_new.trial{1}   = x;
        data_new.time{1}    = t;
        data_new.sampleinfo = [1 length(t)];
        save(['data_AR_RS_' sprintf('%04d',iTetrode)],'data_new');
        
        ymax = 5e3;
        figure;
        subplot(1,2,1);
        plot(data_input.time{1},data_input.trial{1});
        ylim([-ymax ymax]);
        subplot(1,2,2);
        plot(t,x);
        ylim([-ymax ymax]);
        pause(0.5);
        
    end
end

