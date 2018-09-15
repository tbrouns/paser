function output = psr_ft_convert2fieldtrip(input,parameters)

% PSR_FT_CONVERT2FIELDTRIP - Convert PASER structure to FieldTrip structure
%
% Syntax:  output = psr_ft_convert2fieldtrip(input,parameters)
%
% Inputs:
%    input - Can be an LFP or spike structure. 
%            
%            ## In case of LFP data, the input structure should contain the
%               following fields:
% 
%               "data"       : Matrix containing the LFP time series data,
%                              with shape: 
%                              [Number of channels x Number of data points]
%               "timestamps" : Vector containing the timestamp [sec] for
%                              each data point, with shape:
%                              [1 x Number of data points]
%            
%            ## In case of spike data, the input structure should at least
%               contain the following field:
% 
%               "spikes"      : A PASER spike structure (see README)
%            
%               Optionally, it can also contain: 
%  
%               "trialtime"   : A two-column matrix, where the two columns
%                               correspond to relative the on- and offset
%                               of each trial
%               "trialonsets" : Absolute timestamp of [t = 0] for each
%                               trial
%               "popUnit"     : Boolean indicating whether to include a
%                               unit corresponding to the population
%                               activity
%               "probeID"     : Probe number corresponding to the spike
%                               structure
% 
%    parameters - See README 
%
% Outputs:
%    output - Corresponding FieldTrip structure.
%   
%             In case of LFP data, we convert the data to the
%             FT_DATATYPE_RAW format. See FieldTrip documentation:
%             http://www.fieldtriptoolbox.org/reference/ft_datatype_raw
%             http://www.fieldtriptoolbox.org/faq/how_can_i_import_my_own_dataformat
% 
%             In case of spike data, we convert the data to the
%             FT_DATATYPE_SPIKE format. See FieldTrip documentation:
%             http://www.fieldtriptoolbox.org/reference/ft_datatype_spike

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

output = []; % FieldTrip data structure

if (isfield(input,'data') && isfield(input,'timestamps'))
    
    data       = double(input.data);
    timestamps = double(input.timestamps);
        
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
        
    % Parse inputs 
    
    spikes  = input.spikes;
    unitIDs = unique(spikes.assigns);
    probeID = []; 
    popUnit = false; 
    
    if (~isfield(input,'trialtime')   && ~isempty_field(spikes,'spikes.info.trialtime'));  input.trialtime   = spikes.info.trialtime;  end
    if (~isfield(input,'trialonsets') && ~isempty_field(spikes,'spikes.info.trialonset')); input.trialonsets = spikes.info.trialonset; end
    if (~isfield(input,'probeID')     && ~isempty_field(spikes,'spikes.probeID'));         input.probeID     = spikes.probeID;         end

    if (isfield(input,'probeID')); probeID = input.probeID; end
    if (isfield(input,'popUnit')); popUnit = input.popUnit; end
    
    nTrials = size(input.trialtime,1);
    
    % Add a population unit
    
    if (popUnit && all(unitIDs ~= 0)) 
        unitIDs(end+1) = 0;
        unitIDs = sort(unitIDs);
    end
    
    itr = 1;
    for iUnit = unitIDs
        
        if (iUnit == 0); spikeIDs = spikes.assigns >= 0; % Population response
        else,            spikeIDs = ismember(spikes.assigns, iUnit); % Single unit 
        end
                
        % Deal with spikes contained within more than one trial
        
        timestamp = [];
        trialIdxs = [];
        waveforms = [];
        for iTrial = 1:nTrials
            spikeIDs_trial = find(spikes.trials(iTrial,:) & spikeIDs);
            trialIdxs = [trialIdxs,iTrial * ones(size(spikeIDs_trial),'int16')];
            timestamp = [timestamp, spikes.spiketimes(spikeIDs_trial)];
            if (isfield(spikes,'waveforms')); waveforms = [waveforms;spikes.waveforms(spikeIDs_trial,:,:)]; end
        end
        
        time = timestamp;
        for iTrial = unique(trialIdxs)
            spikeIDs_trial = ismember(trialIdxs,iTrial);
            time(spikeIDs_trial) = timestamp(spikeIDs_trial) - input.trialonsets(iTrial);
        end
        
        output.id       {itr} = iUnit;
        output.label    {itr} = ['P' num2str(probeID,'%03d') '-Unit-' num2str(iUnit,'%04d')];
        output.timestamp{itr} = timestamp;
        output.trial    {itr} = trialIdxs;
        output.time     {itr} = time;
        if (~isempty(waveforms)); output.waveforms{itr} = permute(waveforms,[3 2 1]); end
        
        itr = itr + 1;
    end
    
    output.trialtime   = input.trialtime;
    output.trialonsets = input.trialonsets;
end
end