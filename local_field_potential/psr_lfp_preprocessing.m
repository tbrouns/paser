function freqFilt = psr_lfp_preprocessing(freq,parameters)

% PSR_LFP_PREPROCESSING - Wrapper function for FieldTrip's FT_PREPROCESSING
% This function carries out mains hum removal, zero-phase filtering for
% the local field potential data and then downsampling 
% 
% Syntax:  freqFilt = psr_lfp_preprocessing(freq,parameters)
%
% Inputs:
%    freq       - See README
%    parameters - See README and PSR_PARAMETERS_GENERAL
%
% Outputs:
%    freqFilt - Filtered data
%
% See also: FT_PREPROCESSING, FILTFILT

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

% Set parameters
cfg         = [];
cfg.channel = 'all';
cfg.demean  = 'yes';

switch (parameters.lfp.filter.type)
    case 'bp'
        cfg.bpfilter  = 'yes';
        cfg.bpfreq    = [parameters.lfp.filter.bp.lower,parameters.lfp.filter.bp.upper];
        cfg.bpfiltord =  parameters.lfp.filter.bp.order;
        cfg.bpinstabilityfix = 'split';
    case 'lp'
        cfg.lpfilter  = 'yes';
        cfg.lpfreq    = parameters.lfp.filter.lp.upper;
        cfg.lpfiltord = parameters.lfp.filter.lp.order;
        cfg.lpinstabilityfix = 'split';
end

sLength   = size(freq.trial{1},2);
sPoints   = 60 * parameters.general.twin * parameters.Fs;
nSections = ceil(sLength / sPoints);
sPoints   = ceil(sLength / nSections);
itr       = 1;

freqTemp       = [];
freqTemp.trial = [];
freqTemp.time  = [];

for iSection = 1:nSections
    
    I = itr : itr + sPoints - 1; 
    I(I > sLength) = []; 
    freqSection.label      = freq.label;
    freqSection.fsample    = freq.fsample;
    freqSection.trial      = {freq.trial{1}(:,I)};
    freqSection.time       = {freq.time{ 1}(:,I)};
    freqSection.sampleinfo = [1 length(I)];
    
    itr = I(end) + 1;
    
    if (parameters.lfp.filter.eps.run) % Remove mains hum with notch filter
        try
            % Design filter to remove mains hum
            hw = parameters.lfp.filter.eps.bw / 2;    % Half-width
            lf = parameters.lfp.filter.eps.freq + hw; % Lower frequency
            uf = parameters.lfp.filter.eps.freq - hw; % Upper frequency
            od = parameters.lfp.filter.eps.order;
            filterDesign = designfilt(           ...
                'bandstopiir',                   ...
                'FilterOrder',               od, ...
                'HalfPowerFrequency1',       lf, ...
                'HalfPowerFrequency2',       uf, ...
                'DesignMethod',        'butter', ...
                'SampleRate',     parameters.Fs);
            freqSection.trial = filtfilt(filterDesign,cell2mat(freqSection.trial)')';
            freqSection.trial = {freqSection.trial};
        catch ME
            str_1 = ME.message;
            str_2 = 'Error using FILTFILT. Mains hum not removed.';
            warning = {str_1,str_2};
            psr_show_warning(warning,false,mfilename);
        end
    end
    
    % Add mirror padding
    [freqSection,padding] = addPadding(freqSection,parameters);
    
    % Remove NaNs
    freqSection = psr_ft_nan_removal(freqSection);
    missingChanIDs = any(freqSection.missing{1},2);
    
    try
        freqFilt = ft_preprocessing(cfg,freqSection); % Do FieldTrip pre-processing
        freqFilt = removePadding(freqFilt,padding);
    catch ME
        psr_show_warning({ME.message},true,mfilename);
    end
    
    % Downsample filtered signal
    if (~isempty_field(freqFilt,'freqFilt.trial') && ~isempty_field(parameters,'parameters.lfp.Fr'))
        [dataProbe,timestamps] = resample(cell2mat(freqFilt.trial)',cell2mat(freqFilt.time),parameters.lfp.Fr);
        dataProbe(:,missingChanIDs) = NaN; % Insert NaNs
        freqTemp.trial = [freqTemp.trial,dataProbe'];
        freqTemp.time  = [freqTemp.time, timestamps];
    end
end

freq            = [];
freq.fsample    = parameters.lfp.Fr;
freq.trial      = {freqTemp.trial};
freq.time       = {freqTemp.time};
freq.sampleinfo = [1 length(freqTemp.time)];

end


function [data,padding] = addPadding(data,parameters)

% Mirror padding
Fs      = data.fsample;
nLength = data.sampleinfo(2);
padding = round(parameters.Fs * parameters.lfp.filter.padding);
if (padding > nLength); padding = nLength; end
dt =       1 / Fs;
T  = padding / Fs;

timeFront = data.time{1}(1)   - dt : -dt : data.time{1}(1)   - T;
timeBack  = data.time{1}(end) + dt :  dt : data.time{1}(end) + T;
data.time{1} = [fliplr(timeFront),data.time{1},timeBack];

trialFront = data.trial{1}(:,1:padding);
trialBack  = data.trial{1}(:,end-padding+1:end);
data.trial{1} = [fliplr(trialFront),data.trial{1},fliplr(trialBack)];

data.sampleinfo = [1 size(data.trial{1},2)];

end

function data = removePadding(data,padding)

nLength = size(data.trial{1},2);
n = padding+1:nLength-padding;
data.trial{1} = data.trial{1}(:,n);
data.time{1}  = data.time{1}(n);
data.sampleinfo = [1 nLength];

end