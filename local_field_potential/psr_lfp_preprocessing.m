function data = psr_lfp_preprocessing(dataInput,parameters)

data = [];

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

sLength   = size(dataInput.trial{1},2);
sPoints   = 60 * parameters.general.twin * parameters.Fs;
nSections = ceil(sLength / sPoints);
sPoints   = ceil(sLength / nSections);
itr       = 1;

dataTemp = [];
dataTemp.trial = [];
dataTemp.time  = [];

for iSection = 1:nSections
    
    I = itr:itr+sPoints-1; I(I > sLength) = []; 
    dataSection.label   = dataInput.label;
    dataSection.fsample = dataInput.fsample;
    dataSection.trial   = {dataInput.trial{1}(:,I)};
    dataSection.time    = {dataInput.time{1}(:,I)};
    dataSection.sampleinfo = [1 length(I)];
    
    itr = I(end) + 1;
    
    if (parameters.lfp.filter.eps.run) % Remove mains hum with notch filter
        try
            % Design filter to remove mains hum
            hw = parameters.lfp.filter.eps.bw / 2;    % Half-width
            lf = parameters.lfp.filter.eps.freq + hw; % Lower frequency
            uf = parameters.lfp.filter.eps.freq - hw; % Upper frequency
            od = parameters.lfp.filter.eps.order;
            filterDesign = designfilt('bandstopiir',  ...
                'FilterOrder',         od, ...
                'HalfPowerFrequency1', lf, ...
                'HalfPowerFrequency2', uf, ...
                'DesignMethod',        'butter', ...
                'SampleRate', parameters.Fs);
            dataSection.trial = filtfilt(filterDesign,cell2mat(dataSection.trial)')';
            dataSection.trial = {dataSection.trial};
        catch ME
            disp(ME.message);
            disp('Error using "filtfilt". Mains hum not removed.');
        end
    end
    
    % Add mirror padding
    [dataSection,padding] = addPadding(dataSection,parameters);
    
    % Remove NaNs
    dataSection = psr_ft_nan_removal(dataSection);
    missingChanIDs = any(dataSection.missing{1},2);
    
    try
        data = ft_preprocessing(cfg, dataSection); % Do FieldTrip pre-processing
        data = removePadding(data,padding);
    catch ME
        str = 'FieldTrip ERROR:';
        psr_show_warning({str,ME.message});
    end
    
    % Downsample filtered signal
    if (~psr_isempty_field(data,'data.trial') && ~psr_isempty_field(parameters,'parameters.lfp.Fr'))
        [dataProbe,timestamps] = resample(cell2mat(data.trial)',cell2mat(data.time),parameters.lfp.Fr);
        dataProbe(:,missingChanIDs) = NaN; % Insert NaNs
        dataTemp.trial = [dataTemp.trial,dataProbe'];
        dataTemp.time  = [dataTemp.time, timestamps];
    end
end

data.fsample    = parameters.lfp.Fr;
data.trial      = {dataTemp.trial};
data.time       = {dataTemp.time};
data.sampleinfo = [1 length(dataTemp.time)];

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