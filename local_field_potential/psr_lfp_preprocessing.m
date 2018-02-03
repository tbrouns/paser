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

% Add mirror padding
[dataInput,padding] = addPadding(dataInput,parameters);

try
    data = ft_preprocessing(cfg, dataInput); % Do FieldTrip pre-processing
    data = removePadding(data,padding);
catch ME
    fprintf('\n **** FieldTrip ERROR **** \n');
    disp(ME.message);
    fprintf(  ' ************************* \n\n');
end

% Downsample filtered signal
if (~isempty(data) && ~isempty(parameters.lfp.Fr))
    [dataProbe,timestamps] = resample(cell2mat(data.trial)',cell2mat(data.time),parameters.lfp.Fr);
    data.trial = {dataProbe'};
    data.time  = {timestamps};
    data.sampleinfo = [1 length(timestamps)];
end

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