function data = psr_ft_convert2trials(data,parameters,stimtimes,method)

% Initialize

Fr = data.fsample;
stimOnset  = [];
stimOffset = [];
if (nargin < 2); stimtimes = []; end
if (nargin < 3); method    = []; end

% Set onset and offset for every trial

if (strcmp(method,'onset')) % only use onset
    
    tPre  = parameters.analysis.lfp.t_win(1);
    tPost = parameters.analysis.lfp.t_win(2);
            
    stimOnset     = stimtimes(:,1) + tPre;
    stimOffset    = stimtimes(:,1) + tPost;
    stimTimes     = [stimOnset,stimOffset];
    stimDurations = zeros(size(stimtimes,1));
    stimOnset     = round(Fr * stimOnset)  + 1;
    stimOffset    = round(Fr * stimOffset) + 1;
    sPre          = round(Fr * tPre);
    
elseif (strcmp(method,'interval')) % Use both onset and offset to specify interval
        
    sPre = 0;
    stimOnset     = stimtimes(:,1);
    stimOffset    = stimtimes(:,2);
    stimTimes     = [stimOnset,stimOffset];
    stimDurations = stimOffset - stimOnset;
    stimOnset     = round(Fr * stimOnset)  + 1;
    stimOffset    = round(Fr * stimOffset) + 1;
    
end

if (~isempty(stimOnset))
        
    % Adjust trials outside boundary
    
    stimOnset (stimOnset  < data.sampleinfo(1)) = data.sampleinfo(1);
    stimOffset(stimOffset > data.sampleinfo(2)) = data.sampleinfo(2);
    
    cfg = [];
    cfg.trl = [stimOnset,stimOffset,zeros(size(stimOnset))];
    data = ft_redefinetrial(cfg,data);
    
    cfg = [];
    cfg.offset = sPre;
    data = ft_redefinetrial(cfg,data);
    
    % Save additional metrics
    
    data.cfg.trialtime = stimTimes; 
    data.cfg.stimdur   = stimDurations;
    
end

end