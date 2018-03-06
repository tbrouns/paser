function data = psr_ft_convert2trials(data,parameters,stimtimes,method)

% Initialize

Fr = parameters.Fr;
stimOnset  = [];
stimOffset = [];
if (nargin < 2); stimtimes = []; end
if (nargin < 3); method    = []; end

% Set onset and offset for every trial

if (strcmp(method,'onset')) % only use onset
    
    onset  = parameters.lfp.trial.onset;
    offset = parameters.lfp.trial.offset;
        
    onset     = round(Fr * onset);
    offset    = round(Fr * offset);
    stimtimes = round(Fr * stimtimes(:,1));
    
    stimOnset  = stimtimes + onset;
    stimOffset = stimtimes + offset;    
        
elseif (strcmp(method,'interval')) % Use both onset and offset to specify interval
        
    onset = 0;
    stimOnset  = round(Fr * stimtimes(:,1));
    stimOffset = round(Fr * stimtimes(:,2));
    
end

if (~isempty(stimOnset))
        
    % Ignore trials outside boundary
    
    nlength = round(Fr * data.time{1}(end));
    del = stimOnset < 1 | stimOffset > nlength;
    stimOnset (del) = [];
    stimOffset(del) = [];
    
    cfg = [];
    cfg.trl = [stimOnset,stimOffset,zeros(size(stimOnset))];
    data = ft_redefinetrial(cfg,data);
    
    cfg = [];
    cfg.offset = onset;
    data = ft_redefinetrial(cfg,data);
       
end

end