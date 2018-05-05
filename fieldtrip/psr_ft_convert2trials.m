function data = psr_ft_convert2trials(data,parameters,stimtimes,method)

% Initialize

Fr = parameters.Fr;
stimOnset  = [];
stimOffset = [];
if (nargin < 2); stimtimes = []; end
if (nargin < 3); method    = []; end

% Set onset and offset for every trial

if (strcmp(method,'onset')) % only use onset
    
    onset  = parameters.analysis.lfp.t_win(1);
    offset = parameters.analysis.lfp.t_win(2);
        
    onset     = round(Fr * onset);
    offset    = round(Fr * offset);
    stimtimes = round(Fr * stimtimes(:,1)) + 1;
    
    stimOnset  = stimtimes + onset;
    stimOffset = stimtimes + offset;    
        
elseif (strcmp(method,'interval')) % Use both onset and offset to specify interval
        
    onset = 0;
    stimOnset  = round(Fr * stimtimes(:,1)) + 1;
    stimOffset = round(Fr * stimtimes(:,2)) + 1;
    
end

if (~isempty(stimOnset))
        
    % Adjust trials outside boundary
    
    stimOnset (stimOnset  < data.sampleinfo(1)) = data.sampleinfo(1);
    stimOffset(stimOffset > data.sampleinfo(2)) = data.sampleinfo(2);
    
    cfg = [];
    cfg.trl = [stimOnset,stimOffset,zeros(size(stimOnset))];
    data = ft_redefinetrial(cfg,data);
    
    cfg = [];
    cfg.offset = onset;
    data = ft_redefinetrial(cfg,data);
       
end

end