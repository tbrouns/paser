function spikesFT_main = psr_ft_append_spikes(spikesFT_main,spikesFT)

if (~isempty(spikesFT_main)); unitIDs_main = cell2mat(spikesFT_main.id);
else,                         unitIDs_main = [];
end

% Initialize optional fields
nTrials = size(spikesFT.trialonsets,1);
if (isempty_field(spikesFT,'spikesFT.stims')); spikesFT.stims = NaN(nTrials,1); end

% Do the concatenation
if (~isempty(spikesFT))
    nUnits = length(spikesFT.label);
    for iUnit = 1:nUnits
        unitID = spikesFT.id{iUnit};
        I = find(unitIDs_main == unitID, 1);
        if (isempty(I)) % Cluster not present
            if (~isempty_field(spikesFT_main,'spikesFT_main.id')); I = length(spikesFT_main.id) + 1;
            else,                                                  I = 1;
            end
            spikesFT_main.id       {I} = unitID;
            spikesFT_main.label    {I} = spikesFT.label    {iUnit};
            spikesFT_main.time     {I} = spikesFT.time     {iUnit};
            spikesFT_main.timestamp{I} = spikesFT.timestamp{iUnit};
            spikesFT_main.trial    {I} = spikesFT.trial    {iUnit};
            spikesFT_main.waveforms{I} = spikesFT.waveforms{iUnit};
        else
            trialMax = size(spikesFT_main.trialtime,1);
            spikesFT_main.time     {I} = cat(2,spikesFT_main.time     {I},spikesFT.time     {iUnit});
            spikesFT_main.timestamp{I} = cat(2,spikesFT_main.timestamp{I},spikesFT.timestamp{iUnit});
            spikesFT_main.trial    {I} = cat(2,spikesFT_main.trial    {I},spikesFT.trial    {iUnit} + trialMax);
            spikesFT_main.waveforms{I} = cat(3,spikesFT_main.waveforms{I},spikesFT.waveforms{iUnit});
        end
    end
    
    if (~isfield(spikesFT_main,'trialtime'));   spikesFT_main.trialtime   = []; end
    if (~isfield(spikesFT_main,'trialonsets')); spikesFT_main.trialonsets = []; end
    if (~isfield(spikesFT_main,'stims'));       spikesFT_main.stims       = []; end
    
    spikesFT_main.trialtime   = cat(1,spikesFT_main.trialtime,  spikesFT.trialtime);
    spikesFT_main.trialonsets = cat(1,spikesFT_main.trialonsets,spikesFT.trialonsets);
    spikesFT_main.stims       = cat(1,spikesFT_main.stims,      spikesFT.stims);
end

end