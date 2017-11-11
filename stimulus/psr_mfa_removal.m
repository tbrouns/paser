function psr_mfa_removal(spikes,parameters)

    % Add magnetic field artifacts
    
    nCells = length(MFAtimesAll);
    MFAtimes = [];
    for iCell = 1:nCells; MFAtimes = [MFAtimes,MFAtimesAll{iCell}]; end %#ok
    MFAtimes = MFAtimes';
    
    if (~isempty(MFAtimes))
        MFAs = round(Fs * MFAtimes);
        MFAs(MFAs <= sWindowSpike) = sWindowSpike + 1;
        MFAs(MFAs > size(data,2) - sWindowSpike) = size(data,2) - sWindowSpike;
        MFAs = bsxfun(@plus,MFAs,-sWindowSpike:sWindowSpike)';
        MFAs = MFAs(:);
        artifacts = [artifacts; MFAs];
    end

end