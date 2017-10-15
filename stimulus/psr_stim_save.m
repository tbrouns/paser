function psr_mfa_save(files,MFAtimes)

ntets = length(files);

% Save MFA times
    
for iTetrode = 1:ntets
    load(files{iTetrode});
    spikes.stimtimes = MFAtimes; %#ok
    save(files{iTetrode},'spikes','freq','parameters');
end

end