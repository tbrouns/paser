function waveforms = psr_sst_norm_waveforms(waveforms,normFactor)

% Normalize waveforms 
% normFactor: typically the threshold. Given as [1 x Nchannels] vector
% waveforms: [Nspikes x Npoints x Nchannels]

% Check input
[sz1,sz2,sz3] = size(normFactor);
[~,I] = sort([sz1 sz2 sz3]);
normFactor = permute(normFactor,I); % Move vector to 3rd dimension
if (size(normFactor,3) ~= size(waveforms,3)); return; end

sz = [size(waveforms,1) size(waveforms,2) 1];
waveforms = waveforms ./ repmat(normFactor,sz);

end