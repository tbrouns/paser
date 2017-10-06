function PC = ept_pca(spikes,dims,which)

if (nargin < 3); which = 1:size(spikes.waveforms,1); end

PC = pca(spikes.waveforms(which,:)');
if (size(PC,2) > dims)
   PC = PC(:,1:dims); % Take first D principle components
end

end