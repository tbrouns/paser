function PC = psr_pca(spikes,dims,which,parameters)

if (nargin < 3); which = 1:size(spikes.waveforms,1); end % take all spikes
if (nargin < 4); psr_parameter_default; end

if (isa(spikes.waveforms,'int16'))
    spikes.waveforms = psr_single(spikes.waveforms,parameters);
end

PC = pca(spikes.waveforms(which,:)');
if (size(PC,2) > dims)
   PC = PC(:,1:dims); % Take first D principle components
end

end