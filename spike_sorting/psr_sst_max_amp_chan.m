function chanMaxID = psr_sst_max_amp_chan(spikes,clustID,parameters)

% Find maximum amplitude channel with more weight to waveforms that peak in
% the centre of the window and also with an absolute threshold

nPoints = size(spikes.waveforms,2);
[~,ampAbs,ampRel,~,peakLocs] = psr_sst_cluster_amp(spikes,clustID,parameters);

% Weight factor 
nhalf   = 0.5 * nPoints;
ds      = abs(peakLocs - nhalf);
weights = 1 - (ds / nhalf);
ampRel  = weights .* ampRel;
ampAbs  = weights .* ampAbs;

ignore  = ampAbs <= parameters.spikes.min_amp; % Needs to exceed absolute threshold
if (~all(ignore)); ampRel(ignore) = NaN; end % Ignore selected channels
[~, chanMaxID] = max(ampRel);

end