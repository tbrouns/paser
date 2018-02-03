function spikesMain = psr_sst_spike_append(spikesMain,spikes)

% Initialize main array

if (~isfield(spikesMain,    'assigns')); spikesMain.assigns     = []; end
if (~isfield(spikesMain,       'data')); spikesMain.data        = []; end
if (~isfield(spikesMain, 'spiketimes')); spikesMain.spiketimes  = []; end
if (~isfield(spikesMain,  'waveforms')); spikesMain.waveforms   = []; end
if (~isfield(spikesMain,       'info')); spikesMain.info        = []; end
if (~isfield(spikesMain.info,   'dur')); spikesMain.info.dur    = []; end
if (~isfield(spikesMain.info,  'stds')); spikesMain.info.stds   = []; end
if (~isfield(spikesMain.info,'thresh')); spikesMain.info.thresh = []; end

T = sum(spikesMain.info.dur);
if (isempty(T)); T = 0; end

if (isfield(spikes,    'assigns')); spikesMain.assigns    = [spikesMain.assigns,    spikes.assigns];        end
if (isfield(spikes,       'data')); spikesMain.data       = [spikesMain.data,       spikes.data];           end
if (isfield(spikes, 'spiketimes')); spikesMain.spiketimes = [spikesMain.spiketimes, spikes.spiketimes + T]; end
if (isfield(spikes,  'waveforms')); spikesMain.waveforms  = [spikesMain.waveforms;  spikes.waveforms];      end

if (isfield(spikes,'info'))
    if (isfield(spikes.info,  'stds')); spikesMain.info.stds   = [spikesMain.info.stds  ;spikes.info.stds];   end
    if (isfield(spikes.info,'thresh')); spikesMain.info.thresh = [spikesMain.info.thresh;spikes.info.thresh]; end
    if (isfield(spikes.info,   'dur')); spikesMain.info.dur    = [spikesMain.info.dur   ;spikes.info.dur];    end
end

end