function spikesMain = psr_append_spikes(spikesMain,spikes)

% Change to automatic field detection

% Initialize main array

if (psr_isempty_field(spikesMain,'spikesMain.assigns'));         spikesMain.assigns         = []; end
if (psr_isempty_field(spikesMain,'spikesMain.assigns_prior'));   spikesMain.assigns_prior   = []; end
if (psr_isempty_field(spikesMain,'spikesMain.blocks'));          spikesMain.blocks          = []; end
if (psr_isempty_field(spikesMain,'spikesMain.spiketimes'));      spikesMain.spiketimes      = []; end
if (psr_isempty_field(spikesMain,'spikesMain.trials'));          spikesMain.trials          = []; end
if (psr_isempty_field(spikesMain,'spikesMain.waveforms'));       spikesMain.waveforms       = []; end
if (psr_isempty_field(spikesMain,'spikesMain.info.dur'));        spikesMain.info.dur        = []; end
if (psr_isempty_field(spikesMain,'spikesMain.info.trialonset')); spikesMain.info.trialonset = []; end

T = sum(spikesMain.info.dur); if (isempty(T)); T = 0; end
N = max(spikesMain.blocks);   if (isempty(N)); N = 0; end

if (~psr_isempty_field(spikes,'spikes.assigns'));         spikesMain.assigns         = [spikesMain.assigns,         spikes.assigns            ]; end
if (~psr_isempty_field(spikes,'spikes.assigns_prior'));   spikesMain.assigns_prior   = [spikesMain.assigns_prior,   spikes.assigns_prior      ]; end
if (~psr_isempty_field(spikes,'spikes.blocks'));          spikesMain.blocks          = [spikesMain.blocks,          spikes.blocks          + N]; end
if (~psr_isempty_field(spikes,'spikes.spiketimes'));      spikesMain.spiketimes      = [spikesMain.spiketimes,      spikes.spiketimes      + T]; end
if (~psr_isempty_field(spikes,'spikes.trials'));          spikesMain.trials          = [spikesMain.trials,          spikes.trials             ]; end
if (~psr_isempty_field(spikes,'spikes.waveforms'));       spikesMain.waveforms       = [spikesMain.waveforms;       spikes.waveforms          ]; end
if (~psr_isempty_field(spikes,'spikes.info.dur'));        spikesMain.info.dur        = [spikesMain.info.dur;        spikes.info.dur;          ]; end
if (~psr_isempty_field(spikes,'spikes.info.trialonset')); spikesMain.info.trialonset = [spikesMain.info.trialonset; spikes.info.trialonset + T]; end

end