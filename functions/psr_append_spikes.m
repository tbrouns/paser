function spikesMain = psr_append_spikes(spikesMain,spikes)

% Change to automatic field detection

% Initialize main array

if (isempty_field(spikesMain,'spikesMain.assigns'));         spikesMain.assigns         = []; end
if (isempty_field(spikesMain,'spikesMain.assigns_prior'));   spikesMain.assigns_prior   = []; end
if (isempty_field(spikesMain,'spikesMain.blocks'));          spikesMain.blocks          = []; end
if (isempty_field(spikesMain,'spikesMain.spiketimes'));      spikesMain.spiketimes      = []; end
if (isempty_field(spikesMain,'spikesMain.trials'));          spikesMain.trials          = logical([]); end
if (isempty_field(spikesMain,'spikesMain.waveforms'));       spikesMain.waveforms       = []; end
if (isempty_field(spikesMain,'spikesMain.info.dur'));        spikesMain.info.dur        = []; end
if (isempty_field(spikesMain,'spikesMain.info.trialonset')); spikesMain.info.trialonset = []; end

T = sum(spikesMain.info.dur); if (isempty(T)); T = 0; end
N = max(spikesMain.blocks);   if (isempty(N)); N = 0; end

if (~isempty_field(spikes,'spikes.assigns'));         spikesMain.assigns         = [spikesMain.assigns,         spikes.assigns            ]; end
if (~isempty_field(spikes,'spikes.assigns_prior'));   spikesMain.assigns_prior   = [spikesMain.assigns_prior,   spikes.assigns_prior      ]; end
if (~isempty_field(spikes,'spikes.blocks'));          spikesMain.blocks          = [spikesMain.blocks,          spikes.blocks          + N]; end
if (~isempty_field(spikes,'spikes.spiketimes'));      spikesMain.spiketimes      = [spikesMain.spiketimes,      spikes.spiketimes      + T]; end
if (~isempty_field(spikes,'spikes.waveforms'));       spikesMain.waveforms       = [spikesMain.waveforms;       spikes.waveforms          ]; end
if (~isempty_field(spikes,'spikes.info.dur'));        spikesMain.info.dur        = [spikesMain.info.dur;        spikes.info.dur;          ]; end
if (~isempty_field(spikes,'spikes.info.trialonset')); spikesMain.info.trialonset = [spikesMain.info.trialonset; spikes.info.trialonset + T]; end

if (~isempty_field(spikes,'spikes.trials'))         
    nTrialsMain = size(spikesMain.trials,1);
    nSpikesMain = size(spikesMain.trials,2);
    nTrialsNew  = size(spikes.trials,1);
    nSpikesNew  = size(spikes.trials,2);
    spikesMain.trials = [spikesMain.trials,false(nTrialsMain,nSpikesNew)]; 
    spikes.trials     = [false(nTrialsNew,nSpikesMain),spikes.trials];
    spikesMain.trials = [spikesMain.trials;spikes.trials];
end

end