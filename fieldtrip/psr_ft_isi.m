function isi = psr_ft_isi(spikesFT,parameters,clustIDs)

cfg = [];

if (nargin == 3); cfg.spikechannel = spikesFT.label(clustIDs); end

if (~isempty_field(parameters,'parameters.analysis.isi.bins'));       cfg.bins       = parameters.analysis.isi.bins;       end
if (~isempty_field(parameters,'parameters.analysis.isi.outputunit')); cfg.outputunit = parameters.analysis.isi.outputunit; end
if (~isempty_field(parameters,'parameters.analysis.isi.param'));      cfg.param      = parameters.analysis.isi.param;      end

isi = ft_spike_isi(cfg,spikesFT);
    
end