function isi = psr_ft_isi(spikesFT,parameters,unitIDs)

% PSR_FT_ISI - Wrapper function for FieldTrip's FT_SPIKE_ISI
% This function calls the FieldTrip function: FT_SPIKE_ISI, which
% calculates inter-spike interval (ISI) data
% 
% Syntax:  isi = psr_ft_isi(spikesFT,parameters,unitIDs)
%
% Inputs:
%    spikesFT   - A FieldTrip spike structure, see:
%                 http://www.fieldtriptoolbox.org/reference/ft_datatype_spike
%    parameters - See PSR_PARAMETERS_ANALYSIS
%    unitIDs    - (Optional) Unit IDs of the units we want to calculate the
%                 ISI metrics for
%
% Outputs:
%    isi - Output from FT_SPIKE_ISI
%
% See also: FT_SPIKE_ISI

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

cfg = [];

if (nargin == 3); cfg.spikechannel = spikesFT.label(unitIDs); end

if (~isempty_field(parameters,'parameters.analysis.isi.bins'));       cfg.bins       = parameters.analysis.isi.bins;       end
if (~isempty_field(parameters,'parameters.analysis.isi.outputunit')); cfg.outputunit = parameters.analysis.isi.outputunit; end
if (~isempty_field(parameters,'parameters.analysis.isi.param'));      cfg.param      = parameters.analysis.isi.param;      end

isi = ft_spike_isi(cfg,spikesFT);
    
end