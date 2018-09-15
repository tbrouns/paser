function freq = psr_ft_convert2trials(freq,parameters,stimtimes,method)

% PSR_FT_CONVERT2TRIALS - Segment FieldTrip data into trials 
% This function uses input stimulus timings to divide the FieldTrip data
% structure for LFP data into trials aligned to each stimulus onset
% 
% Syntax:  freq = psr_ft_convert2trials(freq,parameters,stimtimes,method)
%
% Inputs:
%    freq       - FieldTrip data structure for LFP data
%    parameters - See PSR_PARAMETERS_ANALYSIS
%    stimtimes  - Array containing stimulus timings [sec]:
% 
%                 ## If method = 'onset', 
%                    then stimtimes should be a column vector containing
%                    each stimulus onset
% 
%                 ## If method = 'interval',
%                    then stimtimes should be a two-column matrix
%                    containing the stimulus onsets in the first column and
%                    the stimulus offsets in the second column
% 
%    method     - Either 'onset' or 'interval' (see above)
% 
% Outputs:
%    freq - Segmented data
%
% See also: PSR_STIMULUS_WINDOW, FT_REDEFINETRIAL, PSR_FT_CONVERT2FIELDTRIP

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

% Initialize

Fr = freq.fsample;
stimOnset  = [];
stimOffset = [];
if (nargin < 2); stimtimes = []; end
if (nargin < 3); method    = []; end

% Set onset and offset for every trial

if (strcmp(method,'onset')) % only use onset
    
    tPre  = parameters.analysis.lfp.t_win(1);
    tPost = parameters.analysis.lfp.t_win(2);
            
    stimOnset     = stimtimes(:,1) + tPre;
    stimOffset    = stimtimes(:,1) + tPost;
    stimTimes     = [stimOnset,stimOffset];
    stimDurations = zeros(size(stimtimes,1));
    stimOnset     = round(Fr * stimOnset)  + 1;
    stimOffset    = round(Fr * stimOffset) + 1;
    sPre          = round(Fr * tPre);
    
elseif (strcmp(method,'interval')) % Use both onset and offset to specify interval
        
    sPre = 0;
    stimOnset     = stimtimes(:,1);
    stimOffset    = stimtimes(:,2);
    stimTimes     = [stimOnset,stimOffset];
    stimDurations = stimOffset - stimOnset;
    stimOnset     = round(Fr * stimOnset)  + 1;
    stimOffset    = round(Fr * stimOffset) + 1;
    
end

if (~isempty(stimOnset))
        
    % Adjust trials outside boundary
    
    stimOnset (stimOnset  < freq.sampleinfo(1)) = freq.sampleinfo(1);
    stimOffset(stimOffset > freq.sampleinfo(2)) = freq.sampleinfo(2);
    
    cfg = [];
    cfg.trl = [stimOnset,stimOffset,zeros(size(stimOnset))];
    freq = ft_redefinetrial(cfg,freq);
    
    cfg = [];
    cfg.offset = sPre;
    freq = ft_redefinetrial(cfg,freq);
    
    % Save additional metrics
    
    freq.cfg.trialtime = stimTimes; 
    freq.cfg.stimdur   = stimDurations;
    
end

end