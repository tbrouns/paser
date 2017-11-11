function [data,artifacts] = psr_lfp_artifact_removal(data_raw,parameters)

Fs = parameters.Fr; 
Fr = parameters.lfp.artifact_freq;

zscore = double(mean(data_raw));
zscore = zscore - mean(zscore) / std(zscore);
zscore = abs(zscore);

t      = (1:size(zscore,2)) / Fs;
data   = gaussfilt(t,zscore,1/Fr);

% dZ   = abs(diff(zscore));
% dV   = abs(diff(data));

[pks_max,locs_max_all] = findpeaks(double(zscore)); % detect peaks in raw data signal
[      ~,locs_min_all] = findpeaks(double( -data)); % detect on- and offset points in resampled data

tsection  = parameters.lfp.artifact_tsection;
threshold = parameters.lfp.artifact_thresh;
acutoff   = findThreshold(zscore,threshold,tsection,Fs);

% Find peaks above threshold

locs_max = locs_max_all(pks_max > acutoff);

% Find adjacent minima

nsamples  = size(data,2);
artifacts = findRange(locs_max,locs_min_all,nsamples);
artifacts = artifacts';
artifacts = unique(artifacts,'rows');

% % Find large derivatives
% 
% locs_new = find(dV < parameters.lfp.artifact_frac * std(dV));
% locs_all = [locs_max_all,locs_min_all,locs_new];
% locs_all = sort(locs_all);
% 
% acutoff = findThreshold(dZ,threshold,tsection,Fs);
% [pks,locs] = findpeaks(double(dZ));
% 
% locs = locs(pks > acutoff);
% 
% artifacts_2 = findRange(locs,locs_all,nsamples);
% artifacts_2 = unique(artifacts_2','rows');
% artifacts   = [artifacts_1;artifacts_2];

% Remove artifacts from data

for iArtifact = 1:size(artifacts,1)
    x = artifacts(iArtifact,1):artifacts(iArtifact,2);
    data_raw(:,x) = NaN;
end

data      = data_raw;
artifacts = artifacts / Fs;

end

function acutoff = findThreshold(data,threshold,tsection,fs)

% Find background noise

nlength  = size(data,2);
nsamples = round(tsection * fs); % cut data in sections of X seconds
nstep    = round(0.1 * (nsamples)); % move window with steps of 0.1*nsection
stdev    = [];

iStart = 1;
iEnd   = iStart + nsamples;
STOP   = 0;

while (~STOP)
    
    if (iEnd > nlength); 
        iStart = nlength - nsamples; 
        iEnd   = nlength;
        STOP   = 1;
    end
    
    if (  iEnd > nlength); iEnd   = nlength; end
    if (iStart <       1); iStart = 1;       end
    
    data_section = data(:,iStart:iEnd);
    stdev  = [stdev;median(abs(data_section)) / 0.6745]; %#ok
    iStart = iStart + nstep;
    iEnd   = iStart + nsamples;
end
    
stdev   = min(stdev);
acutoff = threshold * stdev;

end

function locs_artifact = findRange(locs,locs_all,nsamples)

N             = length(locs);
locs_artifact = zeros(2,N);

for iLoc = 1:N
    loc_max = locs(iLoc);
    id = (loc_max - locs_all) > 0;
    [~,I1] = min(loc_max - locs_all( id));
    [~,I2] = max(loc_max - locs_all(~id));
    n = find(~id,1) - 1;
    if (~isempty(I1)); locs_artifact(1,iLoc) = locs_all(I1);
    else               locs_artifact(1,iLoc) = 1;
    end
    if (~isempty(I2)); locs_artifact(2,iLoc) = locs_all(n + I2);
    else               locs_artifact(2,iLoc) = nsamples;
    end
end

end