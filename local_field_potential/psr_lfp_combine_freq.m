function output = psr_lfp_combine_freq(freq)

% Combine multiple FieldTrip LFP structures into one structure

nBlocks = length(freq);
output  = freq{1};
if (isempty(output)); return; end
nDims   = ndims(output.powspctrm);

% Initialize
sz    = zeros(1,nDims);
sz(1) = length(freq);
for iBlock = 1:nBlocks
    for iDim = 2:nDims
        if (isfield(freq{iBlock},'powspctrm'))
            n = size(freq{iBlock}.powspctrm,iDim);
            if (n > sz(iDim)); sz(iDim) = n; end            
        end
    end
end
output.powspctrm = NaN(sz);

% Insert data
itr = 1;
for iBlock = 1:nBlocks
    if (isfield(freq{iBlock},'powspctrm'))
        x = freq{iBlock}.powspctrm;
        sz1 = size(x,    1);
        sz2 = size(x,nDims);
        I = itr : itr + sz1 - 1;
        if (nDims == 4); output.powspctrm(I,:,:,1:sz2) = x;
        else,            output.powspctrm(I,:,  1:sz2) = x;
        end
        itr = itr + sz1;
    end
end

% Change time vector
if (isfield(output,'time'))
    T  = output.time;
    dt = mean(diff(mean(T,1)));
    n  = size(output.powspctrm,nDims);
    t  = T(1):dt:(n-1)*dt+T(1);
    output.time = t;
end

end