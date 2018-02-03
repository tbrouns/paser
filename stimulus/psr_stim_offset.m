function tPks = psr_stim_offset(data,stimTimes,Fs)

tWinSlope = 0.1; % ms
tWinOff   = 10;  % ms
tWinStim  = 0.1; % sec

sWinSlope  =  (tWinSlope / 1000) * Fs; 
sWinOff    =  (tWinOff   / 1000) * Fs; 
sWinStim   =  -tWinStim  * Fs;

tPks = [];

data = data(1+sWinSlope:end) - data(1:end-sWinSlope);

% Extract stimulus windows

Nmax = size(data,1);
stimTimes = round(Fs * stimTimes{1}(:,1) + 1);
if (size(stimTimes,1) > size(stimTimes,2)); stimTimes = stimTimes'; end

ids = bsxfun(@plus,stimTimes,(sWinStim:0)');
ids(ids < 1)    = 1;
ids(ids > Nmax) = Nmax;

data = abs(data(ids));
data = mean(data,2);

[pks,loc] = findpeaks(data);
[~,I] = sort(pks,'descend');
loc   = loc(I);
npks  = length(loc);

if (npks >= 2)
    
    locMax = loc(1);
    itr = 1;
    d   = 0;
    
    while (d < sWinOff && itr < npks)
        itr = itr + 1;
        d = abs(locMax - loc(itr));
    end
    
    if (itr <= npks)
        loc  = loc([1,itr]);
        tPks = (loc - 1) / Fs;
        tPks = tPks - tWinStim;
    end
end

end