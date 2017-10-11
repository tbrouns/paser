function [spikeTimes_1,spikeTimes_2] = ept_stim_detection(data_all,threshold,parameters,Fs)

window_size = (parameters.mfa.p2ptime / 1000) * Fs;  % in samples

nchans = length(data_all);
spikeTimesAll = cell(nchans,1);

for ichan = 1:nchans

    data = data_all{ichan};
    data = data';
    
    % Find peaks and troughs
    [pks_min,loc_min] = findpeaks(double(-data));
    [pks_max,loc_max] = findpeaks(double( data));
    pks_min           = -pks_min;

    % Calculate p2p amplitudes of adjacent peak+trough pair

    Nmin = size(pks_min,2);
    Nmax = size(pks_max,2);

    if (Nmin > Nmax)
        pks_min = pks_min(1:Nmax);
        loc_min = loc_min(1:Nmax);
    elseif (Nmax > Nmin)
        pks_max = pks_max(1:Nmin);
        loc_max = loc_max(1:Nmin);
    end

    if (loc_min(1) < loc_max(1))
        loc_1 = loc_min;
        loc_2 = loc_max;
        pks_1 = pks_min;
        pks_2 = pks_max;
    else
        loc_1 = loc_max;
        loc_2 = loc_min;
        pks_1 = pks_max;
        pks_2 = pks_min;
    end

    % Detect putative artifacts based on amplitudes

    X = [2  1  1];
    x = [0  0 -1];
    Y = fliplr(X);
    y = fliplr(x);

    loc_all = [];

    for i = 1:length(X)

        dn  =   abs(loc_1(X(i):end+x(i))-loc_2(Y(i):end+y(i))); % distance between min and max
        loc = mean([loc_1(X(i):end+x(i));loc_2(Y(i):end+y(i))],1); % between-peak locations

        p1  = pks_1(X(i):end+x(i));
        p2  = pks_2(Y(i):end+y(i));

        if (sum(p1) < sum(p2)); pmin = p1; pmax = p2;
        else,                   pmin = p2; pmax = p1;
        end

        % Filter based on p2p time
        id   = (dn <= window_size);
        pmax = pmax(id);
        pmin = pmin(id);
        loc  =  loc(id);

        % Filter based on max and min amplitudes
        id      = pmax > threshold & pmin < -threshold;
        loc_all = [loc_all, loc(id)]; %#ok, artifact locations

    end

    spikeTimes = loc_all / Fs;
    spikeTimes = sort(spikeTimes);

    periods = diff(spikeTimes);
    periodSorted = sort(periods);
    I = floor(0.5 * length(periods));
    periodLower = median(periodSorted(1:I));
    periodUpper = median(periodSorted(I:end));
    [~,periodIDs] = min([abs(periods - periodLower);abs(periods - periodUpper)]);

    % Extract dual spikes and only take first

    spikeTimesAll{ichan} = spikeTimes(periodIDs == 1);

end

[X,Y] = meshgrid(spikeTimesAll{2},spikeTimesAll{1});
[D,I] = min(abs(X-Y));
I = I(D < window_size / Fs);
id = zeros(length(spikeTimesAll{1}),1);
id(I) = 1; % part of main signal

spikeTimes_1 = spikeTimesAll{1}(find( id)); %#ok
spikeTimes_2 = spikeTimesAll{1}(find(~id)); %#ok

period_1 = mean(diff(spikeTimes_1));
period_2 = mean(diff(spikeTimes_2));

N1 = floor((length(data) / Fs) * (1 / period_1));
N2 = floor((length(data) / Fs) * (1 / period_2));

disp(['Signal 1: detected ' num2str(length(spikeTimes_1)) ' of ' num2str(N1) ' magnetic field artifacts at ' num2str(1 / period_1) ' Hz.']);
disp(['Signal 2: detected ' num2str(length(spikeTimes_2)) ' of ' num2str(N2) ' magnetic field artifacts at ' num2str(1 / period_2) ' Hz.']);

end