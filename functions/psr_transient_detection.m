function [onset,offset] = psr_transient_detection(signal,threshold)

    % Normalize

    signal = signal / max(signal);
    signal = single(signal); % convert to single type
        
    % Find transients
    
    signalDiff = diff(signal);    
    
    % On- and offsets
    transitions = find(abs(signalDiff) > threshold);
   
    % Join neighbouring transitions
    d = diff(transitions);
    transitions = transitions(find(d ~= 1) + 1);
    
    % Check if we are dealing with onset or offset
    amplitudes = signalDiff(transitions);
    amplitudes = amplitudes ./ abs(amplitudes);
        
    id = find(diff(amplitudes) == -2);
    onset  = transitions(id);
    offset = transitions(id + 1);
    
end