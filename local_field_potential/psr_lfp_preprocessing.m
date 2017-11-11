function [data,parameters] = psr_lfp_preprocessing(data_input,parameters)

% Filtering sometimes returns following error: "Calculated filter
% coefficients have poles on or outside the unit circle and will not be
% stable. Try a higher cutoff frequency or a different type/order of
% filter." Therefore we must try different filter orders and cut-off
% frequencies until one works. Sub-sampling data also works.

data = [];

cfg                  = [];
cfg.channel          = 'all';
cfg.continuous       = 'yes';
cfg.demean           = 'yes';
cfg.preproc.bpfilter = 'yes';

orders = parameters.lfp.bp_order:-1:1;                           % Try different filter orders (min 3)
step   = 0.5 * parameters.lfp.bp_lower;
F_low  = parameters.lfp.bp_lower:step:parameters.lfp.freq_lower; % Try different frequencies (in Hz)

SUCCESS = 0;

for i = 1:length(orders)
    
    if (SUCCESS); break; end
    
    for j = 1:length(F_low)
        
        cfg.preproc.bpfreq     = [F_low(j), parameters.lfp.bp_upper];
        cfg.preproc.bpfiltord  = orders(i);
        
        try
            data  = ft_preprocessing(cfg, data_input);
            SUCCESS = 1;
            
            parameters.lfp.bp_order = cfg.preproc.bpfiltord;
            parameters.lfp.bp_lower = cfg.preproc.bpfreq;
            
            break;
        catch ME
            disp(ME);
            SUCCESS = 0;
        end
    end
end

end