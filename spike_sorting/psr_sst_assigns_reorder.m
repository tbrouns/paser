function spikes = psr_sst_assigns_reorder(spikes)
    assigns = spikes.assigns;
    assignsIDs = unique(assigns);
    nIDs = length(assignsIDs);
    for id = 1:nIDs
        spikeIDs = assigns == assignsIDs(id);
        assigns(spikeIDs) = id;
    end
    spikes.assigns = assigns;
end