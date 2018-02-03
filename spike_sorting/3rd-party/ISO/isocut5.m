function [dip_score,cutpoint,info] = isocut5(samples,sample_weights,opts)

if (nargin < 2); sample_weights = []; end
if (nargin < 3); opts = struct;       end

if (~isfield(opts,'already_sorted'));  opts.already_sorted  =  0; end
if (~isfield(opts,'return_info'));     opts.return_info     = []; end
if (~isfield(opts,'num_bins_factor')); opts.num_bins_factor =  1; end

info = struct;

[~,N] = size(samples);
if (N == 0)
    error('Error in isocut5: N is zero.');
end;
num_bins_factor = opts.num_bins_factor;
num_bins        = ceil(sqrt(N/2) * num_bins_factor);

if (isempty(sample_weights)), sample_weights=ones(1,N); end;

if (opts.already_sorted)
    X = samples;
else
    [X,sort_inds] = sort(samples);
    sample_weights = sample_weights(sort_inds);
end;

while 1
    num_bins_1 = ceil(num_bins/2);
    num_bins_2 = num_bins-num_bins_1;
    intervals  = [1:num_bins_1,num_bins_2:-1:1];
    alpha      = (N-1)/sum(intervals);
    intervals  = intervals*alpha;
    inds       = floor([1,1+cumsum(intervals)]);
    N_sub      = length(inds);
    if (min(intervals) >= 1)
        break;
    else
        num_bins = num_bins - 1;
    end;
end;

cumsum_sample_weights=cumsum(sample_weights);

X_sub          = X(inds);
spacings       = X_sub(2:end) - X_sub(1:end-1);
multiplicities = cumsum_sample_weights(inds(2:end)) - cumsum_sample_weights(inds(1:end-1));
densities      = multiplicities ./ spacings;

densities_unimodal_fit=jisotonic5(densities,'updown',multiplicities);

[~,peak_density_ind] = max(densities_unimodal_fit); 
peak_density_ind = peak_density_ind(1);
[ks_left, ks_left_index ] = compute_ks5(multiplicities(1:peak_density_ind),     densities_unimodal_fit(1:peak_density_ind)      .* spacings(1:peak_density_ind));
[ks_right,ks_right_index] = compute_ks5(multiplicities(end:-1:peak_density_ind),densities_unimodal_fit(end:-1:peak_density_ind) .* spacings(end:-1:peak_density_ind));
ks_right_index = length(spacings) - ks_right_index + 1;

if (ks_left > ks_right)
    critical_range = 1:ks_left_index;
    dip_score = ks_left;
else
    critical_range = ks_right_index:length(spacings);
    dip_score = ks_right;
end;

densities_resid     = densities - densities_unimodal_fit;
densities_resid_fit = jisotonic5(densities_resid(critical_range),'downup',spacings(critical_range));
[~,cutpoint_ind]    = min(densities_resid_fit); 
cutpoint_ind        = cutpoint_ind(1);
cutpoint_ind        = critical_range(1) + cutpoint_ind - 1;
cutpoint            = (X_sub(cutpoint_ind) + X_sub(cutpoint_ind+1)) / 2;

if (opts.return_info)
    info.spacings           = spacings;
    info.lefts              = X_sub(1:end-1);
    info.rights             = X_sub(2:end);
    info.centers            = (info.lefts+info.rights)/2;
    info.densities          = densities;
    info.densities_unimodal = densities_unimodal_fit;
    info.critical_range     = critical_range;
    info.densities_bimodal  = densities_resid_fit+densities_unimodal_fit(critical_range);
    info.plot_xx            = zeros(1,(N_sub-1)*2);
    info.plot_xx(1:2:end)   = info.lefts;
    info.plot_xx(2:2:end)   = info.rights;
    info.plot_densities     = zeros(1,(N_sub-1)*2);
    info.plot_densities(1:2:end) = info.densities;
    info.plot_densities(2:2:end) = info.densities;
    info.plot_densities_unimodal = zeros(1,(N_sub-1)*2);
    info.plot_densities_unimodal(1:2:end) = info.densities_unimodal;
    info.plot_densities_unimodal(2:2:end) = info.densities_unimodal;
    [counts,info.hist_bins] = hist_with_weights(samples,sample_weights,num_bins*3);
    info.hist_counts        = counts;
    info.hist_densities     = counts/(info.hist_bins(2)-info.hist_bins(1));
else
    info = struct;
end;

function [ks] = compute_ks4(counts1,counts2)
S1 = cumsum(counts1) / sum(counts1);
S2 = cumsum(counts2) / sum(counts2);
ks = max(abs(S1 - S2));
ks = ks * sqrt((sum(counts1) + sum(counts2))/2);

function [best_ks,best_len] = compute_ks5(counts1,counts2)
len     = length(counts1);
best_ks = -inf;
while (len >= 4)|| (len == length(counts1))
    ks = compute_ks4(counts1(1:len),counts2(1:len));
    if (ks > best_ks)
        best_ks  = ks;
        best_len = len;
    end;
    len = floor(len/2);
end

function [counts,bins] = hist_with_weights(X,weights,num_bins)
bin_width = (max(X(:)) - min(X(:))) / num_bins;
bin_ints  = round(X / bin_width);
i1        = min(bin_ints);
i2        = max(bin_ints);
bins      = (i1:i2) * bin_width;
ii        = bin_ints - i1 + 1;
counts    = accumarray(ii',weights',[length(bins),1]);