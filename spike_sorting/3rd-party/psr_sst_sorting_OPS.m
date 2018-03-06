function spikes = psr_sst_sorting_OPS(spikes,data,parameters)

% NOTE: spikes.assigns not assigned

%% Set parameters

Fs = parameters.Fs;
K  = parameters.sorting.ops.dims; % Number of PCA components to use.

params = parameters.sorting.ops;
params.Phi_0 = params.Phi_0f * eye(K);
params.bet   = 1./(params.betf * Fs);

%% Convert data

data = psr_int16_to_single(data,parameters);
if (size(data,2) > size(data,1)); data = data'; end

%% Extract data

waves = psr_int16_to_single(spikes.waveforms(:,:),parameters); 
waves = waves'; % [samples per waveform x number of spikes]

%% Reduce dimensionality

[U,~,~] = svd(waves,'econ');
U = U(:,1:K);
clear waves;

%% run M_OPASS

[z,gam,ngam,muu,lamclus,nu,kappa,Phi,S] = m_opass(data,U,params);

%% Plot non-trivial clusters
% Plot spikes
C   = max(gam);
col = hsv(C);
figure(1);
clf;
hold on
for c = 1:C
    plot(squeeze(S(1,1,gam==c)),squeeze(S(2,1,gam==c)),'.','Color',col(c,:),'markersize',20)
end
%     plot(S(1,zic>0),S(2,zic>0),'k.','markersize',10);
hold off
xlabel('PCA Component 1','FontSize',16)
ylabel('PCA Component 2','FontSize',16)
title('Inferred y_k Values for Detected Spikes','FontSize',18)
snames = cell(C+1,1);
for c = 1:C
    snames{c} = num2str(c);
end
snames{C+1} = 'IC';
legend(snames);

end