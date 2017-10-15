function assigns = psr_sst_sorting_OPS(spikes,parameters)

Fs = spikes.params.Fs;
K  = 3; % Number of PCA components to use.

waves = spikes.waveforms(:,:); 
waves = waves'; % [samples per waveform x number of spikes]
data  = spikes.data; 

%% Reduce dimensionality

[U,~,~] = svd(waves,'econ');
U = U(:,1:K);

%% Set parameters:
params.alph    = 1e-1;
params.kappa_0 = 0.01;
params.nu_0    = 0.1;
params.Phi_0   = 0.1*eye(K);
params.a_pii   = 1;
params.b_pii   = 1e7;
params.bet     = 1./(30*Fs);

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