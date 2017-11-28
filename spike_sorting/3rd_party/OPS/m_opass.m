function [z,gam,ngam,muu,lamclus,nu,kappa,Phi,S] = m_opass(x,U,params)
% Runs the M_OPASS algorithm
% x is size N x numCh
% A is size P x K
% of spikes in the data.
%%

apii = params.a_pii; % hyperparameters on the probability of seeing a spike
bpii = params.b_pii; % hyperparameters on the probability of seeing a spike
alph = params.alph;  % parameter of the CRP
Phi0 = params.Phi_0; % prior cluster covariance*nu_0
nu_0 = params.nu_0;  % prior precision of Wishart part of NW distribution

%% Internal Parameters
Cmax      = 50;
curndx    = 0;
lookahead = 500;
rang      = 40;

%%
[N,numCh] = size(x);
[P,K]     = size(U);

%% Calculate precision matrix
[acf] = autocorr(x(:,1),1);
if abs(acf(2)) < 1e-3
    acf(2) = 0;
end
lambi = zeros(P);
for p = 1:P
    lambi(p,:) = 1-p:P-p;
end
sig = acf(2).^abs(lambi)*cov(x(1:1e5,1));
sig(1:(P+1):P^2) = cov(x(1:1e5,1));
lamda = inv(sig);
%%
pii     = apii./bpii;
lamb    = inv(sig);
detlamb = det(lamb);
nu      = repmat(nu_0,Cmax,1);
Phi     = cell(Cmax,numCh);
lamclus = cell(Cmax,numCh);
for ch=1:numCh
    for c=1:Cmax
        Phi{c}{ch}     = Phi0;
        lamclus{c}{ch} = inv(Phi0)*nu_0;
    end
end
kappa_0 = 0.1;
kappa   = kappa_0*ones(Cmax,1);
muu     = zeros(K,numCh,Cmax);
%%
xpad = cell(numCh,1);
for ch = 1:numCh
    xpad{ch}=[x(:,ch);zeros(P,1)];
end
xwind = cell(numCh,1);

%%

C     = 0;
nz    = 0;
z     = zeros(N,1);
gam   = zeros(N,1);
ngam  = zeros(Cmax,1);
S     = zeros(K,numCh,N);
thr   = log(pii./(1-pii));
mT    = N;
sz    = 0;
tlastspike = zeros(Cmax,1);

%%

while curndx < (N - P - rang)
    
    %% set up parameters
    ndx        = (curndx+1:min(mT-P-rang,curndx+lookahead));
    n          = numel(ndx);
    ndxwind    = bsxfun(@plus,ndx,(0:P-1)');
    lthet      = log(ngam./(alph+nz));
    lthet(C+1) = log(alph./(alph+nz));
    for ch = 1:numCh
        xwind{ch} = xpad{ch}(ndxwind);
    end
    
    %% calc llk
    lnone = zeros(1,n,numCh);
    lon   = zeros(C+1,n,numCh);
    for ch = 1:numCh;
        lnone(1,:,ch)=-P/2*log(2*pi)+.5*log(detlamb)-.5*sum((xwind{ch}.*((lamb)*xwind{ch})));
        for c = 1:C
            Q = sig+(1+kappa(c))./kappa(c)*U*(lamclus{c}{ch}\U');
            xwindm = bsxfun(@minus,xwind{ch},U*muu(:,ch,c));
            lon(c,:,ch) = -P/2*log(2*pi)-sum(log(diag(chol(Q))))-.5*sum(xwindm.*(Q\xwindm));
        end
        Q = sig+(1+kappa(C+1))./kappa(C+1)*U*(lamclus{C+1}{ch}\U');
        lon(C+1,:,ch)=-P/2*log(2*pi)-sum(log(diag(chol(Q))))-.5*sum(xwind{ch}.*(Q\xwind{ch}));
    end
    lon      = sum(lon,3);
    lon(:,:) = bsxfun(@plus,lthet(1:C+1,:),lon);
    lnone    = sum(lnone,3);
    H        = bsxfun(@minus,lon,max(lon));
    Hadj     = log(sum(exp(H)));
    lthr     = lnone-max(lon)-Hadj;
    
    %% Find new spike
    Q = find(lthr<thr,1,'first');
    % no spike
    if (numel(Q) == 0) || Q > lookahead-rang
        curndx = curndx + lookahead - rang;
        continue
    end
    % new spike
    [~,offset] = min(lthr(Q:min(Q+rang,numel(lthr))));
    Q     = Q + offset - 1;
    nz    = nz + 1;
    Qt    = Q + curndx;
    z(Qt) = 1;
    [~,Cnew] = max(lon(:,Q));
    if Cnew > C
        C = Cnew;
    end
    tlastspike(Cnew) = Qt;
    ngam(Cnew) = ngam(Cnew)+1;
    gam(Qt)    = Cnew;
    nu(Cnew)   = nu(Cnew)+1;
    for ch = 1:numCh
        curndx = Qt + 1;
        Qmat   = U'*lamda*U+lamclus{Cnew}{ch};
        yhat   = Qmat\(U'*lamda*xwind{ch}(:,Q)+lamclus{Cnew}{ch}*muu(:,ch,Cnew));
        muu(:,ch,Cnew)    = (muu(:,ch,Cnew)*kappa(Cnew)+yhat)./(kappa(Cnew)+1);
        Phi{Cnew}{ch}     = Phi{Cnew}{ch}+kappa(Cnew)./(kappa(Cnew)+1)*(yhat-muu(:,ch,Cnew))*(yhat-muu(:,ch,Cnew))'+inv(Qmat);
        lamclus{Cnew}{ch} = inv(Phi{Cnew}{ch})*nu(Cnew);
        S(:,ch,Qt)        = yhat;
        xpad{ch}(Qt:Qt+P-1) = xpad{ch}(Qt:Qt+P-1)-U*S(:,ch,Qt);
    end
    
    kappa(Cnew) = kappa(Cnew)+1;
    sz = sz+1;
    
    if (C + 1 > Cmax); break; end 
end
