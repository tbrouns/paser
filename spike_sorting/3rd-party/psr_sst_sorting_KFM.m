function assigns = psr_sst_sorting_KFM(spikes,parameters)

% % Define parameters
% 
dims     = parameters.sorting.kfm.dims;
J        = parameters.sorting.kfm.nC; % number of clusters
em_iters = 20;   % number of Em iterations for the GMM
EM_iters = 20;   % number of EM iterations for the KFMM

inspk = psr_sst_wavelet_features(spikes,parameters);
inspk = inspk';

% spiketimes = spikes.spiketimes; 
% T = round(parameters.Fs * spiketimes(end));
% spiketimes = round(parameters.Fs * spiketimes);
% obs_id             = false(T,1);
% obs_id(spiketimes) = true;
% Y = realmax('single') * ones(T, dims, 'single');
% Y(spiketimes,:) = single(inspk);

T = size(spikes.waveforms,1);
obs_id = true(T,1);
Y = inspk; clear inspk;

% Structure to store the results
P = struct([]); P_new = struct([]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameter initialization: use a Gaussian Mixture Model (MoG) to make a first
% guess of the clusters' ids.
disp('Initializing : running MoG ...');

[initP, initp] = GaussianMixtureModel(Y, J, obs_id, em_iters);

Cu = zeros(dims,dims);
Cu(1:dims+1:end) = 1e-3;

for j = 1 : J
    P(j).Q = Cu;                        % system covariance
    for t = 1 : T
        P(j).x(:,t) = initP(j).u;
        P(j).V(:, :, t) = initP(j).Cv;  % initial guess for the state covariance
        P(j).R = initP(j).Cv;           % observation covariance
        P(j).a = 1/J;                   % clusters weights
    end
end

% Auxiliary structures
for j = 1 : J
    P(j).xs  = zeros(dims, T);     % smoothed positions
    P(j).Vs  = zeros(dims, dims, T); % smoothed covariances
end

% initial guess for the responsibilities
p     = initp;                       % probability that the t-th observation belongs to the j-th cluster.
% Initially we guess these using the MoG implementation


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EM RECURSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('running MoK ...');
%conv_KF = [];

for iter = 1 : EM_iters
    
    % Assign cluster ids
    cl_id = zeros(size(obs_id));
    for k = 1 : T
        if obs_id(k)
            [~, I] = max(p(:,k));
            for j = 1 : J
                if (I == j) ~= 0
                    cl_id(k) = j;
                end
            end
        else
            cl_id(k) = -1;
        end
    end
    
    %cl_idKF = cl_id;
    
    % Forward step for updating the means
    for j = 1 : J
        for t = 2 : T
            % if there is an observation at t-1 ...
            if obs_id(t-1) == 1
                P(j).V(:, :, t) = inv(inv(P(j).V(:, :, t - 1) + P(j).Q) + (p(j, t - 1) * inv(P(j).R)));
                P(j).x(:, t) = P(j).V(:, :, t) * (inv(P(j).V(:, :, t - 1) ...
                    + P(j).Q) * P(j).x(:, t - 1) + p(j, t - 1) * (inv(P(j).R) * Y(t - 1, :)'));
            else
                % if there is NO observation ...
                P(j).V(:, :, t) = inv(inv(P(j).V(:, :, t - 1) + P(j).Q));
                P(j).x(:, t) = P(j).V(:, :, t) * (inv(P(j).V(:, :, t - 1) + P(j).Q) * P(j).x(:, t - 1));
            end
        end
    end
    
    % Backward step for updating the means
    for j = 1 : J
        P(j).xs(:, T) = P(j).x(:, T);
        P(j).Vs(:, :, T) = P(j).V(:, :, T);
        for t = T - 1 : -1 : 1
            K = P(j).V(:, :, t) * inv(P(j).V(:, :, t) + P(j).Q);
            P(j).xs(:, t) = P(j).x(:, t) + K * (P(j).xs(:, t + 1) - P(j).x(:, t));
            P(j).Vs(:, :, t) = P(j).V(:, :, t) + K * (P(j).Vs(:, :, t + 1) - (P(j).V(:, :, t) + P(j).Q)) * K';
        end
        for t = T  : -1 : 1
            P(j).x(:, t) = P(j).xs(:, t);
            P(j).V(:, :, t) = P(j).Vs(:, :, t);
        end
    end
    P(2).xs = P(2).x;
    
    % Update observation covariance
    for j = 1 : J
        P_new(j).R = 0 * P(j).R;
        for t = 1 : T
            P_new(j).R = P_new(j).R + p(j, t) * (Y(t, :)' - P(j).x(:, t)) * (Y(t, :) - P(j).x(:, t)');
        end
        P_new(j).R = P_new(j).R / sum(p(j, 1 : T));
    end
    for j = 1 : J
        P(j).R = P_new(j).R;
    end
    
    % Update state covariance
    for j = 1 : J
        for t = 1 : T
            P(j).V(:, :, t) = P(j).R;
        end
    end
    
    % Estimate probabilities
    for t = 1 : T
        if obs_id(t) == 1
            % if there is an observation at t ...
            normalization = 0;
            for j = 1 : J
                p(j, t) = exp(- 0.5 * (log(det(P(j).R + P(j).Q)) + (Y(t, :) ...
                    - P(j).x(:, t)') * inv(P(j).R + P(j).Q) * (Y(t, :)' - P(j).x(:, t))));
                normalization = normalization + p(j, t);
            end
            p(:, t) = p(:, t) / normalization;
        else
            % if there is no observation at t...
            p(:, t) = [0, 0];
        end
    end
    
end

[~,assigns] = max(p,[],2);

end

function [P, p, cl_id] = GaussianMixtureModel(V, J, obs_id, EMit)
% This functions assigns cluster ids to the simulated data 'DataSet' using a 
% mixture of Gaussians model. 

    spkt = obs_id;
    T = length(V(:,1));         % experiment length
    D = length(V(1,:));         % dimensionality of data
    Tobs = sum(spkt);
        
    % Parameter Initialization: use simple a k-means clustering algorithm to
    % determine the centers u1,...,uJ of J components. Set a1,..,aJ = 1/J;
    % Cv_1,...,Cv_J = eye(D). 
    
    % Use only the times that have an observation for the kmeans algorithm
    i = 0; Nobs = sum(obs_id==1); Vobs = zeros(Nobs, D);
    for t = 1 : T
        if obs_id(t) == 1
            i = i+1;
            Vobs(i, :) = V(t, :);          
        end
    end
    
    [~, U] = kmeans(Vobs, J);
    
    % initial guesses for neuron cluster
    P = struct([]);   
    for j = 1 : J
        P(j).a = 1/J;       % clusters' probabilities
        P(j).u = U(j,:);    % clusters means
        P(j).Cv = eye(D);   % observation noise
    end
    
    % auxiliary structures
    w = zeros(J, T);
    p = zeros(J, T);   
    
    for idx = 1 : EMit
        
        % Compute a first guess of pj(i) = p(z(i)=j|guessed parameters), where z 
        % is the identity of the ith waveform and j represents a particular cluster. Then we 
        % have p(j,i) = p(z(i)=j|guessed parameters).
        for i = 1 : T
            if obs_id(i) == 1
                normalization = 0;
                mm = zeros(1,J);
                for j = 1 : J
                    mm(j) = log(P(j).a) - 0.5 * (log(det(P(j).Cv)) + (V(i, :) - P(j).u) * inv(P(j).Cv) * (V(i, :) - P(j).u)');
                end
                    m = max(mm);
                for j = 1 : J
                    w(j, i) = exp(log(P(j).a) - 0.5 * (log(det(P(j).Cv)) + (V(i, :) - P(j).u) * inv(P(j).Cv) * (V(i, :) - P(j).u)') - m);
                    normalization = normalization + w(j, i);
                end
                p(:, i) = w(:, i) / normalization;
            else
                p(:, t) = zeros(J,1);
            end
        end


        % With these conditional probabilities update the values of the parameters    
        for j = 1 : J
            % update the cluster probabilities
            P(j).a = sum(p(j, :)) / Tobs;
        end
        for j = 1 : J            
            % compute the new covariance
            Cv_new = zeros(size(P(j).Cv));
            for i = 1 : T
                if obs_id(i) == 1
                    Cv_new = Cv_new + p(j, i) * (V(i, :) - P(j).u)' * (V(i, :) - P(j).u);
                end
            end   
            Cv_new = Cv_new / sum(p(j, :));
        end
        for j = 1 : J  
            % update the means of neurons clusters
            u_new = zeros(size(P(j).u));        
            for i = 1 : T
                if obs_id(i) == 1
                    u_new = u_new + p(j, i) * V(i, :); 
                end
            end 
            u_new = u_new / sum(p(j, :));
            P(j).u = u_new;
        end
        for j = 1 : J
            % update the covariance
            P(j).Cv = Cv_new;
        end    
     end
    
    % assign cluster ids
    cl_id = zeros(size(obs_id));
    for k = 1 : T
        if obs_id(k) == 1
            [~, I] = max(p(:,k));
            for j = 1 : J
                if (I == j) ~= 0 
                    cl_id(k) = j;
                end
            end
        else
            cl_id(k) = -1;
        end
    end     
end