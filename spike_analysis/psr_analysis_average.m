close all

params = psr_analysis_parameters();
t      = -params.t_pre : params.t_bin : params.t_post;
t_off  = t(1:end-1) + 0.5 * params.t_bin;

nstim     = size(GPSTH_all,1);
nclusters = size(GPSTH_all,2);
% TOTAL     = [];
N         = [];
for iStim = 2:nstim
%     total = [];
    for iCluster = 1:nclusters
        n = mean(GPSTH_all{iStim,iCluster},2);
%         if (isempty(total))
%             total = zeros(size(N));
%         end
%         total = total + N;
        N = [N,n];
    end
%     total = total / nclusters;
%     TOTAL = [TOTAL,total]; %#ok
    
end

% TOTAL = mean(TOTAL,2);

mu = mean(N,2);
sd = std(N,[],2);
errorbar(mu,sd);
ylim([floor(min(mu)) ceil(max(mu))]);

% figure;
% plot(t_off,TOTAL,'LineWidth',1.5);

PSTH = mean(PSTH_all,2);

figure;
plot(t_off,PSTH,'LineWidth',1.5);
ylim([floor(min(PSTH)) ceil(max(PSTH))]);