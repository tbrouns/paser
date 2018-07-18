%
%online sorting
%
%thresholdMethod:
%1 -> use raw data. stdEstimate is threshold in this case.
%2-> use covariance. stdEstimate is chi2 value in this case. covariance
%matrix Cinv needs to be given. Cinv is the inverted covariance matrix.
%
%
%urut/nov04; major revision urut/nov05
function [assigned, nrAssigned, baseSpikes, baseSpikesID] = psr_sst_ost_sortSpikesOnline( spikes, stdEstimate, sortTill, thresholdMethod, Cinv, transformedSpikes )
if nargin<=4
    thresholdMethod=1;
    Cinv=[];
    transformedSpikes=[];
end

defineSortingConstants;

%ID of noise cluster
noiseCluster=CLUSTERID_NOISE_CLUSTER;

%how many datapoints per spike
nrDatapoints=size(spikes,2);

baseSpikes=spikes(1,:);		% mean spikes of each cluster (here called baseSpikes)
baseSpikesCounter=1;		% number of spikes associated to each cluster
baseSpikesID=[1 1];		% mapping of the ID (unique) of each cluster to index in base Spikes. the ID is used to assign spikes in the variable assigned.
IDs=1;				% next free ID (used for creation of new clusters)
nrBaseSpikes=0;			% how many base spikes are known.

%-- threshold
weights=ones(1,nrDatapoints);

%TODO testing
if nrDatapoints==256
    weights(1:20)=0;
    weights(201:256)=0;
end

if thresholdMethod==1
    thres = stdEstimate^2 * length(find(weights));
    disp(['Threshold method: approximation. Threshold is: ' num2str(thres)]);
    %correction factor for bursts (see paper for calculation)
    %  thres=thres*1.2;
    
    thres=thres*1.2;
    thresMerge=thres;
    
    %thresMerge=thres/2;
    
else
    thres=stdEstimate;
    thresMerge=4;  %how many std appart for projection test
    disp(['Threshold method: exact Chi2. Threshold is: ' num2str(thres)]);
end

weights=weights';

%-- ID of each spike
assigned=[];

if sortTill > size(spikes, 1)
    sortTill=size(spikes,1);
end

disp(['Number spikes to be sorted: ' num2str(sortTill)]);
for i=2:sortTill
    if mod(i,100)==0
        disp(['sorted ' num2str(i)]);
    end
    
    if mod(i,1000)==0
        disp(['sorted ' num2str(i) ' of ' num2str(size(spikes,1))]);
        
        %prune of noise clusters recursively
        lowerLimit=2;
        indsRemove = find(baseSpikesCounter<lowerLimit);
        
        while ( length(indsRemove)>0 )
            indToRemove = indsRemove(1);
            
            %assign spikes assigned to this basespike to noise
            IDtoRemove = baseSpikesID( find( baseSpikesID(:,2)==indToRemove), 1);
            assigned(find(assigned==IDtoRemove)) = noiseCluster;
            
            %remove this base spike
            indToKeep = setdiff(1:size(baseSpikes,1),indToRemove);
            baseSpikes = baseSpikes( indToKeep,:);
            baseSpikesCounter = baseSpikesCounter( indToKeep );
            nrBaseSpikes=nrBaseSpikes-1;
            
            %remove ID of base spike in mapping table
            baseSpikesID=baseSpikesID(find(baseSpikesID(:,1)~=IDtoRemove),:);
            
            %correct mapping table for shift caused by this
            indsToShift = find( baseSpikesID(:,2) > indToRemove );
            if length(indsToShift)>0
                baseSpikesID(indsToShift,2) = baseSpikesID(indsToShift,2)-1;
            end
            
            indsRemove = find(baseSpikesCounter<lowerLimit);
        end
        disp('pruning finished');
    end
    
    testSpike=spikes(i,:);
    
    D=calcDistMacro( thresholdMethod, baseSpikes, testSpike, weights,Cinv);
    
    
    if min(D)>thres
        %add as new base spike
        nrBaseSpikes=nrBaseSpikes+1;
        baseSpikes(nrBaseSpikes,:)=testSpike;
        baseSpikesCounter(nrBaseSpikes)=1;
        IDs=IDs+1;
        baseSpikesID(nrBaseSpikes,1:2)=[IDs nrBaseSpikes];
        
        assigned(i)=IDs;
        
        disp(['new cluster created ' num2str(IDs) ]);
    else
        Dsorted=sortrows(D);
        
        for j=1:length(Dsorted)
            ind=find(Dsorted(j)==D,1);
            
            if D(ind)>thres
                assigned(i)=noiseCluster;
                break;
            end
            
            %[min(D) max(D) length(D) length(ind) ind size(baseSpikes,1)]
            
            %test envelope
            %             if length(find(abs(baseSpikes(ind,:)-testSpike) >envSize))>0
            %                  continue;
            %             end
            
            %env at peak
            
            %             if abs ( baseSpikes(ind,95)-testSpike(95) ) > envSize
            %                 continue;
            %             end
            
            %assign the spike to this cluster and increase counter
            IDassigned=baseSpikesID(find(baseSpikesID(:,2)==ind)  ,1);
            assigned(i)=IDassigned;
            baseSpikesCounter(ind)=baseSpikesCounter(ind)+1;
            
            %disp(['assigned spike to existing cluster ' num2str(IDassigned)]);
            
            %--update base spike as mean of last 100 spikes assigned to this cluster
            indsThisCluster = find(assigned==IDassigned) ;
            runningAverageLength=100;
            if length(indsThisCluster)>runningAverageLength
                indsThisCluster = indsThisCluster(end-runningAverageLength:end);
            end
            
            %after a cluster has X spikes assigned,don't update waveform
            %anymore (has converged)
            %if length(indsThisCluster)<10
            
            if thresholdMethod==2
                baseSpikes(ind,:)= mean ( spikes ( indsThisCluster,: ));
            else
                if length(indsThisCluster)<50
                    baseSpikes(ind,:)= mean ( spikes ( indsThisCluster,: ));
                end
            end
            %end
            %--
            
            %-- merging
            %if following break is uncommented: no merging (if you want to do this manually, it will overcluster but will otherwise be ok and it might be possible to differentiate difficult cases manually)
            %break;
            
            %merging only necessary for std threshold method
            %if thresholdMethod==2
            %    break;
            %end
            
            
            
            %ID1 is same as IDassigned
            ID1 = baseSpikesID(find(baseSpikesID(:,2)==ind),1);
            newBaseSpike=baseSpikes(ind,:);
            
            
            
            %after merge --> are we now closer than thres to something?
            
            %DD=calcDistMacro( thresholdMethod, baseSpikes, newBaseSpike, weights,Cinv);
            
            %[Dmerge,Tmerge] = calculateDistanceT2(baseSpikes, newBaseSpike, ind, Cinv, baseSpikesCounter);
            %DD = Tmerge-Dmerge;
            
            if thresholdMethod==1
                DD=calcDistMacro( thresholdMethod, baseSpikes, newBaseSpike, weights,Cinv);
            else
                %DD = calculateDistanceNorm(transformedSpikes, baseSpikes, baseSpikesID, ID1, assigned);
                DD = calculateDistanceNorm2( baseSpikes, newBaseSpike);
                %DD=calculateDistanceNorm(baseSpikes, newBaseSpike);
            end
            
            
            while min(DD(find(DD>0))) < (thresMerge)
                disp(['merg candidate ' num2str(min(DD(find(DD>0))))]);
                %one will be zero,exclude
                DDD=DD(find(DD>0));
                %[min(DDD) size(baseSpikes)]
                if min(DDD)< (thresMerge)
                    ind2=find(min(DDD)==DD,1);
                    
                    if length(ind2)>0
                        
                        IDtoDelete = baseSpikesID(find(baseSpikesID(:,2)==ind2),1);
                        
                        indsToKeep=setdiff(1:size(baseSpikes,1),ind2);
                        
                        %update before remove
                        %[ind ind2]
                        
                        c1=baseSpikesCounter(ind);
                        c2=baseSpikesCounter(ind2);
                        totNr = c1 +  c2;
                        
                        if c1<2 || c2<2
                            disp(['do not merge clusters because too few spikes assigned ID1/ID2=' num2str(ID1) '/' num2str(IDtoDelete)]);
                            DD(ind2) = -999;
                            continue;
                        end
                        
                        disp(['merging n/ID=' num2str(c1) '/' num2str(ID1) ' w. n/ID=' num2str(c2) '/' num2str(IDtoDelete)]);
                        
                        %           	 figure(897)
                        %              subplot(2,2,1)
                        %              plot(1:size(baseSpikes,2),spikes(find(assigned==ID1),:), 'r', 1:size(baseSpikes,2),spikes(find(assigned==IDtoDelete),:), 'b');
                        %              subplot(2,2,2)
                        %              plot(1:size(baseSpikes,2), baseSpikes(ind2,:), 'b',1:size(baseSpikes,2), baseSpikes(ind,:),'r');
                        %              subplot(2,2,3)
                        %              plot(1:size(baseSpikes,2),mean(transformedSpikes(find(assigned==ID1),:)), 'r', 1:size(baseSpikes,2),mean(transformedSpikes(find(assigned==IDtoDelete),:)), 'b');
                        %              subplot(2,2,4)
                        %              [pc,score,latent,tsquare] = princomp( [spikes(find(assigned==ID1),:); spikes(find(assigned==IDtoDelete),:); baseSpikes(ind2,:); baseSpikes(ind,:)] );
                        %              plot(score(1:c1,1), score(1:c1,2),'.r',score(c1+1:end-2,1), score(c1+1:end-2,2),'.b', score(end-1,1), score(end-1,2),'db', score(end,1), score(end,2),'dr');
                        %
                        %              norm( mean(transformedSpikes(find(assigned==ID1),:)) - mean(transformedSpikes(find(assigned==IDtoDelete),:)) )
                        %              norm( mean(spikes(find(assigned==ID1),:)) - mean(spikes(find(assigned==IDtoDelete),:)) )
                        %
                        %              norm( baseSpikes(ind2,:) - baseSpikes(ind,:) )
                        %             title(['will merge ID=' num2str(ind2) ' n=' num2str(c2) ' w. ID=' num2str(ind) ' n=' num2str(c1) ' d=' num2str(min(DDD))]);
                        %
                        
                        baseSpikes(ind,:)=(baseSpikes(ind,:)*c1+baseSpikes(ind2,:)*c2)/totNr;
                        baseSpikesCounter(ind)=baseSpikesCounter(ind)+baseSpikesCounter(ind2);
                        
                        baseSpikesOld=baseSpikes;
                        baseSpikes=[];
                        baseSpikes= baseSpikesOld(indsToKeep,:);
                        baseSpikesCounter= baseSpikesCounter(indsToKeep);
                        nrBaseSpikes=nrBaseSpikes-length(ind2);
                        
                        assigned(find(assigned==IDtoDelete))=ID1;
                        %assigned(find(assigned==IDtoDelete)) = noiseCluster;
                        
                        
                        baseSpikesID=baseSpikesID(find(baseSpikesID(:,1)~=IDtoDelete),:);
                        
                        
                        %shift inds
                        if size(baseSpikes,1)+1>ind2
                            indsToShift = ind2+1:size(baseSpikes,1)+1;
                            indCorrected=false;
                            
                            for k=1:length(indsToShift)
                                ii = find( baseSpikesID(:,2)==indsToShift(k) );
                                
                                baseSpikesID(ii,2) = baseSpikesID(ii,2)-1;
                                
                                if ind==indsToShift(k) && indCorrected==false
                                    ind=ind-1;
                                    indCorrected=true;
                                end
                            end
                        end
                    end
                else
                    break;
                end
                %DD=calculateDistance(baseSpikes, baseSpikes(ind,:), weights);
                %DD=calcDistMacro( thresholdMethod, baseSpikes, baseSpikes(ind,:), weights,Cinv);
                if thresholdMethod==1
                    DD=calcDistMacro( thresholdMethod, baseSpikes, baseSpikes(ind,:), weights,Cinv);
                else
                    %DD=calculateDistanceNorm(baseSpikes, baseSpikes(ind,:));
                    %DD = calculateDistanceNorm(transformedSpikes, baseSpikes, baseSpikesID, ind, assigned);
                    DD = calculateDistanceNorm2( baseSpikes, baseSpikes(ind,:));
                end
                
                %break;
                %[Dmerge,Tmerge] = calculateDistanceT2(baseSpikes, baseSpikes(ind,:), ind, Cinv, baseSpikesCounter);
                %DD = Tmerge-Dmerge;
            end
            break;
        end
        
    end
end


%rank
nrAssigned=zeros(size(baseSpikesID));
nrAssigned(:,1)=baseSpikesID(:,1);
for i=1:size(baseSpikesID,1)
    indOfSpikes=find(assigned==baseSpikesID(i,1));
    nrAssigned(i,2) = length(indOfSpikes);
end

nrAssigned=sortrows(nrAssigned,2);

%----internal functions
