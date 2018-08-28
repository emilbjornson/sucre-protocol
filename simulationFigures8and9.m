%This Matlab script can be used to generate Figures 8 and 9, in the article:
%
%Emil Bjornson, Elisabeth de Carvalho, Jesper H. Sorensen, Erik G. Larsson,
%Petar Popovski, "A Random Access Protocol for Pilot Allocation in Crowded
%Massive MIMO Systems," IEEE Transactions on Wireless Communications,
%To appear.
%
%Download article: http://arxiv.org/pdf/1604.04248
%
%This is version 1.01 (Last edited: 2018-08-28)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.


%Initialization
close all;
clear;



%Toggle between:
%TRUE = Adjacent cells cause interference by sending random data (Gives (a) in Figure 8 and 9)
%FALSE = Adjacent cells are silent (Gives (b) in Figure 8 and 9)
includeIntercellInterference = false;

%Set number of channel models that are considered. The following are added
%one at a time as the variable is increased
%1) Uncorrelated Rayleigh fading
%2) Correlated Rayleigh fading
%3) Line-of-sight propagation
%4) Line-of-sight propagation with random power backoff
numberOfChannelModels = 4;


%Set the number of Monte-Carlo realizations
nbrOfRealizations = 10000;

%Range of BS antennas
Mvalues = [1 5:5:20 25:25:200];

%Extract maximum number of BS antennas
Mmax = max(Mvalues);

%Maximum number of colliding users in the simulation
maxUsers = 11;

%Range of delta values considered in Figure 9
biasDelta = [0 -2:0.5:2];




%%Define simulation scenario

rhoOriginal = 1; %Transmit power of UEs in the cell
q = 1;           %Transmit power of the BS
sigma2 = 1;      %Noise variance
taup = 30;       %Length of RA pilot sequence


%Standard deviation of shadow fading in line-of-sight and non-line-of-sight
shadowFadingStddB_nlos = 10;
shadowFadingStddB_los = 4;

%Generate shadow fading realizations
shadowFadingRealizations = randn(1,nbrOfRealizations,maxUsers);



%Set the inter-cell interference parameters
if includeIntercellInterference == true
    
    qNeighbor = taup*q; %Transmit power of neighbouring BSs
    rhoIntercell = rhoOriginal; %Set transmit power of users in adjacent cells
    
else
    
    qNeighbor = 0; %No transmit power of neighbouring BSs
    rhoIntercell = 0; %No transmit power of users in adjacent cells
    
end


%Set cell radius (in meter)
cellRadius = 250;

%Generate user locations within the cell
userLocations = generatePointsHexagon([1,nbrOfRealizations,maxUsers],cellRadius,0.1*cellRadius);

%Extract user distances from the BS
userDistance = abs(userLocations);

%Extract angles between users and the BS
userAngles = angle(userLocations);



%Create matrix with individual user transmit powers
rhoUsers = rhoOriginal*ones(1,nbrOfRealizations,maxUsers);

%Random backoff in transmit power in the line-of-sight case (dB scale)
rhoModValues = -30*rand(1,nbrOfRealizations,maxUsers);



%Compute locations of all neighboring BSs
neighboringBSs = 2*cellRadius*sqrt(3)/2 * exp(1i*(0:5)*pi/3);

%Generate shadow fading realizations for downlink inter-cell interference
shadowFadingRealizationsIntercellDownlink = randn([length(neighboringBSs),nbrOfRealizations,maxUsers]);


%Go through all users and make sure that they are always served by the BS
%in the own hexagonal cell (even when shadow fading is taken into account)
for k = 1:maxUsers
    
    notReady = 1;
    
    while sum(notReady)>0
        
        notReady = zeros(1,nbrOfRealizations);
        
        for j = 1:length(neighboringBSs)
            %Check which of the users that are served by the right BS
            notReady = notReady + (shadowFadingStddB_nlos*shadowFadingRealizations(1,:,k) - 38*log10(userDistance(1,:,k)) < shadowFadingStddB_nlos*shadowFadingRealizationsIntercellDownlink(j,:,k) - 38*log10(abs(userLocations(1,:,k)-neighboringBSs(j))) );
        end
        
        %Make new random shadow fading realizations for those users that
        %have a better channel to the neighboring BS
        shadowFadingRealizations(1,notReady>0,k) = randn(1,sum(notReady>0));
        shadowFadingRealizationsIntercellDownlink(:,notReady>0,k) = randn(length(neighboringBSs),sum(notReady>0));
        
    end
    
end


%Set number of active UEs in the neighboring cells
K_neighboringcells = 10;

%Generate user locations in neighboring cells
userLocationsNeighboring = generatePointsHexagon([K_neighboringcells,nbrOfRealizations,length(neighboringBSs)],cellRadius,0.1*cellRadius);

%Generate shadow fading realizations of users in neighboring cells, both
%within that cell and to the cell under study
shadowFadingRealizationsWithinOwnCellUplink = randn([K_neighboringcells,nbrOfRealizations,length(neighboringBSs)]);
shadowFadingRealizationsIntercellUplink = randn([K_neighboringcells,nbrOfRealizations,length(neighboringBSs)]);


%Go through all users in neighboring cells
for k = 1:K_neighboringcells
    
    notReady = 1;
    
    while sum(notReady)>0
        
        notReady = zeros(length(neighboringBSs),nbrOfRealizations);
        
        for j = 1:length(neighboringBSs)
            
            %Check which of the users that are served by the right BS
            notReady(j,:) = (shadowFadingStddB_nlos*shadowFadingRealizationsWithinOwnCellUplink(k,:,j) - 38*log10(abs(userLocationsNeighboring(k,:,j))) < shadowFadingStddB_nlos*shadowFadingRealizationsIntercellUplink(k,:,j) - 38*log10(abs(userLocationsNeighboring(k,:,j)-neighboringBSs(j))) );
            
            %Make new random shadow fading realizations for those users in
            %the neighborin cell that have a better channel to the center BS
            shadowFadingRealizationsWithinOwnCellUplink(k,notReady(j,:)>0,j) = randn(1,sum(notReady(j,:)>0));
            shadowFadingRealizationsIntercellUplink(k,notReady(j,:)>0,j) = randn(1,sum(notReady(j,:)>0));
            
        end
        
    end
    
end


%Compute the total inter-cell interference in the uplink and downlink in
%each channel realization
interCellVarianceDownlink = zeros(size(userLocations));
interCellVarianceUplink = zeros(1,nbrOfRealizations);

for j = 1:length(neighboringBSs)
    
    %Note: 27 dBm represents the transmit power and -98.65 dBm represents
    %the noise variance
    interCellVarianceDownlink = interCellVarianceDownlink + qNeighbor*10.^( (27 + 98.65 - 34.53 - 38*log10(abs(userLocations-neighboringBSs(j))) + shadowFadingStddB_nlos*shadowFadingRealizationsIntercellDownlink(j,:,:)  )/10   );
    interCellVarianceUplink = interCellVarianceUplink + sum(rhoIntercell*10.^( (27 + 98.65 - 34.53 - 38*log10(abs(userLocationsNeighboring(:,:,j) + neighboringBSs(j))) + shadowFadingStddB_nlos*shadowFadingRealizationsIntercellUplink(:,:,j)  )/10   ),1) ;
    
end

%Compute the average uplink inter-cell interference, which is used in the
%bias terms
interCellBias = mean(interCellVarianceUplink);


%Compute average signal gain for line-of-sight and non-line-of-sight
%propagation (27 dBm represents the transmit power and -98.65 dBm
%represents the noise variance)
betas_nlos = 10.^( (27+ 98.65 - 34.53 - 38*log10(userDistance) + shadowFadingStddB_nlos*shadowFadingRealizations  )/10   );
betas_los = 10.^( (27+ 98.65 - 30.18 - 26*log10(userDistance) + shadowFadingStddB_los*shadowFadingRealizations  )/10   );


%Generate uncorrelated Rayleigh fading channel realizations
h_iid = (randn(Mmax,nbrOfRealizations,maxUsers)+1i*randn(Mmax,nbrOfRealizations,maxUsers));
h_iid = repmat(sqrt(betas_nlos/2),[Mmax 1 1 ]) .* h_iid;


%Define distance between antennas in the ULA, measured in wavelengths
interAntennaDistance = 0.5;

%Generate line-of-sight channel realizations using Eq. (41)
h_los = exp( repmat((0:Mmax-1)',[1 nbrOfRealizations maxUsers]) .* repmat(-2i*pi*interAntennaDistance*sin(userAngles),[Mmax 1 1]) );
h_los = repmat(sqrt(betas_los),[Mmax 1 1 ]) .* h_los;



%Select correlation factor in the exponential correlation model
correlationFactor = 0.7;


%Vector to store correlated Rayleigh fading
h_corr = zeros(Mmax,nbrOfRealizations,maxUsers);


%Go through all channel realizations
for j = 1:nbrOfRealizations
    
    %Display simulation progress at every 100 channel realizations
    if mod(j,100) == 0
        disp(['Prepare channel: ' num2str(j) ' out of ' num2str(nbrOfRealizations)]);
    end
    
    %Go through all users
    for userInd = 1:maxUsers
        
        %Make the Rayleigh fading channel realizations correlated according to Eq. (40)
        R = toeplitz((correlationFactor*exp(1i*sin(userAngles(1,j,userInd)))).^(0:Mmax-1));
        h_corr(:,j,userInd) = sqrtm(R)*h_iid(:,j,userInd);
        
    end
    
end



%Generate noise realizations at the BS in the center cell
n = sqrt(sigma2/2)*(randn(Mmax,nbrOfRealizations)+1i*randn(Mmax,nbrOfRealizations));

%Generate inter-cell interference realizations at the BS
w = sqrt(repmat(interCellVarianceUplink/2,[Mmax 1])).*(randn(Mmax,nbrOfRealizations)+1i*randn(Mmax,nbrOfRealizations));

%Generate noise realizations at the users
eta = sqrt(sigma2/2)*(randn(1,nbrOfRealizations,maxUsers)+1i*randn(1,nbrOfRealizations,maxUsers));

%Generate inter-cell interference realizations at the users
omega = sqrt(interCellVarianceDownlink/2).*(randn(1,nbrOfRealizations,maxUsers)+1i*randn(1,nbrOfRealizations,maxUsers));


%Define matrix to solve probabilities to resolve collisions, false
%negatives, and false positives
resolved = zeros(length(Mvalues),numberOfChannelModels,maxUsers,length(biasDelta));
falseNegative = zeros(length(Mvalues),numberOfChannelModels,maxUsers,length(biasDelta));
falsePositive = zeros(length(Mvalues),numberOfChannelModels,maxUsers,length(biasDelta));


%Go through all number of colliding users
for userInd = 1:maxUsers
    
    %Display simulation progress
    disp(['beta: ' num2str(userInd) ' out of ' num2str(maxUsers)]);
    
    %Go through all number of antennas
    for mInd = 1:length(Mvalues)
        
        %Go through all channel models and extract the channel realizations
        for cInd = 1:numberOfChannelModels
            
            
            if cInd == 1 %Uncorrelated Rayleigh fading
                
                h = h_iid(1:Mvalues(mInd),:,1:userInd);
                betas = betas_nlos(1,:,1:userInd);
                rho = rhoUsers;
                
            elseif cInd == 2 %Correlated Rayleigh fading
                
                h = h_corr(1:Mvalues(mInd),:,1:userInd);
                betas = betas_nlos(1,:,1:userInd);
                rho = rhoUsers;
                
            elseif cInd == 3 %Line-of-sight with full power
                
                h = h_los(1:Mvalues(mInd),:,1:userInd);
                betas = betas_los(1,:,1:userInd);
                rho = rhoUsers;
                
            elseif cInd == 4 %Line-of-sight with random power backoff
                
                h = h_los(1:Mvalues(mInd),:,1:userInd);
                betas = betas_los(1,:,1:userInd);
                rho = rhoUsers.*10.^(rhoModValues/10);
                
            end
            
            
            %Matrix to store which UEs that retransmit
            retransmit = zeros(userInd,nbrOfRealizations,1,length(biasDelta));
            
            
            
            %Compute the received signal in Eq. (6)
            yt = sqrt(taup) * sum(repmat(sqrt(rho(1,:,1:userInd)),[Mvalues(mInd) 1 1]).*h,3) + w(1:Mvalues(mInd),:) + n(1:Mvalues(mInd),:);
            
            %Compute the precoding vector used by the BS
            v = sqrt(q)*yt./repmat(sqrt(sum(abs(yt).^2,1)),[Mvalues(mInd) 1]);
            
            
            for k = 1:userInd
                
                %Compute the received DL signal at user k in Eq. (13)
                z = sqrt(taup)*sum(conj(h(1:Mvalues(mInd),:,k)).*v,1) + omega(1,:,k) + eta(1,:,k);
                
                %Compute estimate of alpha_t at user k using Approx2 in Eq. (36)
                alphaEst_approx2 = exp(gammaln(Mvalues(mInd)+1/2)-gammaln(Mvalues(mInd)))^2*q*taup^2*rho(1,:,k).*(betas(1,:,k)./real(z)).^2-sigma2;
                indicesTooLow = alphaEst_approx2<rho(1,:,k).*betas(1,:,k)*taup;
                alphaEst_approx2(indicesTooLow) = rho(1,indicesTooLow,k).*betas(1,indicesTooLow,k)*taup;
                
                %Go through all different bias terms
                for b = 1:length(biasDelta)
                    
                    %Check if the UEs will retransmit
                    retransmit(k,:,b) = rho(1,:,k).*betas(1,:,k)*taup>(alphaEst_approx2-interCellBias)/2+rho(1,:,k).*betas(1,:,k)*taup*biasDelta(b)/sqrt(Mvalues(mInd));
                    
                end
                
            end
            
            
            %Compute the probability that exactly one of the users retransmit
            resolved(mInd,cInd,userInd,:) = sum(sum(retransmit,1)==1,2)/nbrOfRealizations;
            
            falseNegative(mInd,cInd,userInd,:) = sum(sum(retransmit,1)==0,2)/nbrOfRealizations;
            
            falsePositive(mInd,cInd,userInd,:) = sum(sum(retransmit,1)>1,2)/nbrOfRealizations;

        end
        
    end
    
end



%%User the probabilities of resolving collisions for different number of
%%colliding UEs and multiply them with the probability that each of these
%%cases occur.


%Number of inactive UEs in the cell
K0 = 5000;

%Probability that a UE wants to become active in a block
pA = 0.005;

%Probability of sending an RA pilot if a UE wants to become active
pP = 0:0.01:1;

%Number of RA pilot signals
taup = 10;


%Range of number of colliding users per RA pilot
userValues = 0:maxUsers;

%Define vector to store probabilities of having different number of
%colliding UEs
probabilities = zeros(length(userValues),length(pP));

for kind = 1:length(userValues)
    
    %Compute probability according to the binomial distribution in Eq. (1),
    %by using logarithms of each term to gain numerical stability
    probabilities(kind,:) =  exp(  gammaln(K0+1) - gammaln(userValues(kind)+1) - gammaln(K0-userValues(kind)+1) + (userValues(kind)).*log(pA*pP/taup) + (K0-userValues(kind)).*log((1-pA*pP/taup)) );
    
end


%Use the probabilties of having different number of colliding UEs to obtain
%the probability to resolve collisions, having false negatives and false
%positives.
acceptedUsers = zeros(length(Mvalues),numberOfChannelModels,1,length(biasDelta));
noUsers = zeros(length(Mvalues),numberOfChannelModels,1,length(biasDelta));
remainingCollisions = zeros(length(Mvalues),numberOfChannelModels,1,length(biasDelta));

for k = 1:maxUsers
    acceptedUsers = acceptedUsers +  probabilities(k+1,end)*resolved(:,:,k,:);
    noUsers = noUsers +  probabilities(k+1,end)*falseNegative(:,:,k,:);
    remainingCollisions = remainingCollisions +  probabilities(k+1,end)*falsePositive(:,:,k,:);
end

%Renormalize the probablities to remove cases where no user is
%transmitting, that is cases where there is no interest to transmit
acceptedUsers = acceptedUsers / (1 - probabilities(1,end));
noUsers = noUsers / (1 - probabilities(1,end));
remainingCollisions = remainingCollisions / (1 - probabilities(1,end));


%Compute the probability of resolving collisions in the baseline scheme
%where only pilots that are used by one UE are considered resolved. The
%probability pP is optimized numerically by finding the one that maximizes
%the performance of the baseline scheme
baseline = max(probabilities(2,:)) / (1 - probabilities(1,end));



%Extract the results for M=100 and uncorrelated Rayleigh fading, for
%different values of the bias delta term
acceptedUsersM100 = reshape(acceptedUsers(9,1,1,2:end),[length(biasDelta)-1 1]);
falsePositiveM100 = reshape(remainingCollisions(9,1,1,2:end),[length(biasDelta)-1 1]);
falseNegativeM100 = reshape(noUsers(9,1,1,2:end),[length(biasDelta)-1 1]);




%%Plot simulation results

figure(1);
hold on; box on;

plot(Mvalues,acceptedUsers(:,4,1,1),'r--','LineWidth',1);
plot(Mvalues,acceptedUsers(:,1,1,1),'k-.','LineWidth',1);
plot(Mvalues,acceptedUsers(:,2,1,1),'b-','LineWidth',1);
plot(Mvalues,acceptedUsers(:,3,1,1),'r:','LineWidth',1);
plot([min(Mvalues) max(Mvalues)],baseline*ones(1,2),'k:','LineWidth',1);

xlabel('Number of Antennas (M)');
ylabel('Probability to Resolve Collision');
legend('SUCRe: LoS (random power)', 'SUCRe: Uncorr Rayleigh','SUCRe: Corr Rayleigh','SUCRe: LoS','Baseline','Location','SouthEast');
ylim([0 1]);


figure(2);
hold on; box on;

plot(biasDelta(2:end),acceptedUsersM100(:,1),'r-','LineWidth',1);
plot(biasDelta(2:end),falseNegativeM100(:,1),'k--','LineWidth',1);
plot(biasDelta(2:end),falsePositiveM100(:,1),'b-.','LineWidth',1);

xlabel('Bias term (in standard deviations)');
ylabel('Probability');
legend('Resolved collision','False negative (no UEs)','False positive (multiple UE)','Location','best');
ylim([0 1]);
