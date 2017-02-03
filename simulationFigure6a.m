%This Matlab script can be used to generate Figure 6a, in the article:
%
%Emil Bjornson, Elisabeth de Carvalho, Jesper H. Sorensen, Erik G. Larsson,
%Petar Popovski, "A Random Access Protocol for Pilot Allocation in Crowded
%Massive MIMO Systems," IEEE Transactions on Wireless Communications,
%To appear.
%
%Download article: http://arxiv.org/pdf/1604.04248
%
%This is version 1.0 (Last edited: 2017-02-03)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.


%Initialization
close all;
clear;



%Set the number of Monte-Carlo realizations
nbrOfRealizations = 10000;

%Range of BS antennas
Mvalues = [100 300 500];

%Extract maximum number of BS antennas
Mmax = max(Mvalues);


%Define simulation scenario

%Set the pathloss of the two users
betaDifferencedB = -6:0.5:6;

beta1 = ones(size(betaDifferencedB)); %Channel variances of user 1
beta2 = beta1.*10.^(betaDifferencedB/10); %Channel variances of user 2

rho = 1;    %Transmit power of both users
q = 1;      %Transmit power of BS
sigma2 = 1; %Noise variance
taup = 30;  %Length of pilot sequence



%Generate random channel realizations of the two users
h1 = sqrt(1/2)*(randn(Mmax,nbrOfRealizations)+1i*randn(Mmax,nbrOfRealizations));
h2 = sqrt(1/2)*(randn(Mmax,nbrOfRealizations)+1i*randn(Mmax,nbrOfRealizations));

%Generate noise realizations at the BS
n = sqrt(sigma2/2)*(randn(Mmax,nbrOfRealizations)+1i*randn(Mmax,nbrOfRealizations));

%Generate noise realizations at the two users
eta1 = sqrt(sigma2/2)*(randn(1,nbrOfRealizations)+1i*randn(1,nbrOfRealizations));
eta2 = sqrt(sigma2/2)*(randn(1,nbrOfRealizations)+1i*randn(1,nbrOfRealizations));


%Define matrices for storing retransmission probabilities
retransmit1_montecarlo = zeros(length(beta1),length(Mvalues));
retransmit1_analytic = zeros(length(beta1),length(Mvalues));
retransmit2_montecarlo = zeros(length(beta1),length(Mvalues));
retransmit2_analytic = zeros(length(beta1),length(Mvalues));


%Go through all beta values
for betaInd = 1:length(beta1)
    
    %Display simulation progress
    disp(['beta: ' num2str(betaInd) ' out of ' num2str(length(beta1))]);
    
    %Go through all number of antennas
    for mInd = 1:length(Mvalues)
        
        %Compute the resulting value of alpha_t in Eq. (15)
        alphat = rho*taup*(beta1(betaInd)+beta2(betaInd));
        
        %Compute the received signal in Eq. (6)
        yt = sqrt(rho*taup*beta1(betaInd))*h1(1:Mvalues(mInd),:) + sqrt(rho*taup*beta2(betaInd))*h2(1:Mvalues(mInd),:) + n(1:Mvalues(mInd),:);
        
        %Compute the precoding vector used by the BS
        v = sqrt(q)*yt./repmat(sqrt(sum(abs(yt).^2,1)),[Mvalues(mInd) 1]);
        
        %Compute the received DL signal at the two users in Eq. (13)
        z1 = sqrt(taup*beta1(betaInd))*sum(conj(h1(1:Mvalues(mInd),:)).*v,1) + eta1;
        z2 = sqrt(taup*beta2(betaInd))*sum(conj(h2(1:Mvalues(mInd),:)).*v,1) + eta2;
        
        
        %Compute estimate of alpha_t at user 1 using Approx2 in Eq. (36)
        alphaEst1_approx2 = exp(gammaln(Mvalues(mInd)+1/2)-gammaln(Mvalues(mInd)))^2*q*rho*beta1(betaInd)^2*taup^2./real(z1).^2-sigma2;
        alphaEst1_approx2(alphaEst1_approx2<rho*beta1(betaInd)*taup) = rho*beta1(betaInd)*taup;
        
        %Compute estimate of alpha_t at user 2 using Approx2 in Eq. (36)
        alphaEst2_approx2 = exp(gammaln(Mvalues(mInd)+1/2)-gammaln(Mvalues(mInd)))^2*q*rho*beta2(betaInd)^2*taup^2./real(z2).^2-sigma2;
        alphaEst2_approx2(alphaEst2_approx2<rho*beta2(betaInd)*taup) = rho*beta2(betaInd)*taup;
        
        
        %Apply the retransmission decision rule for each user
        retransmit1Realizations = rho*beta1(betaInd)*taup>alphaEst1_approx2/2;
        retransmit2Realizations = rho*beta2(betaInd)*taup>alphaEst2_approx2/2;
        
        %Compute the probability of retransmission by Monte-Carlo methods
        retransmit1_montecarlo(betaInd,mInd) = mean(retransmit1Realizations);
        retransmit2_montecarlo(betaInd,mInd) = mean(retransmit2Realizations);
       
        
        %Compute zeta_k from Eq. (38) for each of the user
        zeta1 = exp(gammaln(Mvalues(mInd)+1/2)-gammaln(Mvalues(mInd)))^2*q*rho*beta1(betaInd)^2*taup^2./(sigma2 + 2*(rho*beta1(betaInd)*taup));
        zeta2 = exp(gammaln(Mvalues(mInd)+1/2)-gammaln(Mvalues(mInd)))^2*q*rho*beta2(betaInd)^2*taup^2./(sigma2 + 2*(rho*beta2(betaInd)*taup));
        
        %Compute the probability of retransmission by using Theorem 2
        retransmit1_analytic(betaInd,mInd) = 1 - computeZreCDF(alphat,sqrt(zeta1),rho,q,beta1(betaInd),taup,sigma2,Mvalues(mInd)) + computeZreCDF(alphat,-sqrt(zeta1),rho,q,beta1(betaInd),taup,sigma2,Mvalues(mInd));
        retransmit2_analytic(betaInd,mInd) = 1 - computeZreCDF(alphat,sqrt(zeta2),rho,q,beta2(betaInd),taup,sigma2,Mvalues(mInd)) + computeZreCDF(alphat,-sqrt(zeta2),rho,q,beta2(betaInd),taup,sigma2,Mvalues(mInd));
        
    end
    
end



%%Plot simulation results

figure(1);
hold on; box on;

plot(betaDifferencedB,retransmit1_montecarlo(:,1),'k-.','LineWidth',1);
plot(betaDifferencedB,retransmit1_montecarlo(:,3),'b-','LineWidth',1);

plot(betaDifferencedB,retransmit2_montecarlo(:,1),'k-.','LineWidth',1);
plot(betaDifferencedB,retransmit2_montecarlo(:,3),'b-','LineWidth',1);

plot(betaDifferencedB,retransmit1_analytic(:,1),'ko','LineWidth',1);
plot(betaDifferencedB,retransmit1_analytic(:,3),'bo','LineWidth',1);

plot(betaDifferencedB,retransmit2_analytic(:,1),'kd','LineWidth',1);
plot(betaDifferencedB,retransmit2_analytic(:,3),'bd','LineWidth',1);

xlabel('Difference in SNR');
ylabel('Probability of Retransmission');
legend('M=100','M=500','Location','SouthEast');
ylim([0 1]);
