%This Matlab script can be used to generate Figure 5, in the article:
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
nbrOfRealizations = 2500;

%Range of BS antennas
Mvalues = [10:5:50 60:10:100 120:20:200];

%Extract maximum number of BS antennas
Mmax = max(Mvalues);


%Define simulation scenario

%Set the pathloss of the two users
beta1 = 1; %Channel variance of user 1
beta2 = 1; %Channel variance of user 2

rho = 1;    %Transmit power of both users
q = 1;      %Transmit power of BS
sigma2 = 1; %Noise variance
taup = 30;  %Length of RA pilot sequence

%Compute the resulting value of alpha_t in Eq. (15)
alphat = rho*taup*(beta1+beta2);


%Generate random channel realizations of the two users
h1 = sqrt(beta1/2)*(randn(Mmax,nbrOfRealizations)+1i*randn(Mmax,nbrOfRealizations));
h2 = sqrt(beta2/2)*(randn(Mmax,nbrOfRealizations)+1i*randn(Mmax,nbrOfRealizations));

%Generate noise realizations at the BS
n = sqrt(sigma2/2)*(randn(Mmax,nbrOfRealizations)+1i*randn(Mmax,nbrOfRealizations));

%Generate noise realization for user 1
eta1 = sqrt(sigma2/2)*(randn(1,nbrOfRealizations)+1i*randn(1,nbrOfRealizations));


%Define vectors for storing bias values of the three estimators
bias_approx1 = zeros(length(Mvalues),1);
bias_approx2 = zeros(length(Mvalues),1);
bias_ML = zeros(length(Mvalues),1);

%Define vectors for storing NMSE values of the three estimators
NMSE_approx1 = zeros(length(Mvalues),1);
NMSE_approx2 = zeros(length(Mvalues),1);
NMSE_ML = zeros(length(Mvalues),1);



%Go through all number of antennas
for mInd = 1:length(Mvalues)
    
    %Display simulation progress
    disp(['M: ' num2str(mInd) ' out of ' num2str(length(Mvalues))]);
    
    %Compute the received signal in Eq. (6)
    yt = sqrt(rho*taup)*h1(1:Mvalues(mInd),:) + sqrt(rho*taup)*h2(1:Mvalues(mInd),:) + n(1:Mvalues(mInd),:);
    
    %Compute the precoding vector used by the BS
    v = sqrt(q)*yt./repmat(sqrt(sum(abs(yt).^2,1)),[Mvalues(mInd) 1]);
    
    %Compute the received DL signal at user 1 in Eq. (13)
    z1 = sqrt(taup)*sum(conj(h1(1:Mvalues(mInd),:)).*v,1) + eta1;
    
    
    
    %Compute estimate of alpha_t using Approx1 in Eq. (17)
    alphaEst1_approx1 = Mvalues(mInd)*q*rho*beta1^2*taup^2./real(z1).^2-sigma2;
    alphaEst1_approx1(alphaEst1_approx1<rho*beta1*taup) = rho*beta1*taup;
    
    
    %Compute and store the resulting normalized bias and NMSE with Approx 1
    bias_approx1(mInd) = bias_approx1(mInd) + mean(alphaEst1_approx1-alphat)/alphat;
    NMSE_approx1(mInd) = NMSE_approx1(mInd) + mean(abs(alphaEst1_approx1-alphat).^2)/alphat^2;
    
    
    %Compute estimate of alpha_t using Approx2 in Eq. (36)
    alphaEst1_approx2 = exp(gammaln(Mvalues(mInd)+1/2)-gammaln(Mvalues(mInd)))^2*q*rho*beta1^2*taup^2./real(z1).^2-sigma2;
    alphaEst1_approx2(alphaEst1_approx2<rho*beta1*taup) = rho*beta1*taup;
    
    %Compute and store the resulting normalized bias and NMSE with Approx 2
    bias_approx2(mInd) = bias_approx2(mInd) + mean(alphaEst1_approx2-alphat)/alphat;
    NMSE_approx2(mInd) = NMSE_approx2(mInd) + mean(abs(alphaEst1_approx2-alphat).^2)/alphat^2;
    
    
    
    %Compute ML estimate of alphat using Theorem 1
    alphaEst1_ML = zeros(size(alphaEst1_approx1));

    for indRel = 1:nbrOfRealizations %Go through all realizations
        
       alphaEst1_ML(indRel) = fminsearch(@(xx) -computeZPDF(xx,z1(indRel),rho,q,beta1,taup,sigma2,Mvalues(mInd)),alphaEst1_approx1(indRel),optimset('Display','off'));

    end
    
    alphaEst1_ML(alphaEst1_ML<rho*beta1*taup) = rho*beta1*taup;
    
    
    %Compute and store the resulting normalized bias and NMSE with ML
    NMSE_ML(mInd) = NMSE_ML(mInd) + mean(abs(alphaEst1_ML-alphat).^2)/alphat^2;
    bias_ML(mInd) = bias_ML(mInd) + mean(alphaEst1_ML-alphat)/alphat;
    
    
end



%%Plot simulation results

figure(1);
hold on; box on;

plot(Mvalues,bias_approx1,'k-.','LineWidth',1);
plot(Mvalues,bias_approx2,'b','LineWidth',1);
plot(Mvalues,bias_ML,'r--','LineWidth',1);

xlabel('Number of Antennas (M)');
ylabel('Estimation Bias');
legend('Approx1','Approx2','ML','Location','NorthEast');
ylim([0 0.5]);


figure(2);
hold on; box on;

plot(Mvalues,NMSE_approx1,'k-.','LineWidth',1);
plot(Mvalues,NMSE_approx2,'b','LineWidth',1);
plot(Mvalues,NMSE_ML,'r--','LineWidth',1);
set(gca,'YScale','log');

xlabel('Number of Antennas (M)');
ylabel('Normalized Mean Squared Error');
legend('Approx1','Approx2','ML','Location','NorthEast');
ylim([1e-2 1]);
