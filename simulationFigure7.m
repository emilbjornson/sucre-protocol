%This Matlab script can be used to generate Figure 7, in the article:
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



%Number of inactive UEs in the cell
K0 = 5000;

%Probability that a UE wants to become active in a block
pA = 0.005;

%Probability of sending an RA pilot if a UE wants to become active
pP = 1;

%Number of RA pilot signals
taup = 10;


%Range of number of colliding users per RA pilot
userValues = 0:10;


%Define vector to store probabilities
probabilities = zeros(length(userValues),1);

%Go through all user numbers
for kind = 1:length(userValues)

    %Compute probability according to the binomial distribution in Eq. (1),
    %by using logarithms of each term to gain numerical stability
    probabilities(kind) =  exp(  gammaln(K0+1) - gammaln(userValues(kind)+1) - gammaln(K0-userValues(kind)+1) + (userValues(kind)).*log(pA*pP/taup) + (K0-userValues(kind)).*log((1-pA*pP/taup)) );
    
end



%%Plot simulation results

figure;
box on;

bar(userValues(2:end),probabilities(2:end)/(1-probabilities(1)));
xlabel('Number of UEs per RA pilot');
ylabel('Probability');
ylim([0 0.3]);
