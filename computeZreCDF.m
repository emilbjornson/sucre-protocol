%This Matlab function is used in the simulations of the following article:
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
%
%This function computed the CDF of real part of the the observation z_k
%based on expression in Eq. (39). To achieve a numerically stable
%implementation, the logarithms of some of the terms are considered, as
%further described in Footnote 4.


function cdf = computeZreCDF(alphat,b,rho,q,beta1,taup,sigma2,M)

%Compute lambda1 and lambda2 as defined in Eqs. (33) and (34)
lambda1 = rho*q*beta1^2*taup^2./(alphat+sigma2);
lambda2 = sigma2+q*beta1*taup-lambda1;

%Compute first term in Eq. (39)
cdf = qfunc(-b*sqrt(2/lambda2));


%Compute term that appears repeatedly in Eq. (39)
B = b/sqrt(lambda2)*sqrt(lambda1/(lambda1+lambda2));

%Go through the two sums in Eq. (39)
for k = 0:M-1
    
    factor1ln = -(b.^2/lambda2)*(1-lambda1/(lambda1+lambda2)) - k*log(lambda1)-log(pi*lambda2)/2-gammaln(k+1) - (k+1/2)*log(1/lambda1+1/lambda2);
    
    for n = 0:2*k

        cdf(b>0) = cdf(b>0) - exp( gammaln(2*k+1) - gammaln(2*k+1-n) - gammaln(n+1) + (2*k-n)*log(B(b>0)) + gammaln((n+1)/2) + log( 1 + (-1)^n*gammainc(B(b>0).^2,(n+1)/2) ) + factor1ln(b>0) - log(2));
        
        cdf(b<0) = cdf(b<0) - exp( gammaln(2*k+1) - gammaln(2*k+1-n) - gammaln(n+1) + (2*k-n)*log(B(b<0)) + gammaln((n+1)/2) + log( 1 - gammainc(B(b<0).^2,(n+1)/2) ) + factor1ln(b<0) - log(2));
        
        if n == (2*k)
            cdf(b==0) = cdf(b==0) - exp( gammaln(2*k+1) - gammaln(2*k+1-n) - gammaln(n+1) + gammaln((n+1)/2) + factor1ln(b==0) - log(2));
        end
        
    end
    
end

%Remove any imaginary part that appeared due to lack of numerical precision
cdf = real(cdf);
