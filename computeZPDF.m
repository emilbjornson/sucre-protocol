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
%This function computed the PDF of the observation z_k that is used in 
%Eq. (29). To achieve a numerically stable implementation, the logarithms
%of some of the terms are considered, as further described in Footnote 4.


function pdf = computeZPDF(alphatValues,z,rho,q,beta1,taup,sigma2,M)

%Compute lambda1 and lambda2 as defined in Eqs. (33) and (34)
lambda1 = rho*q*beta1^2*taup^2./(alphatValues+sigma2);
lambda2 = sigma2+q*beta1*taup-lambda1;

%Compute the logarithm of f2 in Eq. (33)
f2ln = - (imag(z).^2./lambda2) -log(pi*lambda2)/2;

%Compute the logarithm of the part of f1 in Eq. (32) that does not depend
%on the summation index n
f1ln = -(real(z).^2/lambda2)*(1-lambda1/(lambda1+lambda2)) - M*log(lambda1)-log(pi*lambda2)/2-gammaln(M) - M*log(1/lambda1+1/lambda2);


%Compute term that appears repeatedly in Eq. (32)
B = real(z)/sqrt(lambda2)*sqrt(lambda1/(lambda1+lambda2));

%Go through the summation in Eq. (33)

pdf = zeros(size(f2ln));

for n = 0:2*M-1
    
    pdf(B>0) = pdf(B>0) + exp( gammaln(2*M) - gammaln(2*M-n) - gammaln(n+1) + (2*M-1-n)*log(B(B>0)) + gammaln((n+1)/2) + log( 1 + (-1)^n*gammainc(B(B>0).^2,(n+1)/2) ) + f2ln(B>0) + f1ln(B>0));
    
    pdf(B<0) = pdf(B<0) + exp( gammaln(2*M) - gammaln(2*M-n) - gammaln(n+1) + (2*M-1-n)*log(B(B<0)) + gammaln((n+1)/2) + log( 1 - gammainc(B(B<0).^2,(n+1)/2) ) + f2ln(B<0) + f1ln(B<0));
    
    if n == (2*M-1)
        pdf(B==0) = pdf(B==0) + exp( gammaln(2*M) - gammaln(2*M-n) - gammaln(n+1) + log(( gamma((n+1)/2) )) + f2ln(B==0) + f1ln(B==0));
    end
    
end

%Remove any imaginary part that appeared due to lack of numerical precision
pdf = real(pdf);
