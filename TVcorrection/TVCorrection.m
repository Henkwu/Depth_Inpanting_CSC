function y = TVCorrection(x,gamma,n)
% Total variation implemented using the approximate (exact in 1D) equivalence between the TV norm and the l_1 norm of the Haar (heaviside) coefficients.
% load('UDWT')
% [n,J] = quadlength(x);
% addpath('./UDWT')
qmf = MakeONFilter('Haar');

[ll,wc] = mrdwt(x,qmf,n);

wc = SoftThresh(wc,gamma);

y = mirdwt(ll,wc,qmf,n);
