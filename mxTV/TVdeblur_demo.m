%TVdeblur_demo  Demo script for TV deblurring algorithm
%
% This script illustrates the use of the TV deblurring algorithm
% implemented in the function TVdeblur; it produces Fig. 4 in the
% accompanying paper. The user can easily modify the script for
% other images, point spread functions, and noise levels.
%
% The scripts loads a clean image, then blurs it and adds noise, and
% finally computes a TV reconstruction.  Note that the computation of
% the TV reconstruction takes several minutes.

% J. Dahl^1, P.C. Hansen^2, S.H. Jensen^1 & T.L. Jensen^1
% CSI project, (1) Aalborg University, (2) Technical University of Denmark
% April 28, 2009.

% The algorithm underlying the blurring is described in P.C. Hansen,
% J.G. Nagy, and D.P. O'Leary, Deblurring Images: Matrices, Spectra,
% and Filtering, SIAM, Philadelphia, 2006.

clear, clc
close all
disp('Starting TVdeblur_demo')

% Set parameters.
im_name   = 'Pirate.tif';   % Name of clean image.
noise_std = 8;              % Standard deviation of image noise.
gauss_std = 3.0;            % St. dev. of Gaussian PSF.
tau = 0.45;                 % Fudge factor in residual bound.

% Load the clean image.
Xc = double(imread(im_name));
[m,n] = size(Xc);

% Construct a Gaussian point spread function.
wx = exp(-((1:m)-ceil(m/2)).^2/(2*gauss_std^2))';
wy = exp(-((1:n)-ceil(n/2)).^2/(2*gauss_std^2))';
PSF = wx*wy'; PSF = PSF/sum(sum(PSF));      % Normalized PSF.

% Compute the DCT spectrum of the doubly symmetric PSF.
e1 = zeros(m,n); e1(1,1) = 1.0;
A = dctshift(PSF,[ceil(m/2),ceil(n/2)]);
S = dcts2(A)./dcts2(e1);

% Blur the clean image and add Gaussian noise; make sure the pixels of the
% noisy image are in the range 0,...,255.
B = idcts2(dcts2(Xc).*S) + randn(m,n)*noise_std;
B(B<0) = 0; B(B>255) = 255;

figure(1), clf, colormap(gray)
imagesc(Xc); axis image off
title('Original image')

figure(2), clf, colormap(gray);
imagesc(B); axis image off
title('Noisy and blurred image')

% Set the residual bound delta and compute reconstruction.
disp(['Computing TV reconstruction for tau = ',num2str(tau),...
      ', this takes several seconds ...'])
delta = tau*sqrt(m*n)*noise_std;
[X,info] = TVdeblur(B,PSF,delta);

figure(3), clf, colormap(gray);
imagesc(X); axis image off
title(['TV deblurred image,  \tau = ',num2str(tau)])
