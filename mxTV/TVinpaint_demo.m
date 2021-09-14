%TVinpaint_demo  Demo script for TV inpainting algorithm
%
% This script illustrates the use of the TV inpainting algorithm
% inplemented in the function TVinpaint; it produces Fig. 3 in the
% accompanying paper. The user can easily modify the script for
% other images, masks, and noise levels.
%
% The script loads a clean image and a mask for the corrupted pixels,
% then adds noise to the ¨image, and finally computes a TV reconstrution.
% Note that the computation of the TV reconstruction takes about a minute.

% J. Dahl^1, P.C. Hansen^2, S.H. Jensen^1 & T.L. Jensen^1
% CSI project, (1) Aalborg University, (2) Technical University of Denmark
% April 28, 2009.

clear, clc
close all
disp('Starting TVinpaint_demo')

% Set parameters.
im_name   = 'Pirate.tif';   % Name of image.
mask_name = 'Text.bmp';     % Name of image with mask; the nonzero
                            % pixels define the mask.
noise_std = 0;             % Standard deviation of image noise.
tau       = 0.85;           % Fudge factor in residual bound.

% Load the clean image and add Gaussian noise; make sure the pixels of the
% noisy image are in the range 0,...,255.
Xc = double(imread('disp.png'));
B = Xc + noise_std*randn(size(Xc));
B(B<0) = 0; B(B>255) = 255; % Pixels in range 0,...,255.

% Load the mask, make sure mask values are 0 and 1, and set all pixels
% that correspond to the mask to white.
% M = imread(mask_name);
      % add inpainting
M=zeros(size(Xc));
M(:,65:70)=1;
M(2:67,50:100)=1;
M(65:67,50:100)=1;
M(95:100,50:100)=1;
M(295:300,5:200)=1;
M(315:320,5:200)=1;
M(15:200,15:200)=1;
M(95:97,95:97)=1;  


% noisyI{image} = M{image}.*I{image};





if ndims(M)==3, M = sum(M,3); end   % Color image -> matrix.
if size(B) ~= size(M), error('Image and mask have different size'), end
M(M>0) = 1;                         % Set mask values to 1.
B(M>0) = 255;                       % Set corrupted pixels to white.

figure(1), clf, colormap(gray)
imagesc(Xc); axis image off
title('Original image')

figure(2), clf, colormap(gray)
subplot(1,2,1)
  imagesc(B); axis image off
  title('Noisy and corrupted image')

% Set the residual bound delta, and compute reconstruction.
disp(['Computing TV reconstruction for tau = ',num2str(tau),...
      ', this takes about a minute ...'])
delta = tau*sqrt(numel(find(M(:)==0)))*noise_std;

[X,info] = TVinpaint(B,M,delta);

subplot(1,2,2)
  imagesc(X); axis image off
  title(['TV inpainted image,  \tau = ',num2str(tau)])
