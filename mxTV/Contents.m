% mxTV. 
% Version 1.0  28-april-09. 
% Copyright (c) 2008 by J. Dahl, P.C. Hansen, S.H. Jensen, and T.L. Jensen
%
% Requires Matlab version 7.5 or later versions.
%    
% These functions accompany the paper
%   J. Dahl, P.C. Hansen, S.H. Jensen, and T.L Jensen, "Algorithms and
%   Software for Total Variation Image Reconstruction via First-Order
%   Methods," Numerical Algorithms, xxx.
%
% Demonstration.
%   Pirate.tif     - Test image
%   Text.bmp       - Image with text mask (for inpainting)
%   TVdeblur_demo  - Demonstrates Total Variation image deblurring
%   TVdenoise_demo - Demonstrates Total Variation denoising of an image
%   TVinpaint_demo - Demonstrates Total Variation inpainting of an image
%  
% TV image reconstruction functions.
%   TVdeblur  - Total Variation deblurring
%   TVdenoise - Total Variation denoising
%   TVinpaint - Total Variation inpainting (and denoising)
%
% Files for computations with image blurring.  These files are from
% P.C. Hansen, J.G. Nagy, and D.P. O'Leary, Deblurring Images: Matrices,
% Spectra, and Filtering, SIAM, Philadelphia, 2006.
%   dcts2     - two-dimensional discrete cosine transform
%   dcts      - discrete cosine transform
%   dctshift  - shift operation on PSF matrix
%   idcts2    - two-dimensional inverse discrete cosine transform
%   idcts     - inverse discrete cosine transform
%    
% Auxilary files.
%   info_type_deblur    - Type-like definitions for deblur
%   info_type_denoise   - Type-like definitions for denoise
%   info_type_inpaint   - Type-like definitions for inpaint
%   install_linux       - Installation script for Linux
%   install_windows     - Installating script for Windows
%    
% Directories.
%   c         - The C functions for TV image reconstruction
%   externlib - Library files for the Windows Lcc compiler and fftw3.h
