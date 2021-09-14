
%   TGV2-L2 Depth Image Upsampling Testfile
%
%   Author: David Ferstl
%
%   If you use this file or package for your work, please refer to the
%   following papers:
% 
%   [1] David Ferstl, Christian Reinbacher, Rene Ranftl, Matthias R√ºther 
%       and Horst Bischof, Image Guided Depth Upsampling using Anisotropic
%       Total Generalized Variation, ICCV 2013.
%
%   License:
%     Copyright (C) 2013 Institute for Computer Graphics and Vision,
%                      Graz University of Technology
% 
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation; either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see
%     <http://www.gnu.org/licenses/>.

% use noisy of clean input data
noisy = true;
downsample_method = 'bicubic';

% upsampling factor 1/2/3/4
uf = 0;
factor=1   
% read high-res rgb and depth
% disp_gt = double(imread('test_imgs/art_big.png'));
% gray_gt_ = im2double(rgb2gray(imread('test_imgs/view1.png')));
% for i=[1 3 7 9 10 12 13 14 16 17]
% for i=[26,28,114,60] 
 close  all
    i=1
   
disp_gt = double(GT_depths_1(:,:,i));
% disp_gt = double(raw_depths_1(:,:,i));
gray_gt_ = im2double(rgb2gray(Images_1(:,:,:,i)));


 [N1,N2]=size(disp_gt);
    MASK=ones(N1,N2);
    MASK(rand(N1,N2)<0.7)=0;
   
    disp_gt=MASK.*disp_gt;
    figure,imagesc(disp_gt)
%     raw_depth=imresize(disp_gt,[427,561],'nearest');
% 
% 
%     disp_gt= raw_depth; 
gtd_depth_ori=GT_depths_1(:,:,i);
 GT_depth=gtd_depth_ori;


smin = min(disp_gt(:));
smax = max(disp_gt(:));
    
[M, N] = size(disp_gt);
[Mmb, Nmb] = size(gray_gt_);

dd = [Mmb - M; Nmb - N];
    
gray_gt = gray_gt_(dd(1)/2+1:end-dd(1)/2, dd(2)/2+1:end-dd(2)/2);
disp_res_n=disp_gt;
ours=disp_gt;
bicubic = imresize(disp_res_n, 2^uf, 'bicubic');
bilinear = imresize(disp_res_n, 2^uf, 'bilinear');
nearest_nb = imresize(disp_res_n, 2^uf, 'nearest');

weights = zeros(M,N);
weights(ours > 0) = 1;

% calculate rmse for common interpolation methods
mse_bicubic = sum(sum(sqrt((double(disp_gt)-double(bicubic)).^2))) / numel(disp_gt);
mse_bilinear = sum(sum(sqrt((double(disp_gt)-double(bilinear)).^2))) / numel(disp_gt);
mse_nearest = sum(sum(sqrt((double(disp_gt)-double(nearest_nb)).^2))) / numel(disp_gt);

% normalize input depth map 
d_min = min(ours(ours>0));
d_max = max(ours(ours>0));
    
ours_norm = (ours-d_min)/(d_max-d_min);
ours_norm(ours_norm < 0) = 0;
ours_norm(ours_norm > 1) = 1;

%% tgv l2

timestep_lambda = 1;

tgv_alpha = [17 1.2];

tensor_ab = [9 0.85];

lambda_tgvl2 = 40;

maxits = 1000;
disp(' ---- ');

check = round(maxits/100);
% check = maxits+10;
     
upsampling_result_norm = upsamplingTensorTGVL2(ours_norm, ours_norm, ...
weights.*lambda_tgvl2, gray_gt, tensor_ab, tgv_alpha./factor, timestep_lambda, maxits, ...
check, 0.1, 1);

test_depth_inpaint = upsampling_result_norm*(d_max-d_min)+d_min;
   figure,imagesc(test_depth_inpaint)
PSNR=20*log10((255*sqrt(numel(GT_depth))) / norm(reshape((test_depth_inpaint) - (GT_depth),1,[])))

% ∆¿º€÷∏±Í
    r1 = [];r2 = [];r3 = [];r4 = [];r5 = [];
    N_total=N1*N2;
    theda=1.25;
    depth_pred=test_depth_inpaint;
    depth_gt=GT_depth;
    max_Number=max(max(disp_gt));
    PSNR=20*log10((max_Number*sqrt(numel(GT_depth))) / norm(reshape((test_depth_inpaint) - (GT_depth),1,[])))
    rel_error = sum(sum(abs(depth_gt - depth_pred)./depth_gt)) / N_total;
    rmse_log_error = sum(sum((abs(log10(depth_gt) - log10(depth_pred))))) / N_total;
    rmse_error = sqrt(sum(sum((depth_gt - depth_pred).*(depth_gt - depth_pred)))) / N_total;
%     sigma_mat = max(depth_gt./depth_pred,depth_pred./depth_gt);
    sigma_mat = abs(depth_gt-depth_pred)./depth_gt+eps;
    temp = sigma_mat;temp(temp < theda) = 1;temp(temp >= theda) = 0; 
    sigma_1_error = sum(sum(temp)) / (N_total);
    temp = sigma_mat;temp(temp < theda*theda) = 1;temp(temp >= theda*theda) = 0; 
    sigma_2_error = sum(sum(temp)) / (N_total);
    temp = sigma_mat;temp(temp < theda*theda*theda) = 1;temp(temp >= theda*theda*theda) = 0; 
    sigma_3_error = sum(sum(temp)) / (N_total);
    error1=[PSNR,rel_error,rmse_log_error,rmse_error,sigma_1_error,sigma_2_error,sigma_3_error];
    error=[error1]