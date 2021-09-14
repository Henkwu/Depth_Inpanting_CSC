% A Demo script for the image inpainting using the algorithm.
close all
addpath('./functions')
addpath ./mexfiles;
addpath ./image_helpers;
addpath('./vlfeat/toolbox');
addpath('./utilities');
addpath(genpath('./spams-matlab'));
addpath(genpath('./mxTV'))
addpath(genpath('./TVcorrection'))
addpath(genpath('./Dataset'));
load('D_initial.mat')
load('NYU_depth_v2.mat')
load('NYU_75_75.mat')
vl_setup();
TV=1;
bilatral_filter=1;
n = 8; % patch size 8
m = 81; % number of filters

D_initial=D_LoBCoD;

for i=1:100  
    i
    close all
    clear I
    I{1}=double(GT_mask_75_75(:,:,i));
    Image=double(rgb2gray(Images_1(:,:,:,i)));  
    input=GT_mask_75_75;
    [I,noisyI,M,noisyImean,lmnI] = Create_Zearo_Mask(input,n,i,I{1});

    D = D_initial;
 
    params = [];
    params.Y_Original = I;
    params.Y_noisy = noisyI;
    params.lmnI = lmnI;
    params.noisyImean = noisyImean; 
    params.M = M;
    params.lambda = 0.01;
    params.MAXITER = 10; % maximun number of epochs 5
    params.D = D;
    params.Train_on = true(1); % false(1)
    %% Run Convolutional Sparse Coding 
    [cleanI_inpainting,objective_inpainting,avgpsnr_inpainting,sparsity_inpainting,...
    totTime_inpainting,alpha_inpainting,D_inpainting] = Inpainting_LoBCoD(params,TV);


%% Plot

%     figure(1); 
%     subplot(1,3,1);
%     plot(totTime_inpainting,objective_inpainting,'.-b');
%     xlabel('Time [Seconds]','fontsize',10)
%     ylabel('Objective','fontsize',10)
%     xlim([0 400])
%     legend('LoBCoD inpainting')
%     grid on
% 
%     subplot(1,3,2);
%     plot(totTime_inpainting,avgpsnr_inpainting,'.-b');
%     xlabel('Time [Seconds]','fontsize',10)
%     ylabel('Average PSNR','fontsize',10)
%     xlim([0 400])
%     legend('LoBCoD')
%     grid on
% 
%     subplot(1,3,3);
%     plot(totTime_inpainting,sparsity_inpainting,'.-b');
%     xlabel('Time [Seconds]','fontsize',10)
%     ylabel('Sparsity','fontsize',10)
%     xlim([0 400])
%     legend('LoBCoD')
%     grid on
% 
%     figure(2); 
%     showDictionary(D_inpainting);
%     title('The trained dictionary on the corrupted images');
% 
% 
%     figure(3);
%     subplot(1,3,1);
%     % imagesc(I{1}+ lmnI{1}); colormap gray; axis off
%     imagesc(I{1}+ lmnI{1});  axis off
%     title('Original image')
% 
%     subplot(1,3,2);
%     imagesc(noisyI{1}+ M{1}.*noisyImean{1}); axis off
%     title('Corrupted image')
% 
%     subplot(1,3,3);
%     imagesc(cleanI_inpainting{1}+ noisyImean{1});  axis off
%     title(['Clean image PSNR = ',num2str(avgpsnr_inpainting(end),4),'dB']
    
%% nan
    test_depth_inpaint=cleanI_inpainting{1}+ noisyImean{1};
    test_depth_inpaint(find(isnan(test_depth_inpaint)==1)) = 0;
    GT_depth=I{1}+ lmnI{1}; 
    %% Hierarchical joint bilateral filter
    if bilatral_filter
    raw_depth_rec=cleanI_inpainting{1}+ noisyImean{1};
    raw_depth_rec(find(raw_depth_rec<0.5))=0;
    Mask_CSC=double(raw_depth_rec==0);
    raw_depth_ori=raw_depth_rec;  
        for j=1:40
            
            kernel=7;
            test_depth_inpaint = jbf_inpaint(raw_depth_ori, Image, kernel, 15, 0.8);  
            test_depth_inpaint=raw_depth_rec+Mask_CSC.*test_depth_inpaint;
            test_N=find(test_depth_inpaint==0);
            if isempty(test_N)                
                break
            else
                raw_depth_ori=test_depth_inpaint;
                continue 
            end
        end
    end
%% Evaluation of the benchmark
    r1 = [];r2 = [];r3 = [];r4 = [];r5 = [];    
    theda=1.25;
    depth_gt = double(GT_depths_1(:,:,i));
    N_total=length(find(M{1}==0));
    depth_pred=double(~M{1}).*test_depth_inpaint;
    depth_gt=double(~M{1}).*double(depth_gt);  
    sigma_mat = max(depth_gt./(depth_pred+eps),depth_pred./(depth_gt+eps));
    temp = sigma_mat;temp(temp < theda) = 1;temp(temp >= theda) = 0; 
    sigma_1_error = sum(sum(double(~M{1}).*temp)) / (N_total);
    temp = sigma_mat;temp(temp < theda*theda) = 1;temp(temp >= theda*theda) = 0; 
    sigma_2_error = sum(sum(double(~M{1}).*temp)) / (N_total);
    temp = sigma_mat;temp(double(~M{1}).*temp < theda*theda*theda) = 1;temp(temp >= theda*theda*theda) = 0; 
    sigma_3_error = sum(sum(double(~M{1}).*temp)) / (N_total); 
    max_Number=max(max(depth_pred));
    rel_error = sum(sum(abs(depth_gt - depth_pred)./(depth_gt+eps))) / N_total;
    rmse_log_error = sum(sum((abs(log10(depth_gt+eps) - log10(depth_pred+eps))))) / N_total;
    rmse_error = sqrt(sum(sum((depth_gt - depth_pred).*(depth_gt - depth_pred))) / N_total);
    error(:,:,i)=[rel_error,rmse_log_error,rmse_error,sigma_1_error,sigma_2_error,sigma_3_error];
end
   
%% Average 
error_avag=0;
for i=1:100
    error_avag=error(:,:,i)+error_avag;
end
error_avag;
