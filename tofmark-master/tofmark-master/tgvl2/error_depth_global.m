error=0;
for i=1:654
    i
    gtd_depth=GT_depths_1(:,:,i);
%     raw_depth=gtd_depth;
%    raw_depth2=imresize(raw_depth,[87,112],'nearest');
%    [N1,N2]=size(raw_depth2);
%     MASK=ones(N1,N2);
%     x=floor(N1/2);
%     y=floor(N2/2);
% %     miss_size=10;
%     miss_x=x-miss_size:x+miss_size;
%     miss_y=y-miss_size:y+miss_size;
%     MASK(miss_x,miss_y)=0;   
%     raw_depth2=MASK.*raw_depth2;
%     raw_depth=imresize(raw_depth2,[427,561],'nearest');
%     
    
    
%     raw_depth=raw_depths_1(:,:,i);
    test_depth_inpaint1=depth(:,:,i);
% Mask = ones(size(raw_depth));
% [N_total,M_total]=size(find(raw_depth==0));
% coord=find(raw_depth==0);
% Mask(raw_depth == 0) = 0;
% test_depth_inpaint1(Mask == 1) = gtd_depth(Mask == 1);    %raw_depth用GT_depth替换

%global 2
r1 = [];r2 = [];r3 = [];r4 = [];r5 = [];
% for i=1:11
    M=gtd_depth~=0;        %GT中有零元素 要加的
    depth_pred = test_depth_inpaint1.*M;
    depth_gt = gtd_depth;
    depth_pred(depth_pred <= 0) = eps;
    depth_gt(depth_gt <= 0) = eps;
    [n1,n2]=size(depth_gt);
    N_total=n1*n2;
    theda=1.25;
    rel_error = sum(sum(abs(depth_gt - depth_pred)./depth_gt)) / N_total;

    rmse_log_error = sum(sum((abs(log10(depth_gt) - log10(depth_pred)))) / N_total);
    rmse_error = sqrt(sum(sum((depth_gt - depth_pred).*(depth_gt - depth_pred))) / N_total);
    sigma_mat = max(depth_gt./depth_pred,depth_pred./depth_gt);
    temp = sigma_mat;temp(temp < theda) = 1;temp(temp >= theda) = 0; 
    sigma_1_error = sum(sum(temp)) / (N_total);
    temp = sigma_mat;temp(temp < theda*theda) = 1;temp(temp >= theda*theda) = 0; 
    sigma_2_error = sum(sum(temp)) / (N_total);
    temp = sigma_mat;temp(temp < theda*theda*theda) = 1;temp(temp >= theda*theda*theda) = 0; 
    sigma_3_error = sum(sum(temp)) / (N_total);
    error1=[rel_error,rmse_log_error,rmse_error,sigma_1_error,sigma_2_error,sigma_3_error];
    error=error+error1;
end
error_final=error/654
