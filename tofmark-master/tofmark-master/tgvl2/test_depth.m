for i=1:654
disp_gt = double(raw_depths_1(:,:,i));
ours=disp_gt;
d_min = min(ours(ours>0));
d_max = max(ours(ours>0));
a=upsampling_result_norm(:,:,i);
TGV_depth(:,:,i) = a*(d_max-d_min)+d_min;
end
save('TGV_result.mat','TGV_result')