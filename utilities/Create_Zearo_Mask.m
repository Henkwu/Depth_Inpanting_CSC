function [I,noisyI,M,noisyImean,lmnI] = Create_Zearo_Mask(input,n,i,GT1)

% This function creates images with holes and substructs its mean using
    NYUV2=1;
    GT_mask=input;


    if NYUV2
         
        I=GT1;  
        I1 = cell(1,1);
        I1{1}=I;
        I=I1;
    end
    k2 = (1/(n^2))*ones(n,n);
    k = ones(n,n);
    sz =  cell(1,length(I));
    noisyImean = cell(1,length(I));
    lmnI = cell(1,length(I));
    for image=1:length(I)
    if NYUV2
        M{image}=double((GT_mask(:,:,i)< 0.005) == 0); % masks 50% of the pixels  
        noisyI{image} = double(GT_mask(:,:,i));
    end

       fprintf('Remove mean: %10d\r',image);
       temp = noisyI{image};
       lmn = rconv2(temp,k);
       imn_size = rconv2(M{image},k);
       noisyImean{image} = lmn./(imn_size+eps);  
       temp = M{image}.*(temp - noisyImean{image});
       noisyI{image} = temp;
       
       tempI = I{image};
       lmnI{image} = rconv2(tempI,k2);
       tempI = tempI - lmnI{image};
       I{image} = tempI;
    end
end