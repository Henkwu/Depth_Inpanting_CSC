function B = jbf_inpaint(D,C,w,sigma_d,sigma_r)% r:range 颜色 d:distance > 1 惩罚参数越大惩罚程度越小
%    D should be a double precision matrix of size NxMx1 (i.e., grayscale)  with normalized values in the 
%    closed interval [0,1]. 
%    C should be similar to D, from which the weights are calculated, with normalized values in the
%   closed interval [0,1]. 
% Pre-compute Gaussian distance weights.
% D = double(D);
% D = (D - min(min(D)))/(max(max(D))-min(min(D)));
D = double(D);
C = double(C);
% C = (C - min(min(min(C))))/(max(max(max(C)))-min(min(min(C))));
C = (C - min(min(C)))/(max(max(C))-min(min(C)));
[X,Y] = meshgrid(-w:w,-w:w);
G = exp(-(X.^2+Y.^2)/(2*sigma_d^2));
% Apply bilateral filter.
dim = size(D);
B = zeros(dim);
for i = 1:dim(1)
   for j = 1:dim(2)      
         % Extract local region.
         iMin = max(i-w,1);
         iMax = min(i+w,dim(1));
         jMin = max(j-w,1);
         jMax = min(j+w,dim(2));
         I = D(iMin:iMax,jMin:jMax); 
         % To compute weights from the color image
         J = C(iMin:iMax,jMin:jMax);
         % Compute Gaussian intensity weights according to the color image
         H = exp(-(J-C(i,j)).^2/(2*sigma_r^2));
         % Calculate bilateral filter response.
         F = H.*G((iMin:iMax)-i+w+1,(jMin:jMax)-j+w+1);
         N=find(I==0);
         if isempty(N)
             B(i,j) = sum(F(:).*I(:))/sum(F(:));  
         else
             a=I(:)~=0;
         B(i,j) = sum(F(:).*I(:))/sum(F(:).*a);   
         end
   end
end
B(find(isnan(B)==1))=0;
B = single(B);
% B = single(B);