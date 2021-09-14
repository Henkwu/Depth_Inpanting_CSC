close all 
a=imread('fingerprint.png');
x=double(a);
y=uint8(TVCorrection(x,2,1));
figure,imshow(a)
figure,imshow(y)