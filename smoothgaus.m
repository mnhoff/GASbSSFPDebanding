function [img_sm] = smoothgaus(img,expfact,gauskern)%151204 mnh
% computes expansion smoothing by convolving with gaussian spatial filter 

%% expand data by expfact(1) x expfact(2) 
img(img==0)=0.001;%replace any null values with a very small number
[bigimg] = expand(img,expfact); %expand image
gauswin = gauskern*expfact;% Gaussian kernel expanded dimension11114
range = (gauswin-1)/2;
[nr,nc,nsl]=size(bigimg);

%% pad expanded image with zeros to prepare for convolution
bigpad = zeros(nr+gauswin(1)-1,nc+gauswin(2)-1,nsl);
bigpad(range(1)+1:nr+range(1),range(2)+1:nc+range(2),:) = bigimg;

%% form a 2D gaussian, where customgauss(gsize, sigmax, sigmay, theta, 
%% offset, factor, center); 
gwin = customgauss(gauswin, gauswin(2)/5, gauswin(1)/5, 0, 0, 1, [0 0]);
gwin = gwin/sum(gwin(:));%normalize 

%% convolve with super image to smooth
% bigimg_s = bigpad; %this should actually be smaller than bigpad by
% (2*range(1), 2*range(2),:,:) due to the convolution
for k = 1:nsl
    bigimg_sm(:,:,k) = conv2(bigpad(:,:,k),gwin,'valid'); 
end
   
%% rebin to original size prior to output
img_sm = rebin(expfact,bigimg_sm);
   