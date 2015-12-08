function [noisim,nois] =  addgausnois(im,avg,stddev)%151203
% add gaussian white noise with mean "avg" and deviation "stddev" (or root
% of variance) to each component of complex data
% outputs noisy image "noisim" and noise alone "nois"

[numrow,numcol,numimage] = size(im);
imim = imag(im); reim = real(im);
for n = 1:numimage
    nsre(:,:,n) = stddev*randn(numrow,numcol);    
    nsim(:,:,n) = stddev*randn(numrow,numcol);
    ns_reim(:,:,n) = reim(:,:,n) + nsre(:,:,n) +avg;
    ns_imim(:,:,n) = imim(:,:,n) + nsim(:,:,n) +avg;
end
noisim = ns_reim + 1i*ns_imim;
nois = nsre + 1i*nsim;

