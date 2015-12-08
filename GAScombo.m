function [soln,varmap,w] = GAScombo_multi(d1,d2,varkern,EXPorGOLD) %151204
% GAScombo combines two images using noise-variance-weighted average, where 
%noise variance is found using difference of image and either a gold 
% standard or smoothed version of the image image. Supervariable EXPorGOLD 
% creating the smoothed image for experimental data
gauskern = 5;
[nr,nc,nsl] = size(d1);
[nr1,nc1,nsl1] = size(d2);
if (nr ~= nr1) || (nc ~= nc1) || (nsl ~= nsl1)
    error('Images are not the same size!');
end
% set the GAS to the GS initially
soln = d1;
%% difference maps from gold standard (sim) or smoothed data (exp)
if size(EXPorGOLD)>1 %sim: EXPorGOLD is the GOLD standard (since it's an image)
    err1 = d1-EXPorGOLD; 
    err2 = d2-EXPorGOLD;%
else %exp: EXPorGOLD is an EXPansion factor, since it's only one number
    err1 = d1-smoothgaus(d1,[EXPorGOLD EXPorGOLD],gauskern);
    err2 = d2-smoothgaus(d2,[EXPorGOLD EXPorGOLD],gauskern);
end

%% calculate variance, weights, and average for each pixel
varmap = zeros(nr,nc,nsl,2);% variance maps
w = zeros(nr,nc,nsl,2);% weights
for sl = 1:nsl
    for r = 1:nr
        for c = 1:nc
            % set aside ROIs - define first/last row/column
            fr = max(1,r-varkern); lr = min(nr,r+varkern);
            fc = max(1,c-varkern); lc = min(nc,c+varkern);
            ROI1 = err1(fr:lr,fc:lc);  ROI2 = err2(fr:lr,fc:lc);          

            %calculate variances
            varmap(r,c,sl,1) = var(ROI1(:));
            varmap(r,c,sl,2) = var(ROI2(:));

            %calculate weights
            w(r,c,sl,1)= varmap(r,c,sl,2)^2/(varmap(r,c,sl,1)^2+varmap(r,c,sl,2)^2);
            w(r,c,sl,2) = 1-w(r,c,sl,1);
   
            %solution is the weighted average
            soln(r,c,sl) = w(r,c,sl,1)*d1(r,c,sl) + w(r,c,sl,2)*d2(r,c,sl);
        end
    end
end
