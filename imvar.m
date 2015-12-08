function [varmap] = imvar(inim,varkern,EXPorGOLD)%151205 mnh
% similar to GAScombo, except calculates variance of noise, which is found
% using difference of image and either a gold standard or smoothed version 
% of the image
gauskern = 5;
[nr,nc,nsl] = size(inim);

if size(EXPorGOLD)>1 %sim: EXPorGOLD is the gold standard
    err = inim-EXPorGOLD;
else %exp: EXPorGOLD is the expansion factor
    err = inim-smoothgaus(inim,[EXPorGOLD EXPorGOLD],gauskern);
end
varmap = zeros(nr,nc,nsl);%variance map
for sl = 1:nsl
    for r = 1:nr
        for c = 1:nc
            % set aside ROIs - define first/last row/column
            fr = max(1,r-varkern); lr = min(nr,r+varkern);
            fc = max(1,c-varkern); lc = min(nc,c+varkern);        
            ROI = err(fr:lr,fc:lc);
            varmap(r,c,sl) = var(ROI(:));%calculate and set variance
        end
    end
end