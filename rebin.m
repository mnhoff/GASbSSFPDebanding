function [small] = rebin(fact,big)
% This will rebin a 'big' matrix into a 'small' one. 
%'fact' = 1x2 array of the shrinking factor
%'nrb'  = the number of rows in the big array
%'ncb'  = the number of cols in the big array
%'nrs'  = the number of rows in the small array
%'ncs'  = the number of cols in the small array

[nrb,ncb,numslice,numblip]=size(big);
nrs = ceil(nrb/fact(1));% this rounds up
ncs = ceil(ncb/fact(2));
small = zeros(nrs,ncs,numslice,numblip);
for b = 1:numblip
    for s = 1:numslice
        for r = 1:nrs
            for c = 1:ncs
                if r == nrs || c == ncs
                    if r == nrs && c == ncs
                        small(r,c,s,b)=mean(mean(big(fact(1)*(r-1)+1 ...
                            :nrb,fact(2)*(c-1)+1:ncb,s,b)));
                    elseif r == nrs
                        small(r,c,s,b)=mean(mean(big(fact(1)*(r-1)+1 ...
                            :nrb,fact(2)*(c-1)+1:fact(2)*c,s,b)));
                    else
                        small(r,c,s,b)=mean(mean(big(fact(1)*(r-1)+1 ...
                            :fact(1)*r,fact(2)*(c-1)+1:ncb,s,b)));
                    end
                else
                    small(r,c,s,b)=mean(mean(big(fact(1)*(r-1)+1 ...
                        :fact(1)*r,fact(2)*(c-1)+1:fact(2)*c,s,b)));
                    
                end
            end
        end
    end
end

