 function weight = linmin(guide,comp,pr)%151206 find minimum Regional
 % Differential Energy (RDE) of the differences of a pt on a line and a 
 %target pt over a pixel region.  Only neighbours with a max phase 
%  difference of "isop" contribute.  Output weights.
[nr,nc,~] = size(comp);
% initialize vars
midnr = ceil(nr/2);midnc = ceil(nc/2);
numer = 0; denom = 0;
isop = pi/1;
for r = 1:nr
    for c = 1:nc
        %check if the isophase criteria is met for the GS pix. wrt ctr. pix
        if abs(angle(guide(r,c)*conj(guide(midnr,midnc)))) < isop 
            %equation from min. of 2nd derivative (see Xiang, 2014)
            numer = numer + ( conj( comp(r,c,pr(2))-guide(r,c) )*( comp(r,c,pr(2))-comp(r,c,pr(1)) ) + ...
                ( comp(r,c,pr(2))-guide(r,c) )*conj( comp(r,c,pr(2))-comp(r,c,pr(1)) ) );
            denom = denom + ( comp(r,c,pr(1))-comp(r,c,pr(2)) )*conj( comp(r,c,pr(1))-comp(r,c,pr(2)) );
         else
            continue
        end
    end
end
% error checking
if comp(midnr,midnc,pr(1)) == comp(midnr,midnc,pr(2))
    weight = 1.0;
elseif denom == 0
    error ('weight denominator is zero for unequal complex values');
else
    weight = numer / (2*denom);
end
    