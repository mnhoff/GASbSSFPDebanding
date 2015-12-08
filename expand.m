function [long_out] = expand(short_in,fact)
% this function repeats values "fact" times each along a given direction
%fact(1) is the number of times to expand the number of rows
%fact(2) is the number of times to expand the number of columns 
[nr,nc,nd,ns,nl] = size(short_in);
long_out = zeros(fact(1)*nr,fact(2)*nc,nd,ns,nl);
if fact(1) == 1 && fact(2) == 1
    long_out = short_in;
else
    if fact(1) ~= 1 % row expansion
        for r = 1:nr
            long_out(fact(1)*(r-1)+1:fact(1)*r,1:nc,:,:,:) = ...
                repmat(short_in(r,:,:,:,:),fact(1),1);
        end
        short_in = long_out; %in case you want to expand columns too
    end
    if fact(2) ~= 1 % column expansion
        for c = 1:nc
            long_out(:,fact(2)*(c-1)+1:fact(2)*c,:,:,:) = ...
                repmat(short_in(:,c,:,:,:),1,fact(2));
        end
    end
end

