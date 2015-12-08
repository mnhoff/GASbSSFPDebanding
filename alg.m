function [master,regmap] = alg_multi(compdata)%151204 mnh
% compute 4pt  algebraic solution (GS)
[numrow,numcol,numpc,numslice] = size(compdata);
if numpc == 4
    master = zeros(numrow,numcol,numslice);
    regmap = zeros(numrow,numcol,numslice);
    for s = 1:numslice
        for r = 1:numrow
            for c = 1:numcol
                dat = compdata(r,c,:,s);    
                numalg = (dat(1)*dat(3)+1i*dat(1)*dat(3))*(dat(2)-dat(4))+ ...
                    (dat(2)*dat(4)-1i*dat(2)*dat(4))*(dat(1)-dat(3));
                denalg = dat(1)*dat(2)-dat(3)*dat(4)+1i*(dat(2)*dat(3)-dat(1)*dat(4));
                master(r,c,s) = numalg/denalg;
                
                % regularize with CS
                if denalg == 0 || (abs(master(r,c,s)) > max(abs(dat)))
                    master(r,c,s) = mean(dat);
                    regmap(r,c,s) = 1;% note when there is a regularized pixel
                end
            end
        end
    end
else
    error('Your data is formatted incorrectly');
end