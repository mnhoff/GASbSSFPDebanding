function [master,regmap] = geo(compdata)%151203 mnh
% compute 4pt geometric cross-solution (GS)
[numrow,numcol,numpc,numslice] = size(compdata);
if numpc == 4
    master = zeros(numrow,numcol,numslice);
    regmap = zeros(numrow,numcol,numslice);
    for s = 1:numslice
        for r = 1:numrow
            for c = 1:numcol
                cdat = compdata(r,c,:,s); rdat = real(cdat); idat = imag(cdat);
                numer = (rdat(1)*idat(3)-rdat(3)*idat(1))*(cdat(2)-cdat(4)) ...
                    - (rdat(2)*idat(4)-rdat(4)*idat(2))*(cdat(1)-cdat(3));
                denom = rdat(1)*idat(2)+rdat(2)*idat(3)+rdat(3)*idat(4)+rdat(4)*idat(1)- ...
                    (rdat(1)*idat(4)+rdat(4)*idat(3)+rdat(3)*idat(2)+rdat(2)*idat(1));
                master(r,c,s) = numer/denom;

                % regularize with CS
                if denom == 0 || (abs(master(r,c,s)) > max(abs(cdat)))
                    master(r,c,s) = mean(cdat);
                    regmap(r,c,s) = 1;% note when there is a regularized pixel
                end
            end
        end
    end
else
    error('Your dataset in incorrectly formatted!');
end
        