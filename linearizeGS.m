function [LGS,negmap] = linearizeGS(cdata,GS,reg_siz,pr)%151205 mnh 
% - inputs fully phase cycled data and the GS to perform a 2nd-pass, 
% increasing SNR.  The solution is constrained to lie on a spoke connecting
% deltheta = 180° phase cycled datasets (pr matches alternative phase
% cycles).  The exact location along the spoke is found from the minimum of
% the distance to the GS over a local region of pixels.  Both such
% solutions for each spoke are averaged for the final solution
[nr,nc,~,nsl] = size(cdata);
LGS1 = GS; %initialize line 1 soln to GS
LGS2 = GS; %initialize line 2 soln to GS

%negmap tells us when we have spurious negative weights coming out of the
%linear minimizer
negmap = zeros(nr,nc,nsl,2);
%using reg-size for the range, introduce some ROI indexing variables 
start = ceil(reg_siz/2);
disp = floor(reg_siz/2);
% pixel by pixel, linearize
for sl = 1:nsl
    for r = start:nr-disp
        for c = start:nc-disp
            %for each soln, use the 2nd derivative to find the ROI min from
            %the GS in the form of weights for a weighted sum
            w1 = linmin(GS(r-disp:r+disp,c-disp:c+disp,sl),cdata(r-disp:r+disp,c-disp:c+disp,:,sl),pr(1,:));
            w2 = linmin(GS(r-disp:r+disp,c-disp:c+disp,sl),cdata(r-disp:r+disp,c-disp:c+disp,:,sl),pr(2,:));

            % calculate the weighted sum 
            LGS1(r,c,sl) = w1.*cdata(r,c,pr(1,1),sl)+(1-w1).*(cdata(r,c,pr(1,2),sl));
            LGS2(r,c,sl) = w2.*cdata(r,c,pr(2,1),sl)+(1-w2).*(cdata(r,c,pr(2,2),sl));
            if w1 < 0 %negative weights?
                negmap(r,c,1) = 1;
            end
            if w2 < 0
                negmap(r,c,2) = 1;
            end
        end
    end
end
LGS = (LGS1+LGS2)/2; % final complex solution
