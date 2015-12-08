%% runGAS.m 151203 mnh *******************************************************
% this code implements the Geometric and Algebraic solutions for bSSFP
% signal demodulation (debanding), as well as a hybrid variance-weighted 
% combination of the two, coined the "Geometric-Algebraic Solution" (GAS)
%
% The code can simulate data, or use experimental data if there are 4 phase
% cycles at deltheta = 0, 90, 180, and 270 degrees.  Various image formats
% are accepted, although the code mayt need to be modified depending on how
% the arrays are arranged
%
% Please refer to the paper for more details: 
% "Combined Geometric and Algebraic Solutions for Removal of bSSFP Banding Artifacts with Performance Comparisons "
% URL to be determined
% Written by: Michael N. Hoff, Jalal B. Andre, and Qing-San Xiang, 2015 
% and to 
% "Banding artifact removal for bSSFP imaging with an elliptical signal model."
% http://onlinelibrary.wiley.com/doi/10.1002/mrm.25098/pdf
% Written by: Qing-San Xiang and Michael N. Hoff, 2014 
%
% Copyright, Michael Hoff, University of Washington, 2015
% 
% GAS is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% GAS is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%% ************************************************************************
clearvars
format long g
realsim = input(' Is the data Experimental (e) or Simulated (s) = ','s');
if realsim == 'e', sim = 0;else sim = 1;end 
cyc = [0 1/4 2/4 3/4]';npc = length(cyc); %phase cycles

%% read in or generate data.  This file desires the generation of a complex
% data array "compdata" with (row,colums, phase cycles, slices)

if sim % simulate data
    %generate dataset parameters
    slind = 1;     % slice index, simdata has presumably only 1 slice
    flip = 40;     flipr = flip*pi/180; %in radians
    TR = 4.2; 
    TE = TR/2;
    % generate data AND random noise
    noisSTD = input('Noise level (ENTER gives standard deviation of 50) = ');
    if isempty(noisSTD),noisSTD = 50;end
    % Contiguous is set up 
    numtis = input('Contiguous (0), or "Tissue"-segmented (1) simulations = ');
    if numtis %segmented tissue
        % for now, hardcode the tissue values
        th_pts = [];
        rowrep = 30; %# rows of data for each tissue
        nc = 192;%# columns of data for all tissues (since they're stacked vertically)
        T1 = [420 240 2400];% three tissues: liver, fat, CSF
        T2 = [50 70 1400];% three tissues: liver, fat, CSF
        [amat,bmat,mmat,theta,compdata,rotmat,noisonly] = sim_bssfp_SegmentTissue(cyc,noisSTD,T1,T2,flipr,TR,rowrep,nc);
    else %contiguously varying data - covers all tissue types by setting 
        % a and b to vary for T_1 and T_2 values in the lit, and freezing flip and M0
        % a and ? both vary on x-axis, choose their granularity
       rowrep = [];
       th_pts = 16;% 
        a_pts = 9;% (so nc = th_pts * a_pts)
        b_pts = 140;%
        [amat,bmat,mmat,theta,compdata,rotmat,noisonly,puredata] = sim_bssfp_contiguous(cyc,noisSTD,th_pts,flipr,a_pts,b_pts);
        %% calculate corresponding T1 and T2
        T2 = -TR./log(amat);
        E1 = (amat.*(1+cos(flipr))- bmat.*(1+amat.^2.*cos(flipr)))./(amat.*(1+cos(flipr))-bmat.*(amat.^2+cos(flipr)));
        T1 = -TR./log(E1);
        T1T2 = T1./T2;  
    end
else % real data
    slind = input('Slice range? (enter for all) = ');% slice index
    % currently this file only reads .mat data, insert your own image reader here for other formats
    [fname,pname] = uigetfile('*.mat','Select *.mat file');
    load(fname)
    % may be able to read key parameters from a DICOM header
    parmfile = input('Load parameters from a DICOM file? (''y'' or ''n''): ','s');
    if parmfile == 'y'
        [fDname,pDname] = uigetfile('*.*','Select DICOM to extract PARM');
        parms = dicominfo(fDname);
        TR = parms.RepetitionTime;
        TE = parms.EchoTime;
        flip = parms.FlipAngle;
        flipr = flip * pi/180;
    end
end

%% reorganize array
[nr,nc,npcc,nsl] = size(compdata);
if npc ~= npcc, error('There appears to be a mismatch with the number of phase cycles'); end
if npc ~= 4, error('This script only runs quad-phase cycles'); end
% set aside slice index as the range and declare
if ~isempty(slind)
    disp(['There are ',num2str(nsl),' total slices!! But I''m only keeping slices '...
        ,num2str(slind(1)),' to ',num2str(slind(end))]); 
else
    slind = 1:nsl; % or if no slices were chosen, keep all slices
    disp(['There are ',num2str(nsl),' total slices!! And I''m keeping them all.']);
end
slrange = slind;
%limit FOV?
rows = 1:nr;   cols = 1:nc; % retain all rows and columns
compdata = compdata(rows,cols,:,slrange);
[nr,nc,~,nsl] = size(compdata);  

%% Complex Sum
CS = squeeze(abs(mean(compdata,3)));%compsum 

%% Geometric Solution
tic
[GS,~] = geo(compdata);%GS

%% Algebraic solution
% real/imag components needs reversal for simulated relative to
% experimental data - if there are problems, try commenting the loop below
% out
if ~sim
    cdata = imag(compdata)+1i*real(compdata);
    clear compdata
    compdata = cdata;   
end

[AS,~] = alg(compdata);

%% GAS composite image
% form a variance-weighted combination of the GS and AS
rng = 2;% 2*rng+1 is the size of the variance ROI
if sim
    expfact = [];
else
    mmat = 9; %for real data, smuggle the expansion factor for smoothing in
end
 
% magnitude GAS
[GAS,varmag,w] = GAScombo(abs(GS),abs(AS),rng,mmat);
toc
var_GS = varmag(:,:,:,1); var_AS = varmag(:,:,:,2);
% real GAS
% [rGAS,varre,rw] = GAScombo(real(GS),real(AS),rng,real(rotmat));
% imaginary GAS
% [iGAS,varim,iw] = GAScombo(imag(GS),imag(AS),rng,imag(rotmat));
% cGAS = rGAS + 1i*iGAS;
% calculate GAS magnitude variance
var_GAS = imvar(abs(GAS),rng,mmat);

% %% Linearized GS (second pass) - not in GAS paper, but recommended for improved SNR if only the GS is employed
% reg_siz = 5; linedpts = [1,3;2,4];
% [LGS,~] = linearizeGS(compdata,GS,reg_siz,linedpts);  

if sim
%% Gold Standards
CSgold = mmat.*(1+(bmat-amat).*(1./sqrt(1-bmat.^2)-1)./bmat);
CSgoldrot = rotmat.*(1+(bmat-amat).*(1./sqrt(1-bmat.^2)-1)./bmat);

%% noise fields
GSmagnois = abs(GS) - mmat;   GScompnois = GS - rotmat; 
ASmagnois = abs(AS) - mmat;   AScompnois = AS - rotmat;
GASmagnois = abs(GAS)-mmat;
CSmagnois = abs(CS) - CSgold;    
% CScompnois = CScomp - CSgoldrot;

%% correlation coefficients
if numtis
     for tis = 1:length(T1)
         corrASGSnois(tis) = corr2(ASmagnois((tis-1)*rowrep+1:rowrep*tis,:),GSmagnois((tis-1)*rowrep+1:rowrep*tis,:));
     end
else
    corrASGS = corr2(ASmagnois,GSmagnois);
end
%% Total Relative Error
if numtis
    for tis = 1:length(T1)
        %section tissue recons and their gold standards
        GS_tis = GS((tis-1)*rowrep+1:rowrep*tis,:);
        AS_tis = AS((tis-1)*rowrep+1:rowrep*tis,:);
        GAS_tis = GAS((tis-1)*rowrep+1:rowrep*tis,:);
        CS_tis = CS((tis-1)*rowrep+1:rowrep*tis,:);
        mmat_tis = mmat((tis-1)*rowrep+1:rowrep*tis,:);   
        CSgold_tis = CSgold((tis-1)*rowrep+1:rowrep*tis,:);
%             LGS_tis = LGS((tis-1)*rowrep+1:rowrep*tis,:);

        %compute corresponding TRE and store in array for each tissue
        TRE_GS(tis) = TotRelError(GS_tis,mmat_tis);
        TRE_AS(tis) = TotRelError(AS_tis,mmat_tis);
        TRE_GAS(tis) = TotRelError(GAS_tis,mmat_tis);
        TRE_CS(tis) = TotRelError(CS_tis,CSgold_tis);
%             TRE_LGS(tis) = TotRelError(LGS_tis,mmat_tis);
        clear GS_tis AS_tis GAS_tis CS_tis mmat_tis CSgold_tis 
    end
else   %for the whole dataset
    %  discard edge pixels
%     frow = 10; lrow = nr-9; fcol = 10; lcol = nc-9;
    %or include them all
    frow = 1; lrow = nr; fcol = 1; lcol = nc;
    GSt = GS(frow:lrow,fcol:lcol);ASt = AS(frow:lrow,fcol:lcol);
    GASt = GAS(frow:lrow,fcol:lcol); CSt = CS(frow:lrow,fcol:lcol);
    % LGSt = LGS(frow:lrow,fcol:lcol); 
    mmatt = mmat(frow:lrow,fcol:lcol);
    CSgoldt = CSgold(frow:lrow,fcol:lcol);
    TRE_GS = TotRelError(GSt(:),mmatt(:));
    TRE_AS = TotRelError(ASt(:),mmatt(:));
    TRE_GAS = TotRelError(GASt(:),mmatt(:));
    TRE_CS = TotRelError(CSt(:),CSgoldt(:));
    %         TRE_LGS = TotRelError(LGSt(:),mmatt(:));
end

%% variance
totvarCS = var(CSmagnois(:));
totvarAS = var(ASmagnois(:));
totvarGS = var(GSmagnois(:));
totvarGAS = var(GASmagnois(:));
end

%% does GAS have the least variance?
varicompare = zeros(nr,nc,nsl);%test if solution is the largest variance
cnt = zeros(nsl,1);
for sl = 1:nsl
    for r = 1:nr
        for c = 1:nc
            if var_GAS(r,c,sl) <= var_GS(r,c,sl) && var_GAS(r,c,sl) <= var_AS(r,c,sl)
                varicompare(r,c,sl) = 1; cnt(sl) = cnt(sl) +1;
            end
        end
    end
end

%% finishing plots
pr = input('Do you want to print images? (1 or 0)'); 
if sim
    [maxmg]=imageplotSIM(compdata,cyc,slrange,slind,numtis,th_pts,amat,bmat,noisSTD,flip,TR,T1,T2,GS,AS,GAS,CS,rowrep,CSgold,mmat,TRE_CS,TRE_GS,TRE_GAS,TRE_AS,var_GS,var_AS,var_GAS,varicompare,pr,cnt);
else
    [maxmg]=imageplotEXP(compdata,cyc,slrange,slind,GS,AS,GAS,CS,var_GS,var_AS,var_GAS,varicompare,pr,cnt);
end

%% print data and TRE
if pr  
    if sim
        save('TRE.mat','TRE_GAS','TRE_GS','TRE_AS','TRE_CS');
    end
end
% save('snr.mat','snrCS','snrGS','snrAS','snrGAS','snrD','lsnrCS','lsnrGS','lsnrAS','lsnrGAS','lsnrD','gsnrCS','gsnrGS','gsnrAS','gsnrGAS','gsnrD');