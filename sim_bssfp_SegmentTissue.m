function [amat,bmat,mmat,theta,cdata,rotM,nois] = sim_bssfp_SegmentTissue(cyc,stddev,T1,T2,flip,TR,rowrep,nc)
%% 151203 vertically stack multiple "tissues" of simulated bSSFP data 
% Off-resonance ? was varied from –2? to 2? horizontally.

%logistics
tis = length(T1);%# tissues
npc = length(cyc);%# phase cycles
nr = tis*rowrep;% # total rows
% parameters
% need a slight theta assymmetry for repetition (-2pi -> 2pi-1/2*nc) 
th = linspace(-2*pi,2*pi*(1-2/nc),floor(nc));
E2 = exp(-TR./T2);
E1 = exp(-TR./T1);
M0 = 4000;
Q = (1-E1)./(1-E1.*cos(flip)-E2.^2.*(E1-cos(flip)));
b = E2.*(1+cos(flip)).*Q;
M = M0.*sin(flip).*Q;

%% build parameter matrices  - stacked
bmat = zeros(nr,nc); amat = zeros(nr,nc); mmat = zeros(nr,nc);  
for dd = 1:tis % for each tissue, stack up the parameter values
    bmat(rowrep*(dd-1)+1:rowrep*dd',1:nc) = repmat(b(dd),rowrep,nc);
    amat(rowrep*(dd-1)+1:rowrep*dd',1:nc) = repmat(E2(dd),rowrep,nc);
    mmat(rowrep*(dd-1)+1:rowrep*dd',1:nc) = repmat(M(dd),rowrep,nc);
    theta = repmat(th,nr,1);
end
%% build the datasets
cdata = zeros(nr,nc,npc);
pc = cyc*2*pi;%phase cycles
rotM = mmat.*exp(1i.*theta/2);% generally, M at TE includes an additional theta/2 rotation
for j = 1:npc
    cdata(:,:,j) = mmat.*exp(1i.*theta/2).*(1-amat.*exp(-1i.*(theta+pc(j))))./(1-bmat.*cos(theta+pc(j)));
end

%% add gaussian white noise
[cdata,nois] = addgausnois(cdata,0,stddev);%in:data, mean, noise std dev
