function [amat,bmat,mmat,theta,cdata,rotM,nois,puredata] = sim_bssfp_contiguous(cyc,stddev,th_pts,flip,a_pts,b_pts)
%151203  this file generates bSSFP data that varies contiguously, 
% simulating all possible tissue types subjected to a complete range of
% off-resonance ?. Variation is restricted to T_1 and T_2  
% Typical tissue relaxation value ranges hard coded below
% correspond to: 
% T_1 = 200 ? 3000 ms and 
% T_2 = 40 ? 3000 ms 
% at 1.5T (21-25) and flip = 40 correspond to 
% b = 0.1 ? 0.8, set on the vertical axis 
% a = 0.9 ? 0.9999, set on the horizontal axis...
% for each a
% ? = –? to ?, yielding a total horizontal ? variation of 20?. 

npi = 1;%1*pi for -pi to pi, 2*pi for -2pi to 2pi
npc = length(cyc);% # phase cycles
%% base val vectors: M0, theta (th), a, b
M0 = 4000;
th = linspace(-npi*pi,npi*pi*(1-2/th_pts),floor(th_pts));
% a and b range
a1 = 0.91;%1st 'a' value, standard
al = 0.999;%last 'a' value
b1 = 0.795;%1st 'b' value, standard
bl = 0.1;%last 'a' value
ares = (al-a1)/(a_pts-1); % agran = 19 for the GAS paper, meaning resol = 0.001
bres = (bl-b1)/(b_pts-1);
a = [a1:ares:al 0.9999];
b = (b1:bres:bl)';

%% matrices: theta, amat, bmat, and mmat
nc = a_pts*th_pts;
theta = repmat(th,b_pts,a_pts);
amat = zeros(b_pts,nc);
for t = 1:a_pts
    amat(:,th_pts*(t-1)+1:t*th_pts) = repmat(a(t),b_pts,th_pts);
end
bmat = repmat(b,1,nc);

%% calculate M from a, b, and flip 
mmat = M0*bmat.*sin(flip)./(amat*(1+cos(flip)));

%% build the datasets
puredata = zeros(b_pts,nc,npc);
pc = cyc*2*pi;%phase cycles
rotM = mmat.*exp(1i.*theta/2);% generally, M at TE includes an additional theta/2 rotation
for j = 1:npc
    puredata(:,:,j) = rotM.*(1-amat.*exp(-1i.*(theta+pc(j))))./(1-bmat.*cos(theta+pc(j)));
end
%% add gaussian white noise
if ~isempty(stddev)
    [cdata,nois] = addgausnois(puredata,0,stddev);%in:data, mean, noise std dev
end      
