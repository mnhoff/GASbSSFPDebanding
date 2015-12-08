function TRE = TotRelError(data,gold)% 151206 calculate the total 
% relative error TRE of an image from a gold standard 
% see eq [14] from Chang & Xiang, MedPhys 2006, SPEED-ACE
[nrg,ncg] = size(gold);         [nr,nc] = size(data);
if nrg ~= nr || ncg ~= nc
    error('The image and it''s gold standard are different sizes!');
end
data = abs(data(:));            gold = abs(gold(:));
dif = data - gold; % differences
SOSres = sum(dif.^2);% sum-of-sqares differences
TRE = sqrt(SOSres)/sum(gold); 