function [nC, nR] = getCorRemDim(pn)
% To extract number of "Corner" and "Remaning" nodes for each partition                             Aprl 2,2014/Ajit
% "pn" is subdomain number (partition number)

%%-    Date:       Author:   Comments/Modifications:
%   Aprl/02/2014     AD       Original

%tempMdCR2  = ['../data/dimcrn00',num2str(pn)]; %% meshdim = [nC nR]
tempMdCR2  = ['../dimcrn',pad(num2str(pn),4,'left','0')];
meshdimCR12 = dlmread([tempMdCR2 '.dat']);

nC = meshdimCR12(1);                   %% number of corner nodes
nR = meshdimCR12(2);

end

