function [Air, Aic, Ari, Aci, Arr, Acc, Arc, Acr] = twoLevelVecSchur(AmatUV, nP, nB, nC, nComp)
% To extract the "Submatrices with one level decompositin"                 Aprl 2,2014/Ajit
% "pn" is subdomain number (partition number)

nI = nP-nB;
nR = nB-nC;
nIR = nP-nC;


% Aii = zeros(2*nI,2*nI); Agg = zeros(2*nB,2*nB);
% Aig = zeros(2*nI,2*nB); Agi = zeros(2*nB,2*nI);
Air = zeros(nComp*nI,nComp*nR); Ari = zeros(nComp*nR,nComp*nI);
Aic = zeros(nComp*nI,nComp*nC); Aci = zeros(nComp*nC,nComp*nI);
Arr = zeros(nComp*nR,nComp*nR); Arc = zeros(nComp*nR,nComp*nC);
Acr = zeros(nComp*nC,nComp*nR); Acc = zeros(nComp*nC,nComp*nC);
%Fi = zeros(1,2*nI); Fg = zeros(1,2*nB); 

for j = 1:nComp
%     Fi((j-1)*nI+1:j*nI) = bvecUV(nP*(j-1)+1:nP*(j-1)+nI);                  
%     Fg((j-1)*nB+1:j*nB) = bvecUV(nP*(j-1)+(nI+1):nP*(j-1)+nP);
    for k = 1:nComp
%     Aii((j-1)*nI+1:j*nI, (k-1)*nI+1:k*nI) = AmatUV(nP*(j-1)+1:nP*(j-1)+nI, nP*(k-1)+1:nP*(k-1)+nI);
%     Agg((j-1)*nB+1:j*nB, (k-1)*nB+1:k*nB) = AmatUV(nP*(j-1)+(nI+1):nP*(j-1)+nP, nP*(k-1)+(nI+1):nP*(k-1)+nP);
%     Aig((j-1)*nI+1:j*nI, (k-1)*nB+1:k*nB) = AmatUV(nP*(j-1)+1:nP*(j-1)+nI, nP*(k-1)+(nI+1):nP*(k-1)+nP);
%     Agi((k-1)*nB+1:k*nB ,(j-1)*nI+1:j*nI) = AmatUV(nP*(k-1)+(nI+1):nP*(k-1)+nP, nP*(j-1)+1:nP*(j-1)+nI);    
    
    Air((j-1)*nI+1:j*nI, (k-1)*nR+1:k*nR) = AmatUV(nP*(j-1)+1:nP*(j-1)+nI, nP*(k-1)+(nI+1):nP*(k-1)+nIR);
    Ari((k-1)*nR+1:k*nR, (j-1)*nI+1:j*nI) = AmatUV(nP*(k-1)+(nI+1):nP*(k-1)+nIR, nP*(j-1)+1:nP*(j-1)+nI);
    Aic((j-1)*nI+1:j*nI, (k-1)*nC+1:k*nC) = AmatUV(nP*(j-1)+1:nP*(j-1)+nI, nP*(k-1)+(nIR+1):nP*(k-1)+nP);
    Aci((k-1)*nC+1:k*nC, (j-1)*nI+1:j*nI) = AmatUV(nP*(k-1)+(nIR+1):nP*(k-1)+nP, nP*(j-1)+1:nP*(j-1)+nI);

    Arc((j-1)*nR+1:j*nR, (k-1)*nC+1:k*nC) = AmatUV(nP*(j-1)+(nI+1):nP*(j-1)+nIR, nP*(k-1)+(nIR+1):nP*(k-1)+nP);
    Acr((k-1)*nC+1:k*nC, (j-1)*nR+1:j*nR) = AmatUV(nP*(k-1)+(nIR+1):nP*(k-1)+nP, nP*(j-1)+(nI+1):nP*(j-1)+nIR);
    Arr((j-1)*nR+1:j*nR, (k-1)*nR+1:k*nR) = AmatUV(nP*(j-1)+(nI+1):nP*(j-1)+nIR, nP*(k-1)+(nI+1):nP*(k-1)+nIR);
    Acc((j-1)*nC+1:j*nC, (k-1)*nC+1:k*nC) = AmatUV(nP*(j-1)+(nIR+1):nP*(j-1)+nP, nP*(k-1)+(nIR+1):nP*(k-1)+nP);
    end
end

% Fi = Fi';
% Fg = Fg';


end

