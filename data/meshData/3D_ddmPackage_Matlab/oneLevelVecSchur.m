function [Aii, Aig, Agi, Agg, Fi, Fg] = oneLevelVecSchur(AmatUV, bvecUV, nP, nB, nComp)
% To extract the "Submatrices with one level decompositin"       Dec 27,2017/Ajit
% "pn" is subdomain number (partition number)

%%-    Date:       Author:   Comments/Modifications:
%   Aprl/02/2014     AD       Original

nI = nP-nB;

Aii = zeros(2*nI,2*nI); Agg = zeros(2*nB,2*nB);
Aig = zeros(2*nI,2*nB); Agi = zeros(2*nB,2*nI);
Fi = zeros(1,2*nI); Fg = zeros(1,2*nB); 

for j = 1:nComp
    Fi((j-1)*nI+1:j*nI) = bvecUV(nP*(j-1)+1:nP*(j-1)+nI);                  
    Fg((j-1)*nB+1:j*nB) = bvecUV(nP*(j-1)+(nI+1):nP*(j-1)+nP);
    for k = 1:nComp
    Aii((j-1)*nI+1:j*nI, (k-1)*nI+1:k*nI) = AmatUV(nP*(j-1)+1:nP*(j-1)+nI, nP*(k-1)+1:nP*(k-1)+nI);
    Agg((j-1)*nB+1:j*nB, (k-1)*nB+1:k*nB) = AmatUV(nP*(j-1)+(nI+1):nP*(j-1)+nP, nP*(k-1)+(nI+1):nP*(k-1)+nP);
    Aig((j-1)*nI+1:j*nI, (k-1)*nB+1:k*nB) = AmatUV(nP*(j-1)+1:nP*(j-1)+nI, nP*(k-1)+(nI+1):nP*(k-1)+nP);
    Agi((k-1)*nB+1:k*nB ,(j-1)*nI+1:j*nI) = AmatUV(nP*(k-1)+(nI+1):nP*(k-1)+nP, nP*(j-1)+1:nP*(j-1)+nI);
    end
end

Fi = Fi';
Fg = Fg';

% Aii  = Amat(1:(nP-nB),1:(nP-nB));         
% Aig  = Amat(1:(nP-nB),(nP-nB)+1:nP);
% Agi  = Amat((nP-nB+1):nP,1:(nP-nB));
% Agg  = Amat((nP-nB+1):nP,(nP-nB+1):nP);
% Fi   = bvec(1:(nP-nB));                   
% Fg   = bvec((nP-nB+1):nP);

%%% (A_ii)^s =Aii, (A_itau)^s =Aig, (A_taui)^s =Agi,
%%% (A_tautau)^s =Agg, (f_i)^s =Fi & (f_tau)^s =Fg

end

