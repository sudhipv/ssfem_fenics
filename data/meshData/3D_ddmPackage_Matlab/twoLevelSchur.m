function [Air, Aic, Ari, Aci, Arr, Acc, Arc, Acr] = twoLevelSchur(Amat, nP, nB, nC)
% To extract the "Submatrices with one level decompositin"                 Aprl 2,2014/Ajit
% "pn" is subdomain number (partition number)

%%-    Date:       Author:   Comments/Modifications:
%   Aprl/02/2014     AD       Original


Air  = Amat(1:(nP-nB),(nP-nB+1):(nP-nC));
Aic  = Amat(1:(nP-nB),(nP-nC+1):nP);


Ari  = Amat((nP-nB+1):(nP-nC),1:(nP-nB));
Aci  = Amat((nP-nC+1):nP,1:(nP-nB));


Arr  = Amat((nP-nB+1):(nP-nC),(nP-nB)+1:(nP-nC));
Acc  = Amat((nP-nC+1):nP,(nP-nC+1):nP);


Arc  = Amat((nP-nB+1):(nP-nC),(nP-nC+1):nP);
Acr  = Amat((nP-nC+1):nP,(nP-nB+1):(nP-nC));

%%% (A_ii)^s =Aii, (A_itau)^s =Aig, (A_taui)^s =Agi,
%%% (A_tautau)^s =Agg, (f_i)^s =Fi & (f_tau)^s =Fg

end

