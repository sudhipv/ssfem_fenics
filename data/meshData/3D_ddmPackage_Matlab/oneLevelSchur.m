function [Aii, Aig, Agi, Agg, Fi, Fg] = oneLevelSchur(Amat, bvec, nP, nB)
% To extract the "Submatrices with one level decompositin"                 Aprl 2,2014/Ajit
% "pn" is subdomain number (partition number)

%%-    Date:       Author:   Comments/Modifications:
%   Aprl/02/2014     AD       Original

Aii  = Amat(1:(nP-nB),1:(nP-nB));         
Aig  = Amat(1:(nP-nB),(nP-nB)+1:nP);
Agi  = Amat((nP-nB+1):nP,1:(nP-nB));
Agg  = Amat((nP-nB+1):nP,(nP-nB+1):nP);
Fi   = bvec(1:(nP-nB)); 
Fg   = bvec((nP-nB+1):nP);

%%% (A_ii)^s =Aii, (A_itau)^s =Aig, (A_taui)^s =Agi,
%%% (A_tautau)^s =Agg, (f_i)^s =Fi & (f_tau)^s =Fg

end

