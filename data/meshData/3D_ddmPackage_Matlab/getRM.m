function Rmat = getRM(GlobalBN, nBG, nB, pn)
%%% RESTRICTION MATRIX (SCATTER & GATHER OPERATOR)                         Aprl 2,2014/Ajit
%%% "pn" is subdomain number (partition number)

%%-    Date:       Author:   Comments/Modifications:
%   Aprl/02/2014     AD       Original

%tempBn3  =  ['../data/rbnodes00',num2str(pn)];
%tempBn3  =  ['../data/bnodes00',num2str(pn)];
%tempBn3  =  ['../data/bnodes',pad(num2str(pn),3,'left','0')];

tempBn3  =  ['../nbnodes',pad(num2str(pn),4,'left','0')];
%%% NOTE: we need to use re-arranged (nbones[r,c]) interface nodes  
%%% because, we also re-arrange nodes -> Nnodes[i,r,c] which is 
%%% later used for re-assembling final solution vector
bnode12 = dlmread([tempBn3 '.dat']);
Rmat = zeros(nB, nBG);

for i = 1:nB
    Ri = bnode12(i);
    for j = 1:nBG
        Rj = GlobalBN(j);
        if Ri == Rj
            Rmat(i,j) = 1;
        else
            Rmat(i,j) = 0;
        end
    end
end

end

%%%%%%% END %%%%%%%
