function BcSmat = getBcS(m1, nBGc)

%%-Deterministic restriction operators for Mapping Remaining nodes
%%- m1 = 1;

GlobalCN = dlmread('../corner_nodes.dat', '');

%tempBn4  = ['../data/cnodes00',num2str(m1)];
tempBn4  = ['../cnodes',pad(num2str(m1),4,'left','0')];
cnode12 = dlmread([tempBn4 '.dat']);
BcSmat = zeros(length(cnode12), length(GlobalCN));

%%-nBGc = length(GlobalCN);

for i = 1:length(cnode12)
    Ri = cnode12(i);
    for j = 1:nBGc
        Rj = GlobalCN(j);
        if Ri == Rj
            BcSmat(i,j) = 1;
        else
            BcSmat(i,j) = 0;
        end
    end
end

end

%%%%%%% END %%%%%%%