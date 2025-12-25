function RcSmat = getRcS(m1, nB)

%%-Deterministic restriction operators for Maping Corner nodes  
%%-m1 = 1; 

%tempBn3  = ['nbnodes00',num2str(m1)];% 
tempBn3  =  ['../nbnodes',pad(num2str(m1),4,'left','0')];
bnode12 = dlmread([tempBn3 '.dat']);

%tempBn4  = ['cnodes00',num2str(m1)];
tempBn4  = ['../cnodes',pad(num2str(m1),4,'left','0')];
cnode12 = dlmread([tempBn4 '.dat']);
RcSmat = zeros(length(cnode12), length(bnode12));

%%-nB = length(bnode12);

for i = 1:length(cnode12)
    Ri = cnode12(i);
    for j = 1:nB
        Rj = bnode12(j);
        if Ri == Rj
            RcSmat(i,j) = 1;
        else
            RcSmat(i,j) = 0;
        end
    end
end

end

%%%%%%% END %%%%%%%