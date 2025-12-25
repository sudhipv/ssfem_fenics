function RrSmat = getRrS(m1, nB)

%%-Deterministic restriction operators for Mapping Remaining nodes  
%%-m1 = 1; 

%tempBn3  =  ['nbnodes00',num2str(m1)];
tempBn3  =  ['../nbnodes',pad(num2str(m1),4,'left','0')];
bnode12 = dlmread([tempBn3 '.dat']);

%tempBn4  = ['rnodes00',num2str(m1)];
tempBn4  = ['../rnodes',pad(num2str(m1),4,'left','0')];
rnode12 = dlmread([tempBn4 '.dat']);
RrSmat = zeros(length(rnode12), length(bnode12));

%%-nB = length(bnode12);

for i = 1:length(rnode12)
    Ri = rnode12(i);
    for j = 1:nB
        Rj = bnode12(j);
        if Ri == Rj
            RrSmat(i,j) = 1;
        else
            RrSmat(i,j) = 0;
        end
    end
end

end

%%%%%%% END %%%%%%%