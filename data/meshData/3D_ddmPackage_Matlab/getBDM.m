function Dmat = getBDM(m1, nBG)   %%(partition number, Global boundary nodes)
%m1 = 2; 

%%%%%%% Block Diagonal Matrix %%%%%%%
GlobalBN = dlmread('../boundary_nodes.dat', '');  %nBG = length(GlobalBN);
ndom = dlmread('../num_partition.dat', '');
Count = zeros(1,nBG);

for k = 1:ndom 
    %tempBn2  = ['nbnodes00',num2str(k)]; %  
    tempBn2  = ['../nbnodes',pad(num2str(k),4,'left','0')];
    bnode12 = dlmread([tempBn2 '.dat']);   
    nLB = length(bnode12);

    for i = 1:nBG                        % nBG: number of Global BN
        Ri = GlobalBN(i);
        for j = 1:nLB                    % nLB: number of local BN
            Rj = bnode12(j);
            if Ri == Rj 
            Count(i) = Count(i) + 1;
            end
        end
    end
end

%Testing = zeros(nBG, nBG);
% for k = m1
    %tempBn3  = ['nbnodes00',num2str(m1)]; %['bnodes00',num2str(m1)];
    tempBn3  = ['../nbnodes', pad(num2str(m1),4,'left','0')];
    bnode13 = dlmread([tempBn3 '.dat']);
    Dmat = zeros(length(bnode13), length(bnode13));
    
    for i = 1:length(bnode13)
        Ri2 = bnode13(i);
        for j = 1:length(GlobalBN)
            Rj2 = GlobalBN(j);
            if Ri2 == Rj2
                Dmat(i,i) = 1/Count(j);
            end
        end
    end
    
end

%%%%%%% END %%%%%%%
