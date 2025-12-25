function DrS = getDrS(m1, nBGr)   %%(partition number, Global remainig nodes)
%m1 = 2; 
%nBGr = length(GlobalBN);

%%%%%%% Block Diagonal Matrix %%%%%%%
GlobalRN = dlmread('../remaining_nodes.dat', '');  
ndom = dlmread('../num_partition.dat', '');
Count = zeros(1,nBGr);

for k = 1:ndom 
    tempBn2  = ['../rnodes000',num2str(k)]; %  
    rnode12 = dlmread([tempBn2 '.dat']);
    nBr = length(rnode12);

    for i = 1:nBGr                        % nBG: number of Global BN
        Ri = GlobalRN(i);
        for j = 1:nBr                     % nLB: number of local BN
            Rj = rnode12(j);
            if Ri == Rj 
            Count(i) = Count(i) + 1;
            end
        end
    end
end

%Testing = zeros(nBGr, nBGr);
% for k = m1
    tempBn3  = ['../rnodes000',num2str(m1)]; %['bnodes00',num2str(m1)];
    rnode13 = dlmread([tempBn3 '.dat']);
    DrS = zeros(length(rnode13), length(rnode13));
    
    for i = 1:length(rnode13)
        Ri2 = rnode13(i);
        for j = 1:length(GlobalRN)
            Rj2 = GlobalRN(j);
            if Ri2 == Rj2
                DrS(i,i) = 1/Count(j);
            end
        end
    end
    
end

%%%%%%% END %%%%%%%

%%% TESTING OF DRS MATRIX %%%
% % aa = [];
% % for i=1:16
% % m1 = i;
% % load boundary_nodes.dat;
% % nBG = length(boundary_nodes);
% % R = getRM(m1, nBG);
% % 
% % GlobalRN = dlmread('remaining_nodes.dat', '');
% % nBGr = length(GlobalRN);
% % D = getDrS(m1, nBGr);
% % RDR = R'*D*R;
% % aa = aa + RDR;
% % end