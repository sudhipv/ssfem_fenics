function BrSmat = getBrS(m1, nBGr)

%%-Deterministic restriction operators for Mapping Remaining nodes
% m1 = 2;

GlobalRN = dlmread('../remaining_nodes.dat', '');
FinalCount = dlmread('../FinalCount.dat', '');
Count = FinalCount(:,2);

tempBn4  = ['../rnodes000',num2str(m1)];
rnode12 = dlmread([tempBn4 '.dat']);
BrSmat = zeros(length(rnode12), length(GlobalRN));

%nBGr = length(GlobalRN);

for i = 1:length(rnode12)
    Ri = rnode12(i);
    for j = 1:nBGr
        Rj = GlobalRN(j);
        Rc = Count(j); 
        if Ri == Rj && Rc == m1
            BrSmat(i,j) = 1;
        elseif Ri == Rj && Rc ~= m1
            BrSmat(i,j) = -1;
        else
            BrSmat(i,j) = 0;
        end
    end
end

end

%%%%%%% END %%%%%%%

%%% Tried %%%
% % for i = 1:length(rnode12)
% %     Ri = rnode12(i);
% %     for j = 1:nBGr
% %         Rj = GlobalRN(j);
% %         if Ri == Rj
% %             BrSmat(i,i) = 1;
% % %             BrSmat(i,j) = 1;
% % %         else
% % %             BrSmat(i,j) = 0;
% %         end
% %     end
% % end
% % end
