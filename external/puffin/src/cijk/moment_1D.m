function  moment1D = moment_1D(pcdata1D)
% purpose: computes 1D moments 
% the triple products <psi_i psi_j psi_k> : in case of 3 dimensional RV 
%
% input:
%    pcdata1D: contains basic 1D PC basis data
%
% output: 
%    moment1D: calculate 1D moments 


%%%------------------------------------------------------------------------
% these are extracted for convenience
nord = pcdata1D.nord;
nclp = pcdata1D.nclp;
w    = pcdata1D.w;
psi  = pcdata1D.psi;

apow=zeros(nord + 1, nord + 1, nord + 1);

for k = 1 : nord + 1
    for j = 1 : nord + 1
        for i = 1 : nord + 1
            sum = 0;
            for m = 1 : nclp
                sum = sum + psi(m,i)*psi(m,j)*psi(m,k)*w(m);
            end
            apow(i,j,k) = sum;
        end
    end
end


moment1D = apow;
%moment1D.multiIndex = multiIndex;
%moment1D.nPCTerms = nPCTerms;
%moment1D.ndim = ndim;


% nsparse(1:nPCTerms)=0;
% 
% % tolerance to prune out zero coeffs
% tol = 1e-4;
% 
% for k = 1 : nPCTerms
%     for j = 1 : nPCTerms
%         for i = 1 : nPCTerms
%             aprod = 1;
%             for m = 1 : ndim
%                 l1 = multiIndex(k, m) + 1;
%                 l2 = multiIndex(j, m) + 1;
%                 l3 = multiIndex(i, m) + 1;
%                 aprod=aprod*apow(l1, l2, l3);
%             end    
%             if (aprod > tol)
%                 nsparse(k) = nsparse(k) + 1;
%                 isparse(k, nsparse(k)) = i;
%                 jsparse(k, nsparse(k)) = j;
%                 csparse(k, nsparse(k)) = aprod;
%             end
%         end
%     end
% end
% 
% % setup the output
% pcdata_out = pcdata1D;
% pcdata_out.nsparse = nsparse;
% pcdata_out.isparse = isparse;
% pcdata_out.jsparse = jsparse;
% pcdata_out.csparse = csparse;
