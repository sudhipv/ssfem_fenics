% Main code : By Ajit Desai : July 2015
% purpose: to calculate all the Moments for N-Dimensional Hermite PC basis.
% Cijk =  <Psi_i Psi_j Psi_k>  : Tripple product require for SSFEM
%
% input:
%       ndim = stochastic dimension/number of KLE terms (L)
%       nord = order of PCE (p)
%
% output:
%    pcdata1D: struct containing 1D PC basis information:
%       pcdata1D.ndim   ---   stochastic dimension
%       pcdata1D.nord   ---   PCE order
%       pcdata1D.nclp   ---   number of 1D quadrature points                           
%       pcdata1D.x      ---   1D quadrature nodes
%       pcdata1D.w      ---   1D quadrature weights
%       pcdata1D.multIndex    ---   the multi-indicies 
%       pcdata1D.nPCTerms     ---   number of terms in PCE
% 
%       moment1D: momements of 1D-Hermite polynomials 
%       momentND: Aprod : moments of ND-Hermite polynomials  
%       cijk: <Psi_i Psi_j Psi_k> Tripple product
%       ijk: Non-zero tripple product indices
%
% dependency: 
%       pcdata_1D.m : to setup 1D-Hermite PC basis data
%       quadXW_1D.m : 1D Gauss-Hermite quadrature nodes and weights
%       psi_1D.m    : generates 1-D HERMITE polynomials
%       multiIndex.m: generates the multi-index for the PC basis
%       moment_1D.m : calculates moments for 1D-Hermite PC basis

% NOTE: This package is developed based on Sandia's UQTK tookbox
% "http://www.sandia.gov/UQToolkit/"
% Using reference : "Spectral Methods for Uncertainty Quantification"
% Book by O.P. Le Maître & O.M. Knio
% Appendix-C: Implementation of Product and Moment Formulas

% %%%------------------------------------------------------------------------
% % Basic inputs: 
% ndim = input('enter stochastic dimension (L)  ');
% nord = input('enter order of PCE (p)  ');
% 
% % 1D Hermite PC basis data
% pcdata1D = pcdata_1D(nord, ndim);
% 
% % 1D Moments of Hermite PC basis 
% moment1D = moment_1D(pcdata1D);
% 
% % Required Data
% multiIndex  = pcdata1D.multiIndex;
% nPCTerms  = pcdata1D.nPCTerms;
% apow = moment1D;
% 
% % ND Moments using 1D moments 
% tol = 1e-4;
% 
% nZcijk = [];
% ijk = [];
% cijk = [];
% indexi = 1;
% for i = 1 : nPCTerms
%     nZ = 0;
%     for j = 1 : nPCTerms
%         for k = 1 : nPCTerms
%             aprod = 1;
%             for m = 1 : ndim
%                 l1 = multiIndex(i, m) + 1;
%                 l2 = multiIndex(j, m) + 1;
%                 l3 = multiIndex(k, m) + 1;
%                 aprod=aprod*apow(l1, l2, l3);
%             end  
%             temp1 = aprod;
%             if (temp1 > tol)
%                 nZ = nZ+1;
%                 indexi;
%                 ijk = [ijk, [i;j;k]];
%                 cijk = [cijk, temp1];
%             end
%             indexi = indexi+1;
%         end
%     end
%     nZcijk = [nZcijk;nZ];
% end
% 
% nZijk = [];
% ijk3 = ijk';
% ijk23 = ijk3(:,2:3);
% ijk23 = unique(ijk23,'rows');
% ijk2 = ijk23(:,1);
% for i=1:nPCTerms
%     nZZ = 0;
%     for j=1:length(ijk2)
%         if (i==ijk2(j))
%             nZZ = nZZ+1;
%         end
%     end
%     nZijk = [nZijk;nZZ];
% end
% 
% % setup the output
% fileID = fopen('cijk_approx','w');
% fprintf(fileID,'%d\n',size(ijk,2));
% for i = 1:size(ijk,2)
%    fprintf(fileID,'%d %d %d %17.8E\n',ijk(:,i),cijk(i));
% end
% fclose(fileID);


%% IF you allready have the cijk file use this part of the code to 
% extract non zero structure
ndim = input('enter stochastic dimension (L)  ');
nord = input('enter order of PCE (p)  ');

nPCTerms = factorial(nord+ndim)./(factorial(nord).*factorial(ndim))

str1 = strcat('000',int2str(nord));

if ndim <= 10
    str2 = strcat('000',int2str(ndim));
else
    str2 = strcat('00',int2str(ndim));   
end

str3 = strcat('cijk',str1,str2);
cijk = dlmread(str3);
ijk = cijk(2:end,1:3);

nZijk = [];
ijk23 = ijk(:,2:3);
ijk23 = unique(ijk23,'rows');
ijk2 = ijk23(:,1);
for i=1:nPCTerms
    nZZ = 0;
    for j=1:length(ijk2)
        if (i==ijk2(j))
            nZZ = nZZ+1;
        end
    end
    nZijk = [nZijk;nZZ];
end

str4 = strcat('nZijk',str1,str2,'.dat');
dlmwrite(str4, nZijk, '\t');
%! cp nZijk.dat ../

plot(ijk(:,1),ijk(:,2),'gs',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor',[0.5,0.5,0.5])
set(gca,'YDir','reverse')
axis([0 nPCTerms+1 0 nPCTerms+1])
axis square




