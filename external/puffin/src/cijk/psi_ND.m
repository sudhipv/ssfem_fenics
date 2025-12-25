function Psi = psi_ND(pcdata1D)
% multidimensional Hermite polynomials using 1D Hermite polynomials 
%           evaluated at xi.
% Ii is based on : Psi_k(xi_N) = Prod_1toN{pis(multiIndices}

% input paremeters:
% data1D  : struct containing 1D PC basis data
%     xi  : column vector of points/quad points : to calculate pc_basis 
%
% output:
%    Psi  : a vector with Psi(k) = \Psi_k(xi), k=1, ..., 1+P
%
% note: Here the modes of the PC expansion are 
%       numbered 1, ..., P+1 because Matlab array 
%       index must begin with 1. 
%       Hence, 1 is added to the multiIndices     

%%%------------------------------------------------------------------------
% check input
% xi = pcdata1D.x;
% if length(xi) ~= pcdata1D.ndim
%    fprintf('error: input vector must be %i dimensional\n', pcdata1D.ndim);
%    Psi = 0; 
%    return
% end

% get the 1D psi's: Hermite Polynomials  
% psi = psi_1D(x, pcdat1D.nord);
psi = pcdata1D.psi; 

ndim = pcdata1D.ndim;
multiIndex = pcdata1D.multiIndex;
nPCTerms = pcdata1D.nPCTerms;

% get the mult-dim Psi's : Hermite Polynomials 
Psi = ones(1, nPCTerms);
for i = 1 : ndim
    Psi = Psi .* psi(i, multiIndex(1 : nPCTerms, i) + 1);
end