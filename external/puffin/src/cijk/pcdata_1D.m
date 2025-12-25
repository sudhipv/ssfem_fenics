function pcdata1D = pcdata_1D(nord, ndim)
% purpose: sets up all the 1D Hermite PC basis parameters
%
% input:  
%    nord: order of the PCE : p
%    ndim: the stochastic dimension : L (number of KLE terms)
%
% output: 
%    pcdata: struct containing 1D PC basis information:
%       pcdata1D.ndim   ---   stochastic dimension
%       pcdata1D.nord   ---   PCE order
%       pcdata1D.nclp   ---   number of 1D quadrature points                           
%       pcdata1D.x      ---   1D quadrature nodes
%       pcdata1D.w      ---   1D quadrature weights
%       pcdata1D.multIndex    ---   the multi-indicies 
%       pcdata1D.nPCTerms     ---   number of terms in PCE

% basic data structure 
pcdata1D.ndim = ndim;
pcdata1D.nord = nord;
pcdata1D.nclp = 2 * nord + 1;          % done as 2n+1 quad points 

% computing 1D quadrature points and weights
[x, w] = quadXW_1D(pcdata1D.nclp);
pcdata1D.x = x;
pcdata1D.w = w;

% one-dimensional Hermite polynomials at collocation points(quad-points) 
pcdata1D.psi = psi_1D(x, nord);

% multi-indices and the number of terms in PCE
[pcdata1D.multiIndex, pcdata1D.nPCTerms] = multiIndex(nord, ndim);

