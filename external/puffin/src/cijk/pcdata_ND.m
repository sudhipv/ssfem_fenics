% function pcdataND = pcdata_ND(pcdata1D)
% purpose: to calculate all the ND Hermite basic PC parameters data
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
%    pcdataND: struct containing basic PC information:
%       pcdataND.ndim      ---   stochastic dimension
%       pcdataND.nord      ---   PCE order                              
%       pcdataND.nPCTerms  ---   number of terms in PCE
%       pcdataND.Psi       ---   N-D hermite polynomials
% 
% dependency: 
%       pcdata_1D.m : to setup 1D-Hermite PC basis data
%       quadXW_1D.m : 1D Gauss-Hermite quadrature nodes and weights
%       psi_1D.m    : generates 1-D Hermite polynomials
%       multiIndex.m: generates the multi-index for the PC basis
%       psi_ND.m    : generates N-D Hermite polynomials 

% NOTE: This package is developed base on Sandia's UQTK tookbox
% Using reference : "Spectral Methods for Uncertainty Quantification"
% Book by O.P. Le Maître & O.M. Knio
% Appendix-C: Implementation of Product and Moment Formulas

%%%------------------------------------------------------------------------
% Basic inputs: 
ndim = input('enter stochastic dimension (L)  ');
nord = input('enter order of PCE (p)  ');

% 1D Hermite polynomial PC-Basis data
pcdata1D = pcdata_1D(nord, ndim);

% ND Hermite polynomials at collocation points(quad-points)
pcdataND.Psi = psi_ND(pcdata1D);

pcdataND.ndim = ndim;
pcdataND.nord = nord;
pcdataND.nPCTerms = pcdata1D.nPCTerms;

