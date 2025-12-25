
#### July 2015 : Ajit Desai
%%%%---------------------------------------------------------------------------%%%%


PACKAGE-1: Moments for N-Dimensional Hermite PC basis using 1D-Hermite PC basis

Purpose: to calculate the Moments for N-Dimensional Hermite PC basis using 1D-Hermite PC basis 
Cijk =  <Psi_i Psi_j Psi_k>  : Triple product require for SSFEM : 
This saves a lot of time, since we’re only solving 1D-quadrature instead of ND-quadrature. 

MAIN Code : >> moment_ND.m : To use this package just execute this and follow screen instructions. 

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
   
% NOTE: This package is developed base on Sandia's UQTK toolkit
% Using reference : "Spectral Methods for Uncertainty Quantification"
% Book by O.P. Le Maître & O.M. Knio
% Appendix-C: Implementation of Product and Moment Formulas
%%%%---------------------------------------------------------------------------%%%%


PACKAGE-2: Construct ND-Hermite PC basis using 1D-Hermite PC basis

purpose: to calculate all the ND Hermite PC basic parameters data using 1D-Hermite PC basis.
This is useful for PCE when we want to use higher order PC expansions. Saves lot of time. 

Main-Code : >> pcdata_ND.m  : To use this package just execute this code and follow screen instructions.

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

% NOTE: This package is developed base on Sandia's UQTK toolkit
% Using reference : "Spectral Methods for Uncertainty Quantification"
% Book by O.P. Le Maître & O.M. Knio
% Appendix-C: Implementation of Product and Moment Formulas

%%%%---------------------------------------------------%%%%
############################################################