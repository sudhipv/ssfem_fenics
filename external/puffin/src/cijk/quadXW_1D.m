function [x, w] = quadXW_1D(nclp)
% purpose:  returns the 1D Gauss-Hermite quadrature nodes and weights of the
%
% input:
%    nclp     : number of 1D nodes 
%
% output:
%    x, w     : quadrature nodes and weights 
%%%------------------------------------------------------------------------

% computing 1D quadrature points and weights for Gauss-Hermite quadrature
[x, w] = qrule(nclp, 9, 0, 0);
x = x * sqrt(2);
w = w / sqrt(pi);

