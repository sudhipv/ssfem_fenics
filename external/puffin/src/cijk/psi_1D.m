function psi = psi_1D(x, nord)
% purpose:  Generates 1-D HERMITE polynomials
%
% input:
%     x  : column vector of points/quad points : to calculate pc_basis
%  nord  : order of expansion

%
% output:
%    psi:  polynomial evaluated at points x [length(x) * nord]

%%%------------------------------------------------------------------------
nclp = length(x);
psi = zeros(nclp, nord + 1);
psiD = zeros(nclp, nord + 1);

% % Non-Dimentionalized Hermite PC basis
% % 1-st Order Hermite PC  : it's always fixed to psi(1) = 1
% psi(:, 1) = 1;
% psiD(:, 1) = 1;
% % 2nd order Hermite PC : it's always fixed to psi(2) = x
% if (nord > 0)
%     psi(:, 2) = x;
%     psiD(:, 2) = x;
%     
%     % 3rd and more are solved by recursive formula : 
%     % psi(n) = x psi(n-1) - (n-2) psi(n-2)
%     for i = 3 : nord+1
%         psiD(:, i) = x(:) .* psiD(:, i-1) - (i-2) * psiD(:, i-2)
%         psi(:, i) = psiD(:, i)/sqrt(factorial(i-1))
%     end
% end

%%------------------------------------------------------------------------
%%Standard Hermite PC basis 
% 1-st Order Hermite PC  : it's always fixed to psi(1) = 1
psi(:, 1) = 1;

% 2nd order Hermite PC : it's always fixed to psi(2) = x
if (nord > 0)
    psi(:, 2) = x;
    
    % 3rd and more are solved by recursive formula : 
    % psi(n) = x psi(n-1) - (n-2) psi(n-2)
    for i = 3 : nord+1
        psi(:, i) = x(:) .* psi(:, i-1) - (i-2) * psi(:, i-2);
    end
end

