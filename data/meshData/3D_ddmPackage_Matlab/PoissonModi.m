%% Variational formulation of Poisson's equation: Modified for specific use by Ajit
%%% Copyright (C) 2003 Johan Hoffman and Anders Logg.
%%% Licensed under the GNU GPL Version 2.

function integral = PoissonModi(u, v, w, du, dv, dw, dx, ds, x, d, t, eq)

disp(x)

if eq == 1
  integral = du'*dv*dx + g(x,d,t)*u*v*ds;
else
  integral = f(x,d,t)*v*dx + (g(x,d,t)*gd(x,d,t) - gn(x,d,t))*v*ds;
end

%%%%%%% Conductivity (penalty factor) %%%%%%%
function y = g(x, d, t)
y = 1e7;

% ident = d;
% if ((ident == 1) || (ident == 2) || (ident == 3) || (ident == 4) || (ident == 5))
%          y = 1e7;
% else
%          y = 0;
% end


%%%%%%% Dirichlet boundary condition %%%%%%%
function y = gd(x, d, t)
y = 0;
% y = 1 + x(1)*x(1) + 2*x(2)*x(2);
% y = sin(5*x(1));

%%%%%%% Neumann boundary condition %%%%%%%
function y = gn(x, d, t)
y = 0;
% y = 1 + x(1)*x(1) + 2*x(2)*x(2);
% y = sin(5*x(1));

%%%%%%% Right-hand side, source term %%%%%%%
function y = f(x, d, t)
%y = 5*pi^2*sin(pi*x(1))*sin(2*pi*x(2));

%%%%%%% Extra Soucre Terms Tried %%%%%%%
%y = 5*exp(-50.0*(x(1)*x(2)));
%y = 5*pi^2*cos(pi*x(1))*cos(2*pi*x(2));
%y = 1 + x(1)^2 + 2*x(2)^2;
y = -6;
%y = 10*exp(-((x(1) - 0.5)^2 + (x(2) - 0.5)^2) / 0.02);
%centre(1) = 0.5;               
%centre(2) = 2.5;             
%y = 5*exp(-200*dot((x-centre(1)),(x-centre(2))));
%y = 5*3.14*3.14*sin(3.14*x(1))*sin(2*3.14*x(2));

%%% physical membrane problem
%T = 10.0;  % tension
%A = 1.0;   % pressure amplitude
% R = 0.3;   % radius of domain
% theta = 0.2;
% x0 = 0.6*R*cos(theta);
% y0 = 0.6*R*sin(theta);
% sigma = 0.025;
% % sigma = 50;  % large value for verification
% y = 4*exp(-0.5*(((R*x(1) - x0)/sigma)^2)-0.5*(((R*x(2) - y0)/sigma)^2));

%%--------------------------------------------------------------------------------------------------------
%%************* END *****************
%%--------------------------------------------------------------------------------------------------------   


