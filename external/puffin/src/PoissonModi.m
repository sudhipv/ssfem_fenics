%% Variational formulation of Poisson's equation: Modified for specific use by Ajit
%%% Copyright (C) 2003 Johan Hoffman and Anders Logg.
%%% Licensed under the GNU GPL Version 2.   

%%: Random Process Case
function integral = PoissonModi(k, u, v, w, du, dv, dw, dx, ds, x, d, t, eq)   %ck changed to k

if eq == 1
    if k == 1
        integral = (ck(x,k) * (du'*dv*dx)) + g(x,d,t)*u*v*ds;
    else
        integral = (ck(x,k) * (du'*dv*dx));
    end
else
    integral = f(x,d,t)*v*dx + (g(x,d,t)*gd(x,d,t) - gn(x,d,t))*v*ds;
end

%%%%%%% RP Case: Stochastic Terms (KLE) %%%%%%%
function y = ck(x,k)

offsets = [0.5; 0.5];             
xdash = x-offsets;

%%: For b1 & b2 = 1.0
multipliers = [0.92184, 0.49248, 0.29374, 0.20437, 0.15575];
omegas = [1.30654, 3.67319, 6.58462, 9.63168, 12.72324];

% %%: For b1 & b2 = 0.1
% %%: NOTE: to use this we need to change multipliers indices 
% multipliers = [0.5614, 0.5196, 0.4649];
% omegas = [2.6276, 5.3073, 8.0671];

%%: Mean and Standard Deviation of Underlaying Gaussian
meang = 0.1;      %%: similar to meanc in Fortran package
sigma = 0.3;      %%: similar to sigma in Fortran package 

%%: ndim is number of underlaying random variables: Dimensions 
ndim = 3;

%%: KLE terms of underlaying Gaussian random variables: b1=b2=1.0
g(1) = sigma*multipliers(1)*multipliers(1)*(cos(omegas(1)*xdash(1))*cos(omegas(1)*xdash(2)));
g(2) = sigma*multipliers(1)*multipliers(2)*(cos(omegas(1)*xdash(1))*sin(omegas(2)*xdash(2)));
g(3) = sigma*multipliers(2)*multipliers(1)*(sin(omegas(2)*xdash(1))*cos(omegas(1)*xdash(2)));
g(4) = sigma*multipliers(1)*multipliers(3)*(cos(omegas(1)*xdash(1))*sin(omegas(3)*xdash(2)));
g(5) = sigma*multipliers(3)*multipliers(1)*(sin(omegas(3)*xdash(1))*cos(omegas(1)*xdash(2)));

%%: Mean of LogNormal Process 
%meanl2 = exp(meang + 0.5*sum(g(1)^2+g(2)^2));
%meanl3 = exp(meang + 0.5*sum(g(1)^2+g(2)^2)+g(3)^2);

% %%: For const assembly vector: F=1;
% if (k==0)             %%: it's not necessory, it call eq==2 only 
%     y = 1;
%%: For stochastic assembly matrix: PCE terms: 2-RV case
if (ndim ==2)
    switch (k)
        case 1
            y = 1;
        case 2
            y = g(1);
        case 3
            y = g(2);
        case 4
            %y = g(1)^2/factorial(2);      %% Standard
            y = g(1)^2/sqrt(factorial(2)); %% Non-Diam
        case 5
            y = g(1)*g(2);
        case 6
            %y = g(2)^2/factorial(2);      %% Standard
            y = g(2)^2/sqrt(factorial(2)); %% Non-Diam
    end
    y = meang*y;
    
%%: For stochastic assembly matrix: PCE terms: 3-RV case    
elseif (ndim ==3)
    switch (k)
        case 1
            y = 1;
        case 2
            y = g(1);
        case 3
            y = g(2);
        case 4
            y = g(3);
        case 5
            %y = g(1)^2/factorial(2);      %% Standard
            y = g(1)^2/sqrt(factorial(2)); %% Non-Diam   
        case 6
            y = g(1)*g(2);
        case 7
            y = g(1)*g(3);
        case 8
            %y = g(2)^2/factorial(2);      %% Standard
            y = g(2)^2/sqrt(factorial(2)); %% Non-Diam
        case 9
            y = g(2)*g(3);
        case 10
            %y = g(3)^2/factorial(2);      %% Standard
            y = g(3)^2/sqrt(factorial(2)); %% Non-Diam
    end
    y = meang*y;
    
%%: For stochastic assembly matrix: PCE terms: 5-RV case        
elseif (ndim ==4)
    switch (k)
        case 1
            y = 1;
        case 2
            y = g(1);
        case 3
            y = g(2);
        case 4
            y = g(3);
        case 5
            y = g(4);            
        case 6
            y = g(1)^2/sqrt(factorial(2)); %% Non-Diam   
        case 7
            y = g(1)*g(2);
        case 8
            y = g(1)*g(3);
        case 9
            y = g(1)*g(4);
        case 10
            y = g(2)^2/sqrt(factorial(2)); %% Non-Diam
        case 11
            y = g(2)*g(3);
        case 12
            y = g(2)*g(4);
        case 13
            y = g(3)^2/sqrt(factorial(2)); %% Non-Diam
        case 14
            y = g(3)*g(4);
        case 15
            y = g(4)^2/sqrt(factorial(2)); %% Non-Diam
    end
    y = meang*y;
    
%%: For stochastic assembly matrix: PCE terms: 5-RV case        
elseif (ndim ==5)
    switch (k)
        case 1
            y = 1;
        case 2
            y = g(1);
        case 3
            y = g(2);
        case 4
            y = g(3);
        case 5
            y = g(4);
        case 6
            y = g(5);            
        case 7
            y = g(1)^2/sqrt(factorial(2)); %% Non-Diam   
        case 8
            y = g(1)*g(2);
        case 9
            y = g(1)*g(3);
        case 10
            y = g(1)*g(4);
        case 11
            y = g(1)*g(5);
        case 12
            y = g(2)^2/sqrt(factorial(2)); %% Non-Diam
        case 13
            y = g(2)*g(3);
        case 14
            y = g(2)*g(4);
        case 15
            y = g(2)*g(5);
        case 16
            y = g(3)^2/sqrt(factorial(2)); %% Non-Diam
        case 17
            y = g(3)*g(4);
        case 18
            y = g(3)*g(5);
        case 19
            y = g(4)^2/sqrt(factorial(2)); %% Non-Diam
        case 20
            y = g(4)*g(5);
        case 21
            y = g(5)^2/sqrt(factorial(2)); %% Non-Diam
    end
    y = meang*y;
end

%%%%%%% Conductivity (penalty factor) %%%%%%%
function y = g(x, d, t)
%y = 1e7;
ident = d;
if ((ident == 1) || (ident == 2) || (ident == 3) || (ident == 4))
         y = 1e50;
else
         y = 0;
end

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

% amp = 45.23;    
% centre(1) = 0.5;
% centre(2) = 0.5;
% y = amp*exp(-200*dot((x-centre(1)),(x-centre(2))));
 y = 1;

%%%%%%% Extra Soucre Terms Tried %%%%%%%
%y = 1;
%y = 5*pi^2*sin(pi*x(1))*sin(2*pi*x(2));
%y = 5*exp(-50.0*(x(1)*x(2)));
%y = 5*pi^2*cos(pi*x(1))*cos(2*pi*x(2));
%y = 1 + x(1)^2 + 2*x(2)^2;
%y = 10*exp(-((x(1) - 0.5)^2 + (x(2) - 0.5)^2) / 0.02);
% centre(1) = 0.5;               
% centre(2) = 0.5;             
% y = 5*exp(-200*dot((x-centre(1)),(x-centre(2))));
%
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

%%-------------------------------------------------------------------------
%%************* END *****************
%%------------------------------------------------------------------------- 


