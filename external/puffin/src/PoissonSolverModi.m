%% A simple solver for Poisson's equation: Modified for specific use by Ajit
% Copyright (C) 2003 Johan Hoffman and Anders Logg.
% Licensed under the GNU GPL Version 2.

%%%%%%% Load the mesh %%%%%%%
tic
disp('FEM-Puffin Poisson Solver Initiated...')

%preprocessor
%autoPET
points = dlmread('../data/points.dat','');
triangles = dlmread('../data/triangles.dat','');
edges = dlmread('../data/edges.dat','');
p = points';
e = edges';
t = triangles';

ck = 1;

%%%%%%% Assemble matrix %%%%%%%
A = AssembleMatrix(p, e, t, 'PoissonModi', [], 0 , ck);

%%%%%%% Assemble vector %%%%%%%
b = AssembleVector(p, e, t, 'PoissonModi', [], 0, ck);

%%%%%%% Solve the linear system %%%%%%%
U = A \ b;

%%%%%%% Compute exact solution %%%%%%%
% u = zeros(size(U));
% for i = 1:size(p,2)
%   x = p(:,i);
%   u(i) = sin(pi*x(1)) * sin(2*pi*x(2));
%   %u(i) = 1 + x(1)^2 + (2*x(2)^2); 
% end
toc
%%%%%%% Compute the error %%%%%%%
%error = U - u;
%enorm = max(abs(error));
%disp(['Maximum norm error: ' num2str(enorm)]);

%%%%%%% Plot solution, exact solution, and error %%%%%%%
figure(2); clf
pdesurf(p,t,U)
shading faceted
title('Computed solution')

% figure(3); clf
% pdesurf(p,t,u);
% shading faceted
% title('Exact solution')

%figure(3); clf
%pdesurf(p,t,error)
%shading faceted
%title('Error')

figure(1)
pdemesh(p,e,t)
axis equal
title('Mesh File')

%%--------------------------------------------------------------------------------------------------------
%%************* END *****************
%%--------------------------------------------------------------------------------------------------------   
