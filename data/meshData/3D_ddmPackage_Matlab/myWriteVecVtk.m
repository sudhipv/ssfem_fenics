%% Code to write Vector valued PDE solution to ".vtk" file
%
% Copyright (C) 2017 Ajit Desai, Ph.D. Candidate ajit.ndesai@gmail.com
%
% The output can be visualized using Paraview
%%-----------------------------------------------------------------------------------------------

clearvars

%%: Load points (FEM nodal co-ordinates) data
points = dlmread('../points.dat','');   %dlmread('points.dat','');
x = points(:,1);
y = points(:,2);
z = points(:,3);
np = length(x);

%%: Load tetrahedrons (FEM nodal connectivity) data
tetra = dlmread('../tetrahedrons4c.dat','');  

%%: Vector valued PDE outputs (From DDM solver)
% solution = dlmread('../../solution/final_solutionMatlab.dat','');
solution = dlmread('../../solution/final_solu.dat','');

%%: If have assembly mat-vec you can solve directly here
%load('/Users/ajit/Downloads/FenicsTests/det-fem/2D_linearElasticity/outputs/Ab.mat')
%solution = An\bn';

%%: make up a "vector" field: np = number of nodes: Solution = [Ux, Uy, Uz]
U = solution(1:np);
V = solution(np+1:2*np);
W = solution(2*np+1:3*np);

% setup filename to write output VTK
 filename = '../../vtkOutputs/python_elasticity.vtk';  %% Binary (use paraview)
% filename = '../../vtkOutputs/matlab_elasticity.vtk';  %% Binary (use paraview)
data_title = 'displacements';

% organize data
% .type field must be 'scalar' or 'vector'
% number of vector components must equal dimension of the mesh
data_struct.type = 'vector';
data_struct.name = 'displacements';
data_struct.data = [U,V,W];   %% Need to provide nodal components

% if *all data* is oriented column-wise, i.e. lots of rows, few columns,
% then set this to *false*
flipped = false;
% otherwise, if you want to transpose the data, then set this to *true*

% write the file
tic
stat = vtk_write_tetrahedral_grid_and_data(filename,data_title,points,tetra,data_struct,flipped);
toc

dispout = ['VTK output: ', filename];
disp(dispout)

% END %