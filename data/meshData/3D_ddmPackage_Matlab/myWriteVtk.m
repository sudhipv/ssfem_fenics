% test_vtk_write_binary
clearvars

% x = vtx_coord(:,1);
% y = vtx_coord(:,2);
% z = vtx_coord(:,3);

points = dlmread('../points.dat','');   %dlmread('points.dat','');
x = points(:,1);
y = points(:,2);
z = points(:,3);

tetra = dlmread('../tetrahedrons4c.dat','');   %dlmread('points.dat','');

% solution = dlmread('../../solution/final_solutionMatlab.dat','');

solution = dlmread('../../solution/final_solu.dat','');

% % make up a "velocity" field
% U = x.^2 + y.^3;
% V = cos(2*pi*(y + z));
% W = sin(2*pi*(x - y));

% setup filename
% filename = '../../vtkOutputs/matlab_poisson.vtk';  %% Binary (use paraview)

filename = '../../vtkOutputs/python_poisson.vtk';  


data_title = 'displacement';

% organize data
% .type field must be 'scalar' or 'vector'
% number of vector components must equal dimension of the mesh
data_struct.type = 'scalar';
data_struct.name = 'displacement';
data_struct.data = solution;

% data_struct(2).type = 'vector';
% data_struct(2).name = 'velocity';
% data_struct(2).data = [U,V,W];

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