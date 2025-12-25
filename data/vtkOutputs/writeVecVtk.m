% test_vtk_write_binary
clearvars

% x = vtx_coord(:,1);
% y = vtx_coord(:,2);
% z = vtx_coord(:,3);

points = dlmread('./points.dat','');   %dlmread('points.dat','');
x = points(:,1);
y = points(:,2);
z = points(:,3);
np = length(x);

tetra = dlmread('./tetrahedrons4c.dat','');   %dlmread('points.dat','');

solution = dlmread('mean_solutionVector.dat','');
U = solution(1:np);
V = solution(np+1:2*np);
W = solution(2*np+1:3*np);

% setup filename
filename = 'mean_solutionVector.vtk';  %% Binary (use paraview)
data_title = 'displacement';

% organize data
% .type field must be 'scalar' or 'vector'
% number of vector components must equal dimension of the mesh
% data_struct.type = 'scalar';
% data_struct.name = 'displacement';
% data_struct.data = solution;

data_struct.type = 'vector';
data_struct.name = 'displacements';
data_struct.data = [U,V,W];   %% Need to provide nodal components

% if *all data* is oriented column-wise, i.e. lots of rows, few columns,
% then set this to *false*
flipped = false;
% otherwise, if you want to transpose the data, then set this to *true*

% write the file
tic
stat = vtk_writeFn(filename,data_title,points,tetra,data_struct,flipped);
toc

solutionSD = dlmread('sd_solutionVector.dat','');
filenameSD = 'sd_solutionVector.vtk';  %% Binary (use paraview)
data_struct.data = solutionSD;

U = solutionSD(1:np);
V = solutionSD(np+1:2*np);
W = solutionSD(2*np+1:3*np);

% data_struct(2).type = 'vector';
% data_struct(2).name = 'velocity';
data_struct.data = [U,V,W];

% if *all data* is oriented column-wise, i.e. lots of rows, few columns,
% then set this to *false*
flipped = false;
% otherwise, if you want to transpose the data, then set this to *true*

% write the file
tic
stat = vtk_writeFn(filenameSD,data_title,points,tetra,data_struct,flipped);
toc

pcedata = dlmread('./pcedata.dat','');
nPCE = pcedata(3);

for i=2:nPCE

    readfile  = ['pce_', pad(num2str(i),2,'left','0')];
    readData = dlmread([readfile '.dat'],'');
    %data_struct.data = readData;
    writefile  = ['elasticity3D_',readfile,'.vtk'];

    U = readData(1:np);
    V = readData(np+1:2*np);
    W = readData(2*np+1:3*np);

    % data_struct(2).type = 'vector';
    % data_struct(2).name = 'velocity';
    data_struct.data = [U,V,W];

    % if *all data* is oriented column-wise, i.e. lots of rows, few columns,
    % then set this to *false*
    flipped = false;
    % otherwise, if you want to transpose the data, then set this to *true*

    % write the file
    stat = vtk_writeFn(writefile,data_title,points,tetra,data_struct,flipped);
    disp(writefile)
end


solutionSD = dlmread('mean_GlobSolutionVector.dat','');
filenameSD = 'mean_GlobSolutionVector.vtk';  %% Binary (use paraview)
data_struct.data = solutionSD;

U = solutionSD(1:np);
V = solutionSD(np+1:2*np);
W = solutionSD(2*np+1:3*np);

% data_struct(2).type = 'vector';
% data_struct(2).name = 'velocity';
data_struct.data = [U,V,W];

% if *all data* is oriented column-wise, i.e. lots of rows, few columns,
% then set this to *false*
flipped = false;
% otherwise, if you want to transpose the data, then set this to *true*

% write the file
tic
stat = vtk_writeFn(filenameSD,data_title,points,tetra,data_struct,flipped);
toc

dispout = ['VTK output: ', filename];
disp(dispout)

