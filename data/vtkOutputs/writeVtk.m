% test_vtk_write_binary
clearvars

points = dlmread('./points.dat','');   %dlmread('points.dat','');
x = points(:,1);
y = points(:,2);
z = points(:,3);
tetra = dlmread('./tetrahedrons4c.dat','');   %dlmread('points.dat','');


%%%% MEAN %%%%%%
solution = dlmread('mean_solutionVector.dat','');

% % make up a "velocity" field
% U = x.^2 + y.^3;
% V = cos(2*pi*(y + z));
% W = sin(2*pi*(x - y));

% setup filename
filename = 'poisson3D_mean.vtk';  %% Binary (use paraview)
data_title = 'u';

% organize data
% .type field must be 'scalar' or 'vector'
% number of vector components must equal dimension of the mesh
data_struct.type = 'scalar';
data_struct.name = 'u';
data_struct.data = solution;

% data_struct(2).type = 'vector';
% data_struct(2).name = 'velocity';
% data_struct(2).data = [U,V,W];

% if *all data* is oriented column-wise, i.e. lots of rows, few columns,
% then set this to *false*
flipped = false;
% otherwise, if you want to transpose the data, then set this to *true*

% write the file
stat = vtk_writeFn(filename,data_title,points,tetra,data_struct,flipped);
disp(filename)


%%% SD %%%% 
solutionSD = dlmread('sd_solutionVector.dat','');
filenameSD = 'poisson3D_sd.vtk';  %% Binary (use paraview)
data_struct.data = solutionSD;

% data_struct(2).type = 'vector';
% data_struct(2).name = 'velocity';
% data_struct(2).data = [U,V,W];

% if *all data* is oriented column-wise, i.e. lots of rows, few columns,
% then set this to *false*
flipped = false;
% otherwise, if you want to transpose the data, then set this to *true*

% write the file
stat = vtk_writeFn(filenameSD,data_title,points,tetra,data_struct,flipped);
disp(filenameSD)

pcedata = dlmread('./pcedata.dat','');
nPCE = pcedata(3);

for i=2:nPCE
  
    readfile  = ['pce_', pad(num2str(i),2,'left','0')];
    readData = dlmread([readfile '.dat'],'');
    data_struct.data = readData;
    writefile  = ['poisson3D_',readfile,'.vtk'];
    
    % data_struct(2).type = 'vector';
    % data_struct(2).name = 'velocity';
    % data_struct(2).data = [U,V,W];

    % if *all data* is oriented column-wise, i.e. lots of rows, few columns,
    % then set this to *false*
    flipped = false;
    % otherwise, if you want to transpose the data, then set this to *true*

    % write the file
    stat = vtk_writeFn(writefile,data_title,points,tetra,data_struct,flipped);
    disp(writefile)
end



% END %
