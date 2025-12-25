%
% Copyright (C) 2017 Ajit Desai, Ph.D. Candidate ajit.ndesai@gmail.com
%
%% Decompose the input mesh file into Global Points, Elements, Triangles and Boundary (Interface) nodes
% Previously 'preprocessor1': SCRIPT TO EXTRACT POINTS, ELEMENTS, AND TRIANGLES (P, E, T) FROM 2D".msh" (GMSH) FILE : Ajit/March/2013
% ALL BOUNDARY ELEMENTS (GLOBAL): Global Nodes, GBNodes etc
% NOTE: You need to have the Mesh file "gmsh.msh" ready before executing this script 
% gmsh.msh file need to be created using the GMSH version-2.7.1 and above.  

%%-    Date:       Author:   Comments/Modifications:     
%   Mar/03/2013     AD       Original
%   Apr/23/2013     AD       Exception Handling for no-partitions mesh
%   Oct/23/2013     AD       Modified to read input mesh file in 'main processor' 
%   Feb/06/2014     AD       The method to read "mesh-file(.msh)" need to be modify, this method is very slow for large scale mesh
%   Oct/23/2014     AD       Optimized 

%%-----------------------------------------------------------------------------------------------
clearvars   %% here add line to remove old data files
disp('Global 3D-Mesh decomposition initiated...')

%%%%%%% GMESH/MESH FILE INPUT %%%%%%%
meshFile = sprintf('../foo3D.msh');

R1 = 4; C1 = 0; R2 = 4; C2 = 0;
range1 = [R1 C1 R2 C2];
n = dlmread(meshFile,'',range1);
number_of_nodes = n

%%%%%%% TO EXTRACTS THE POINTS "p" %%%%%%%
R1 = 4+1; C1 = 1; R2 = n+4; C2 = 3;  % from C1 to C3
range2 = [R1 C1 R2 C2]; 
points = dlmread(meshFile,'',range2);
xlims(1:1,1) = min(points(:,1));
xlims(1:1,2) = max(points(:,1));
ylims(1:1,1) = min(points(:,2));
ylims(1:1,2) = max(points(:,2));
zlims(1:1,1) = min(points(:,3));
zlims(1:1,2) = max(points(:,3));
dbounds = vertcat(xlims,ylims,zlims);

%%%%%%% TO EXTRACTS THE CELLS (EDGES "e" TRIANGLE "t" & tetrahedrons) %%%%%
fid=fopen(meshFile,'r');                                 
nhead = n+7;                              
for i=1:nhead,  buffer = fgetl(fid);  end
for i = 1:1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    temp=sscanf(tline,'%f');
    m = temp';
    number_of_elements = m;
end

shrEdg = []; edges = []; triangles = []; 
tetrahedronsShrd = []; tetrahedronsNotShrd = [];
tetrahedrons = []; mtetrahedrons = [];
for i = 1:m
    sE = []; bE = []; bT = []; aa = []; sT =[]; nsT =[]; bH=[]; bH2=[];
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    temp=sscanf(tline,'%f');
    aa = temp';
    
%%% to find the number of partitions %%%
    if i == 1
        naa = length(aa);
    end
    
    if aa(3) > 4 %% Shared vertices & tetrahedons
        sE = [aa(end-3) aa(7); aa(end-2) aa(7); aa(end-1) aa(7); aa(end) aa(7)];
        sT = [aa(end-3), aa(end-2) aa(end-1) aa(end) aa(7)];
    %elseif aa(2)==2 && aa(3) == 4
    %    nsT = [aa(end-3), aa(end-2) aa(end-1) aa(end) aa(7)];
    end
    if aa(2) == 1      %% edges
        bE = [aa(end-1) aa(end) aa(end-1)];
    elseif aa(2) == 2  %% triangles
        bT = [aa(end-2) aa(end-1) aa(end) aa(end-3)];
    elseif aa(2) == 4  %% tetrahedrons 
        bH = [aa(end-3), aa(end-2) aa(end-1) aa(end)];
        bH2 = [aa(end-3), aa(end-2) aa(end-1) aa(end) aa(7)];
    end
    shrEdg = [shrEdg; sE];                                    % Nodes of Shared Triangles
    edges  = [edges; bE];
    triangles = [triangles; bT];
    tetrahedrons = [tetrahedrons; bH];
    mtetrahedrons = [mtetrahedrons; bH2];
    tetrahedronsShrd = [tetrahedronsShrd; sT];                      % Shared Triangles
    tetrahedronsNotShrd = [tetrahedronsNotShrd; nsT];
    clear aa  ans buffer nhead temp tline sE bE sT nsE sT bH;
end
fclose(fid);

number_of_edges = length(edges(:,end));
number_of_triangles = length(triangles(:,end));
number_of_tetrahydrons = length(tetrahedrons(:,end));

%%%%%%% TO EXTRACTS THE GLOBAL BOUNDARY NODES %%%%%%%
Ne = length(shrEdg);
nodes = zeros(Ne,3);
nb = 0; j = 0;
for i = 1:Ne
    cur_node = shrEdg(i,1);
    cur_dom = shrEdg(i,2);
    found = 0;
    for k = 1:j
        if cur_node == nodes(k,1)
            found = 1;
            if cur_dom ~=nodes(k,2)
                nb = nb+1;
                nodes(k,3) = 1;
            end
        end
    end
    if found == 0
        j = j+1;
        nodes(j,1) = cur_node;
        nodes(j,2) = cur_dom;
    end
end

boundary_nodes = [];
nb = 0;
for i = 1:j
    if nodes(i,3) == 1
        nb = nb + 1;
        boundary_nodes = [boundary_nodes ; nodes(i,1)];
    end
end

%%% Exception Handling 
num_partition = max(tetrahedronsShrd(:,5));

if num_partition == 1
   disp('No Partitions in the MESH');
   num_partition = 1;   
else
   dispP = ['number of partitions = ', num2str(num_partition)];
   disp(dispP)
end
ndom = num_partition;
nb_sus = Ne;
meshdim = [number_of_nodes, number_of_edges, number_of_triangles, nb, ndom, nb_sus];
meshdim3D = [number_of_nodes, number_of_edges, number_of_triangles, number_of_tetrahydrons, nb, ndom, nb_sus];

%%%%%%% EXPORT DATA *.dat FILES %%%%%%%
dlmwrite('../boundary_nodes.dat', boundary_nodes, '\t')
dlmwrite('../points.dat', points, '\t')
dlmwrite('../edges.dat', edges, '\t')
dlmwrite('../triangles.dat', triangles, '\t')
dlmwrite('../tetrahedrons.dat', tetrahedrons, '\t')
dlmwrite('../mtetrahedrons.dat', mtetrahedrons, '\t')
dlmwrite('../tetrahedronsShrd.dat', tetrahedronsShrd, '\t')
dlmwrite('../tetrahedrons_on_boundary.dat', shrEdg, '\t')
dlmwrite('../num_partition.dat', num_partition, '\t')
dlmwrite('../dbounds.dat', dbounds, '\t')
dlmwrite('../meshdim.dat', meshdim3D, '\t')  %% both names are same??

%%%%%%% PLOTTING MESH %%%%%%%
% p = points';
% e = edges';
% t = triangles';
% 
% pdemesh(p,e,t)
% axis equal
% title('Mesh File')

%%-IMP NOTES:
% save boundary_nodes.dat boundary_nodes -ascii; 
% '/t' for for tab deliminator 
% % save('boundary_nodes.mat', 'boundary_nodes')

%%--------------------------------------------------------------------------------------------------------
%%-********* END *************
%%--------------------------------------------------------------------------------------------------------
