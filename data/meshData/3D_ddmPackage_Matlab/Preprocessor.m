%% Two-Level Mesh Decomposition
% This preprocessor can be used to decompose the mesh & create data for One-Level or Two-Level-DDM-PCGM Solvers. 
% The interface nodes of each sub-domain are partitioned into a set of 
% *Remaining nodes* (boundary nodes shared only between two adjacent sub-domains) and 
% *Corner nodes* (boundary nodes shared among more than two sub-domains plus the nodes at the ends of interface edges) 
% Interface Problem: "First Level" : Require Interface Nodes
% Coarse Problem   : "Second Level": Require Corner & Remaining Nodes as well 
% Note: Need to have the Mesh file, re-named "gmsh.msh" ready in the same folder before executing this script 

%%-    Date:       Author:   Comments/Modifications:     
%   Mar/23/2013     AD       Original
%   Oct/23/2014     AD       Optimized 
%%-----------------------------------------------------------------------------------------------

%% To clean previous data and clear Matlab workspace 
clear all; close all; clc;
disp('Cleaning previous data...')
!rm ../data/*.dat 

%% Input Mesh file
meshFile = sprintf('../gmsh.msh');

%% Mesh Decomposition into the Points(p), Edges(e), Triangles(t) and inteface-Nodes for the "whole" mesh
disp('Global Mesh decomposition initiated...')

globalDecomposer

disp('number of nodes'); disp(n);
disp('number of elements'); disp(m);
disp('number of partitions');disp(ndom);
disp('number of boundary nodes'); disp(nb);

clear all; 
%% Mesh Decomposition into the Points(p), Edges(e), Triangles(t) and inteface-Nodes for "each" sub-domain
disp('Local Mesh decomposition initiated...')

localDecomposer

%% Disply the number of floating subdomains  
disp('number of floating domains'); 
disp(nFSD); 
if nFSD < 7
    disp(FSD')
end

clear all;
%% Extract the Corner & Remaining Nodes and then Re-arranging 
disp('Corner-Remaining nodes extraction initiated...')

cornerRemaingNodes

clear all;
%% Re-arranging & Re-numbering of Points, Edges, Triangles and Boundary-Nodes
disp('Mesh re-arranging & re-numbering initiated...')

twoLevelRenum

disp('re-arranging as [I, R, C]');
disp('re-numbering as')
disp('points*    --> rpoints*');
disp('edges*     --> redges*');
disp('triangles* --> rtriangles*');
disp('bnodes*    --> rbnodes*');

% %% This is only require for the FETI-DP Solver
% booleanMatData

%%-----------------------------------------------------------------------------------------------