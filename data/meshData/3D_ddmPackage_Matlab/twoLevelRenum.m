%% Rearranging the Nodes, Points as [Internal: Remaining: Corner] = [I R C]
% Then renumbering Points, Edges, Triangles and Bnodes into rPoints, rEdges, rTriangles and rBnodes
% This we need to form the global assembly of stiffness matrix and the solution vector.

%%-    Date:       Author:   Comments/Modifications:     
%   Mar/23/2013     AD       Original/Previously it was embedded with the localDecomposer.m or the preprocessro2.m, 
%   Oct/23/2014     AD       Optimized

%%-----------------------------------------------------------------------------------------------
% clear all; close all; clc;

ndom = dlmread('../data/num_partition.dat', '');
points = dlmread('../data/points.dat','');
GlobalBN = dlmread('../data/boundary_nodes.dat', '');

%%--------------------------------------------------------------------------------------------------------
% ********** Re-arranging of NODES, POINTS and BNODES ************
%%--------------------------------------------------------------------------------------------------------
% disp('Mesh re-arranging initiated...')

for i = 1:ndom 
    tempNodes2  = ['../data/nodes00',num2str(i)];
    nodeTemp1 = dlmread([tempNodes2 '.dat']);
    nN = length(nodeTemp1);
    
    tempBNodes2  = ['../data/bnodes00',num2str(i)];
    bnodeTemp1 = dlmread([tempBNodes2 '.dat']);
    nB = length(bnodeTemp1);
    
    tempRNodes2  = ['../data/rnodes00',num2str(i)];
    rnodeTemp1 = dlmread([tempRNodes2 '.dat']);
    nR = length(rnodeTemp1);
    
    tempCNodes2  = ['../data/cnodes00',num2str(i)];
    cnodeTemp1 = dlmread([tempCNodes2 '.dat']);
    nC = length(cnodeTemp1);
    % nB = nR + nC;
    
    Bnodes = [rnodeTemp1 ; cnodeTemp1];
    
    bnodes2  = ['../data/nbnodes00',num2str(i)];       % ['NBnodes00',num2str(i)];
    dlmwrite([bnodes2 '.dat'], Bnodes, '\t');    
    
    for j=1:nR
        k=(nN-nB) + j;
        nodeTemp1(k) = rnodeTemp1(j);
    end
        
    for j=1:nC
        k=(nN-nC) + j;
        nodeTemp1(k) = cnodeTemp1(j);
    end
            
    tempPoints  = ['../data/Nnodes00',num2str(i)];       %['Nnodes00',num2str(i)];
    dlmwrite([tempPoints '.dat'], nodeTemp1, '\t');
    
    for j=1:nN
        pointsTemp2 = points(nodeTemp1(j),:);
        nameP  = ['../data/rpoints00',num2str(i)];      %['npoints00',num2str(i)];
        dlmwrite([nameP '.dat'], pointsTemp2, '-append');
    end        
end
clear nB

%%--------------------------------------------------------------------------------------------------------
% ********** Re-numbering of EDGES,TRIANGLES & BNODES ************
%%--------------------------------------------------------------------------------------------------------
% disp('Mesh re-numbering initiated...')

for i = 1:ndom
    node_renumET = [];
    tempReN1  = ['../data/Nnodes00',num2str(i)];
    tempRead1 = dlmread([tempReN1 '.dat']);
    
    tempMd2  = ['../data/meshdim00',num2str(i)]; % meshdim = [nPnt nEdg nTra nPar]
    meshdim12 = dlmread([tempMd2 '.dat']);
    
    nP = meshdim12(1);                    % number of points
    nE = meshdim12(2);                    % number of edges
    nT = meshdim12(3);                    % number of triangles
    nB = meshdim12(4);                    % number of global boundary nodes
    
    for j = 1:nP  %nPnt
        tempEg1 = tempRead1(j);
        node_renumET(tempEg1) = j;
    end
    
    %%%--------------------------------------------------------------------
    %%% Exceptional handling
    tempTest2  = ['../data/edges00',num2str(i)];
    test2 = load([tempTest2 '.dat']);
    test3 = isempty(test2);
    
    if test3 == 1
        % disp('floating sub-domain');
        NewEdges = [];
        %%%--------------------------------------------------------------------
    else
        tempReEg2  = ['../data/edges00',num2str(i)];
        tempRead2 = dlmread([tempReEg2 '.dat']);
        
        NewEdges = zeros(nE,7);
        for j = 1:nE  %nEdg
            tempE2 = tempRead2(j);
            tempE3 = tempRead2(j:j,2);
            node_renum1 =  node_renumET(tempE2);
            node_renum2 =  node_renumET(tempE3);
            NewEdges(j,:) = [node_renum1 node_renum2 tempRead2(j:j,3) tempRead2(j:j,4) tempRead2(j:j,5) tempRead2(j:j,6) tempRead2(j:j,7)];
        end
    end
    newNameE  = ['../data/redges00',num2str(i)]; %['nedges00',num2str(i)];
    dlmwrite([newNameE '.dat'], NewEdges);
    
    tempReT2  = ['../data/triangles00',num2str(i)];
    tempRead3 = dlmread([tempReT2 '.dat']);
    
    NewTriangles = zeros(nT,4);
    for j = 1:nT  %nEdg
        tempT2 = tempRead3(j);
        tempT3 = tempRead3(j:j,2);
        tempT4 = tempRead3(j:j,3);
        node_renumT1 =  node_renumET(tempT2);
        node_renumT2 =  node_renumET(tempT3);
        node_renumT3 =  node_renumET(tempT4);
        NewTriangles(j,:) = [node_renumT1  node_renumT2  node_renumT3 tempRead3(j:j,4)];
    end
    newNameT  = ['../data/rtriangles00',num2str(i)];%['ntriangles00',num2str(i)];
    dlmwrite([newNameT '.dat'], NewTriangles);
end

%%--------------------------------------------------------------------------------------------------------
% ********** Renumbering of Boundary Nodes ************
%%--------------------------------------------------------------------------------------------------------
node_renum = [];
nBG = length(GlobalBN);
jj = [];
for j = 1:nBG
    tempBG = GlobalBN(j);
    node_renum(tempBG) = j;
    jj = [jj; j];
end

dlmwrite('../data/nboundary_nodes.dat', jj, '\t')

for i = 1:ndom
    tempBN1  = ['../data/nbnodes00',num2str(i)];        % ['NBnodes00',num2str(i)];
    tempReadB1 = dlmread([tempBN1 '.dat']);
    nB = length(tempReadB1);
    NewNodes = zeros(nB,1);
    
    for j = 1:nB  
        tempB1 = tempReadB1(j);
        NewNodes(j,:) = node_renum(tempB1);
    end
    newNameLB  = ['../data/rbnodes00',num2str(i)];      % ['nbnodes00',num2str(i)];
    dlmwrite([newNameLB '.dat'], NewNodes);
end

%%--------------------------------------------------------------------------------------------------------
% ********** Display Mesh Information ************
%%--------------------------------------------------------------------------------------------------------
% NP = length(points);
% % nCG = load corner_nodes.dat
% % nFG = load floatingDomians.dat
% disp('number of nodes'); disp(NP);
% disp('number of partitions'); disp(ndom);
% disp('number of floating domians'); disp(nFG);
% disp('number of boundary nodes'); disp(nBG);
% disp('number of corner nodes'); disp(nCG);

%%--------------------------------------------------------------------------------------------------------
%%-********* END *************
%%--------------------------------------------------------------------------------------------------------
