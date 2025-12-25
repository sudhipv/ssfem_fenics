%% Then renumbering Points, Edges, Triangles and Bnodes into rPoints, rEdges, rTriangles and rBnodes
% Previously Rearranged the Nodes, Points as [Internal: Boundary]

%%-    Date:       Author:   Comments/Modifications:     
%   Mar/23/2013     AD       Original/Previously it was embedded with preprocessro2

%%-----------------------------------------------------------------------------------------------
% clear all; close all; clc;
clearvars
disp('Local 3D-Mesh re-arranging & re-numbering initiated...')

ndom = dlmread('../num_partition.dat', '');
points = dlmread('../points.dat','');
GlobalBN = dlmread('../boundary_nodes.dat', '');
nzs = 4;  % number of zeros to pad

%%--------------------------------------------------------------------------------------------------------
% ********** Re-arranging of NODES, POINTS and BNODES ************
%%--------------------------------------------------------------------------------------------------------
% disp('One-Level Mesh re-arranging Initiated...')

for i = 1:ndom 
    %tempNodes2  = ['../nodes00',num2str(i)];
    tempNodes2  = ['../nodes', pad(num2str(i),nzs,'left','0')];
    nodeTemp1 = dlmread([tempNodes2 '.dat']);
    nN = length(nodeTemp1);
    
    %tempBNodes2  = ['../bnodes00',num2str(i)];
    tempBNodes2  = ['../bnodes', pad(num2str(i),nzs,'left','0')];
    bnodeTemp1 = dlmread([tempBNodes2 '.dat']);
    nB = length(bnodeTemp1);
    
    for j=1:nB
        k=(nN-nB) + j;
        nodeTemp1(k) = bnodeTemp1(j);
    end
           
    %tempPoints  = ['../Nnodes00',num2str(i)];       %['Nnodes00',num2str(i)];
    tempPoints  = ['../Nnodes', pad(num2str(i),nzs,'left','0')];
    dlmwrite([tempPoints '.dat'], nodeTemp1, '\t');
    
    for j=1:nN
        pointsTemp2 = points(nodeTemp1(j),:);
        %nameP  = ['../rpoints00',num2str(i)];      %['npoints00',num2str(i)];
        nameP  = ['../rpoints', pad(num2str(i),nzs,'left','0')];
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
    %tempReN1  = ['../nodes00',num2str(i)];
    %tempReN1  = ['../Nnodes00',num2str(i)];
    tempReN1  = ['../Nnodes', pad(num2str(i),nzs,'left','0')];
    tempRead1 = dlmread([tempReN1 '.dat']);
    
    %tempMd2  = ['../meshdim00',num2str(i)]; 
    tempMd2  = ['../meshdim', pad(num2str(i),nzs,'left','0')];
    meshdim12 = dlmread([tempMd2 '.dat']);
    
    nP = meshdim12(1);                    % number of points
    nE = meshdim12(2);                    % number of edges
    nT = meshdim12(3);                    % number of triangles
    nH = meshdim12(4);                      % tetrahedrons
    nB = meshdim12(5);                    % number of global boundary nodes
       
   for j = 1:nP  %nPnt
    tempEg1 = tempRead1(j);
    node_renumET(tempEg1) = j;
   end
   
%%% edges:
%     %%%--------------------------------------------------------------------
%     %%% Exceptional handling 
%     tempTest2  = ['../edges00',num2str(i)];
%     test2 = load([tempTest2 '.dat']);
%     test3 = isempty(test2);
%     
%     if test3 == 1
%         % disp('floating subdomain');
%         NewEdges = [];
%     %%%--------------------------------------------------------------------    
%     else        
%         tempReEg2  = ['../edges00',num2str(i)];
%         tempRead2 = dlmread([tempReEg2 '.dat']);
%         
%         NewEdges = zeros(nE,7);
%         for j = 1:nE  %nEdg
%             tempE2 = tempRead2(j);
%             tempE3 = tempRead2(j:j,2);
%             node_renum1 =  node_renumET(tempE2);
%             node_renum2 =  node_renumET(tempE3);
%             NewEdges(j,:) = [node_renum1 node_renum2 tempRead2(j:j,3) tempRead2(j:j,4) tempRead2(j:j,5) tempRead2(j:j,6) tempRead2(j:j,7)];
%         end        
%     end
    NewEdges = []; %%Temporery emapy file untill we resolve shared edges issue 
    %newNameE  = ['../redges00',num2str(i)]; %['nedges00',num2str(i)];
    newNameE  = ['../redges', pad(num2str(i),nzs,'left','0')];
    dlmwrite([newNameE '.dat'], NewEdges);  
    
    %%% Triangles:
    %tempReT2  = ['../triangles00',num2str(i)];
    tempReT2  = ['../triangles', pad(num2str(i),nzs,'left','0')];
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
    %newNameT  = ['../rtriangles00',num2str(i)];  %['ntriangles00',num2str(i)];
    newNameT  = ['../rtriangles', pad(num2str(i),nzs,'left','0')];
    dlmwrite([newNameT '.dat'], NewTriangles); 
    
    %tempReT2  = ['../triangles00',num2str(i)];
    tempReT2  = ['../triangles', pad(num2str(i),nzs,'left','0')];
    tempRead3 = dlmread([tempReT2 '.dat']);
    
    
    %%% Tetrahedrons
    %tempReH2  = ['../mtetrahedrans00',num2str(i)];
    tempReH2  = ['../mtetrahedrans', pad(num2str(i),nzs,'left','0')];
    tempRead4 = dlmread([tempReH2 '.dat']);
    Newtetrahedrons = zeros(nT,4);
    Newmtetrahedrons = zeros(nT,5);
   for j = 1:nH  %nEdg
    tempH2 = tempRead4(j);
    tempH3 = tempRead4(j:j,2);
    tempH4 = tempRead4(j:j,3);
    tempH5 = tempRead4(j:j,4);
    node_renumH1 =  node_renumET(tempH2);
    node_renumH2 =  node_renumET(tempH3);
    node_renumH3 =  node_renumET(tempH4);
    node_renumH4 =  node_renumET(tempH5);
    Newtetrahedrons(j,:) = [node_renumH1  node_renumH2  node_renumH3  node_renumH4 ];
    Newmtetrahedrons(j,:) = [node_renumH1  node_renumH2  node_renumH3  node_renumH4 tempRead4(j:j,5)];
   end
    %newNameH  = ['../rtetrahedrons00',num2str(i)];%['ntriangles00',num2str(i)];
    newNameH  = ['../rtetrahedrons', pad(num2str(i),nzs,'left','0')];
    dlmwrite([newNameH '.dat'], Newtetrahedrons);    
    
    %newNameH2  = ['../rmtetrahedrons00',num2str(i)];
    newNameH2  = ['../rmtetrahedrons', pad(num2str(i),nzs,'left','0')];
    dlmwrite([newNameH2 '.dat'], Newmtetrahedrons);   
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

dlmwrite('../nboundary_nodes.dat', jj, '\t')

for i = 1:ndom
    %tempBN1  = ['../bnodes00',num2str(i)];% ['NBnodes00',num2str(i)];
    tempBN1  = ['../bnodes', pad(num2str(i),nzs,'left','0')];
    tempReadB1 = dlmread([tempBN1 '.dat']);
    nB = length(tempReadB1);
    NewNodes = zeros(nB,1);
   for j = 1:nB  %nPar
    tempB1 = tempReadB1(j);
    NewNodes(j,:) = node_renum(tempB1);
   end   
    %newNameLB  = ['../rbnodes00',num2str(i)];% ['nbnodes00',num2str(i)];
    newNameLB  = ['../rbnodes', pad(num2str(i),nzs,'left','0')];
    dlmwrite([newNameLB '.dat'], NewNodes);  
end

%%--------------------------------------------------------------------------------------------------------
% ********** Re-naming of 'nodes*'-->'Nnodes*' to make compatible with 'get'functions ************
% ********** Re-naming of 'points'-->'rpoints*' to make compatible with 'PostAssembling' ************
%%--------------------------------------------------------------------------------------------------------
%%% NOTE: This is necessory to make one-level and two-level decomposition 
%%% codes outputs compatible with each other and solver can work on both
for i=1:ndom
    %tempPoints  = ['../nodes00',num2str(i)];
    tempPoints  = ['../nodes', pad(num2str(i),nzs,'left','0')];
    pointsTemp1 = dlmread([tempPoints '.dat']);
    nN = length(pointsTemp1);
  
    %tempPoints2  = ['../Nnodes00',num2str(i)];
    tempPoints2  = ['../Nnodes', pad(num2str(i),nzs,'left','0')];
    dlmwrite([tempPoints2 '.dat'], pointsTemp1); 
    
    %tempBnodes  = ['../bnodes00',num2str(i)];
    tempBnodes  = ['../bnodes', pad(num2str(i),nzs,'left','0')];
    bnodesTemp1 = dlmread([tempBnodes '.dat']);
    
    %tempBnodes2  = ['../nbnodes00',num2str(i)];
    tempBnodes2  = ['../nbnodes', pad(num2str(i),nzs,'left','0')];
    dlmwrite([tempBnodes2 '.dat'], bnodesTemp1); 
    
    tempPoints  = ['../points', pad(num2str(i),nzs,'left','0')];
    pointsTemp1 = dlmread([tempPoints '.dat']);
    
    tempPoints2  = ['../rpoints', pad(num2str(i),nzs,'left','0')];
    dlmwrite([tempPoints2 '.dat'], pointsTemp1);
    
%     for j=1:nN
%         pointsTemp2 = points(pointsTemp1(j),:);
%         nameP  = ['../rpoints00',num2str(i)];
%         dlmwrite([nameP '.dat'], pointsTemp2, '-append');
%     end
end
%%--------------------------------------------------------------------------------------------------------
%%-********* END *************
%%--------------------------------------------------------------------------------------------------------
