%% Then renumbering Points, Edges, Triangles and Bnodes into rPoints, rEdges, rTriangles and rBnodes
% Previously Rearranged the Nodes, Points as [Internal: Boundary]

%%-    Date:       Author:   Comments/Modifications:     
%   Mar/23/2013     AD       Original/Previously it was embedded with preprocessro2

%%-----------------------------------------------------------------------------------------------
% clear all; close all; clc;

ndom = dlmread('../num_partition.dat', '');
points = dlmread('../points.dat','');
GlobalBN = dlmread('../boundary_nodes.dat', '');

%%--------------------------------------------------------------------------------------------------------
% ********** Re-numbering of EDGES,TRIANGLES & BNODES ************
%%--------------------------------------------------------------------------------------------------------
% disp('Mesh re-numbering initiated...')

nzs = 4; 

for i = 1:ndom    
    node_renumET = [];
    %tempReN1  = ['../nodes00',num2str(i)];
    tempReN1  = ['../nodes', pad(num2str(i),nzs,'left','0')];
    tempRead1 = dlmread([tempReN1 '.dat']);
    
    %tempMd2  = ['../meshdim00',num2str(i)]; % meshdim = [nPnt nEdg nTra nPar]
    tempMd2  = ['../meshdim',pad(num2str(i),nzs,'left','0')];
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
    %tempTest2  = ['../edges00',num2str(i)];
    tempTest2  = ['../edges',pad(num2str(i),nzs,'left','0')];
    test2 = load([tempTest2 '.dat']);
    test3 = isempty(test2);
    
    if test3 == 1
        % disp('floating subdomain');
        NewEdges = [];
    %%%--------------------------------------------------------------------    
    else        
        %tempReEg2  = ['../edges00',num2str(i)];
        tempReEg2  = ['../edges',pad(num2str(i),nzs,'left','0')];
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
    %newNameE  = ['../redges00',num2str(i)]; %['nedges00',num2str(i)];
    newNameE  = ['../redges',pad(num2str(i),nzs,'left','0')];
    dlmwrite([newNameE '.dat'], NewEdges);  
        
    %tempReT2  = ['../triangles00',num2str(i)];
    tempReT2  = ['../triangles',pad(num2str(i),nzs,'left','0')];
    tempRead3 = dlmread([tempReT2 '.dat']);
    
    NewTriangles = zeros(nT,3);
   for j = 1:nT  %nEdg
    tempT2 = tempRead3(j);
    tempT3 = tempRead3(j:j,2);
    tempT4 = tempRead3(j:j,3);
    node_renumT1 =  node_renumET(tempT2);
    node_renumT2 =  node_renumET(tempT3);
    node_renumT3 =  node_renumET(tempT4);
    NewTriangles(j,:) = [node_renumT1  node_renumT2  node_renumT3];
   end
    %newNameT  = ['../rtriangles00',num2str(i)];%['ntriangles00',num2str(i)];
    %newNameT  = ['../rtriangles',pad(num2str(i),nzs,'left','0')];
    newNameT  = ['../triangles',pad(num2str(i),nzs,'left','0')];
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

dlmwrite('../nboundary_nodes.dat', jj, '\t')

for i = 1:ndom
    %tempBN1  = ['../bnodes00',num2str(i)];% ['NBnodes00',num2str(i)];
    tempBN1  = ['../bnodes',pad(num2str(i),nzs,'left','0')];
    tempReadB1 = dlmread([tempBN1 '.dat']);
    nB = length(tempReadB1);
    NewNodes = zeros(nB,1);
   for j = 1:nB  %nPar
    tempB1 = tempReadB1(j);
    NewNodes(j,:) = node_renum(tempB1);
   end   
    %newNameLB  = ['../rbnodes00',num2str(i)];% ['nbnodes00',num2str(i)];
    newNameLB  = ['../rbnodes',pad(num2str(i),nzs,'left','0')];
    dlmwrite([newNameLB '.dat'], NewNodes);  
end

%%--------------------------------------------------------------------------------------------------------
% ********** Re-naming of 'nodes*'-->'Nnodes*' to make compatible with 'get'functions ************
% ********** Re-naming of 'points'-->'rpoints*' to make compatible with 'PostAssembling' ************
%%--------------------------------------------------------------------------------------------------------
for i=1:ndom
    %tempPoints  = ['../nodes00',num2str(i)];
    tempPoints  = ['../nodes',pad(num2str(i),nzs,'left','0')];
    pointsTemp1 = dlmread([tempPoints '.dat']);
    nN = length(pointsTemp1);
  
    %tempPoints2  = ['../Nnodes00',num2str(i)];
    tempPoints2  = ['../Nnodes',pad(num2str(i),nzs,'left','0')];
    dlmwrite([tempPoints2 '.dat'], pointsTemp1); 
    
    for j=1:nN
        pointsTemp2 = points(pointsTemp1(j),:);
        %nameP  = ['../rpoints00',num2str(i)];
        %nameP  = ['../rpoints',pad(num2str(i),nzs,'left','0')];
        nameP  = ['../points',pad(num2str(i),nzs,'left','0')];
        dlmwrite([nameP '.dat'], pointsTemp2, '-append');
    end
end
%%--------------------------------------------------------------------------------------------------------
%%-********* END *************
%%--------------------------------------------------------------------------------------------------------
