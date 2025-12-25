%% Mesh Decomposition into Points, Elements, Triangles, Nodes and Boundary-Nodes for EACH Sub-domain
% Previously "preprocessor2": FOLLOWED BY preprocessor1.m TO EXTRACT LOCAL DATA:                 AJIT/March/2013
% NOTE: Need to run the 'globalDecomposer.m' before executing this script

%%-    Date:       Author:   Comments/Modifications:     
%   Mar/03/2013     AD       Original
%   Apr/23/2013     AD       Exception Handling for no-partitioned mesh
%   Oct/02/2013     AD       Exception Handling for reading empty "edge00*.dat" file is added
%   Oct/22/2013     AD       Node re-numbering part moved from this script to preprocessor4: 'Renum'
%   Feb/06/2014     AD       WORK: Need to modify this code for SPMD (data-parallel) Capability 
%   Oct/23/2014     AD       Optimized

%%-----------------------------------------------------------------------------------------------
% clear all; clc; 
% disp('Local Mesh decomposition initiated...')

%%%%%%% LOADING ALL GLOBAL DATA %%%%%%%
num_partition = dlmread('../data/num_partition.dat', '');
if num_partition == 1
   disp('No Mesh Partitions: DDM fails, need atleast two partitions'); exit;   
end

points = dlmread('../data/points.dat','');
edges = dlmread('../data/edges.dat','');
triangles = dlmread('../data/triangles.dat','');
GlobalBN = dlmread('../data/boundary_nodes.dat', '');
trianglesNotShrd = dlmread('../data/trianglesNotShrd.dat', '');
trianglesShrd = dlmread('../data/trianglesShrd.dat', '');
meshDim = dlmread('../data/meshdim.dat', '');

%%% Global mesh dimensions:
np = meshDim(1);
ne = meshDim(2);
nt= meshDim(3);
nb= meshDim(4); 
ndom = meshDim(5);
nb_sus = meshDim(6);

%%% Local mesh dimensions:
npi(ndom) = 0;
nbi(ndom) = 0;
nti(ndom) = 0;
nei(ndom) = 0;

%%%%%%%  BOUNDARY NODES %%%%%%%
for i = 1:nb
    bnodes(i) = GlobalBN(i);
end

triangles_on_boundary = dlmread('../data/triangles_on_boundary.dat', '');
for i = 1:nb_sus
cur_node = triangles_on_boundary(i,1);
cur_dom = triangles_on_boundary(i,2);
    for j = 1:nb
        if cur_node == bnodes(j)    
        tempBnodes  = ['tempBnodes00',num2str(cur_dom)];
        dlmwrite([tempBnodes '.dat'], cur_node, '-append');
        nbi(cur_dom) = nbi(cur_dom)+1;
        end
    end
end

for j = 1:ndom
    tempNodes  = ['tempBnodes00',num2str(j)];
    bnodeTemp1 = dlmread([tempNodes '.dat']);
    
    nodes = [];
    nbicur = 0;
 
    for i = 1:nbi(j)
        cur_node = bnodeTemp1(i);
        found = 0;
        for k = 1:nbicur
            if cur_node == nodes(k)
                found = 1;
            end
        end
          
        if found == 0
           nbicur = nbicur+1;
%            if nbicur > size(nodes)
%               disp('not enough allocation for nodes in this script')
%            exit
%            end
           nodes(nbicur) = cur_node;
        end
    end

  nbi(j) = nbicur;
 
  bnodes2  = ['../data/bnodes00',num2str(j)];
  dlmwrite([bnodes2 '.dat'], nodes', '\t');
end
!rm tempBnodes*.dat
    
%%--------------------------------------------------------------------------------------------------------
% ****** TO SEPARATE 'NODES','EDGES','TRIANGLES'&'BOUNDARY_NODES'*******
% ************ FOR EACH PARTITION (LOCAL DATA) **************
%%--------------------------------------------------------------------------------------------------------

%%%%%%% EDGES & TRIANGLES %%%%%%%
FSD = [];    %  floating sub-domain
for k = 1:ndom
    
   %%%--------------------------------------------------------------------
   %%% Exception handling 
   test1 = unique(edges(:,5));
    if k ~= test1
       %disp('floating domain'); %disp(k); 
       domainNumber = k;
       FSD = [FSD; k];
       aaa = [];
       floatingEdge  = ['../data/edges00',num2str(domainNumber)];
       dlmwrite([floatingEdge '.dat'], aaa );  
    end
    %%%--------------------------------------------------------------------
    
    for i = 1:ne
        j = edges(i:i,5:5);             
        if j == k
            tempEdg  = ['../data/edges00',num2str(j)];
            dlmwrite([tempEdg '.dat'], edges(i,1:end), '-append');
            nei(j) = nei(j) + 1;
        % else
        % disp('file not found')
        end       
    end
    
    for l = 1:nt
        m = triangles(l:l,end:end);
        if m == k
            nameT  = ['../data/triangles00',num2str(m)];
            dlmwrite([nameT '.dat'], triangles(l,1:end), '-append');
            nti(m) = nti(m) + 1;
        end
    end        
end

%%%%%%% POINTS & NODES %%%%%%%
for i = 1:ndom
    tempNodes2  = ['../data/bnodes00',num2str(i)];
    nodeTemp1 = dlmread([tempNodes2 '.dat']);
   for j = 1:nbi(i)
       nodes(j) = nodeTemp1(j);
   end

    tempTriangles2  = ['../data/triangles00',num2str(i)];
    triangle2 = dlmread([tempTriangles2 '.dat']);
   
    npicur = nbi(i);
    for j = 1:nti(i)
        trianglei = triangle2(j,1:3);
        for l = 1:3
            found = 0;
            for k = 1:npicur
                if nodes(k) == trianglei(l)
                    found = 1;
                end
            end
            
            if found == 0
                npicur = npicur+1;
                nodes(npicur) = trianglei(l);
                nodes4  = ['../data/nodes00',num2str(i)];
                dlmwrite([nodes4 '.dat'], trianglei(l)', '-append');
            end
        end
    end
    
    for j = 1:nbi(i)   
    nodes4  = ['../data/nodes00',num2str(i)];
    dlmwrite([nodes4 '.dat'], nodes(j)', '-append');
    end
   npi(i) = npicur; 
end
                
for i=1:ndom
    tempPoints  = ['../data/nodes00',num2str(i)];
    pointsTemp1 = dlmread([tempPoints '.dat']);
    
    for j=1:npi(i)
        pointsTemp2 = points(pointsTemp1(j),:);
        nameP  = ['../data/points00',num2str(i)];
        dlmwrite([nameP '.dat'], pointsTemp2, '-append');
    end
end

%%%%%%% LOCAL MESH DIMENSION %%%%%%%

for i=1:ndom
    meshdimTemp = [npi(i), nei(i), nti(i), nbi(i)];
    nameMeshDim  = ['../data/meshdim00',num2str(i)];
    dlmwrite([nameMeshDim '.dat'], meshdimTemp, '\t');
end

dlmwrite('../data/floatingDomians.dat', FSD, '\t');
nFSD = length(FSD);  % number of floating sub-domains

disp('Next: run "wireBasketNodes" ')


%%-NOTE: Renumbering Part has moved to new script "Renum"

%%--------------------------------------------------------------------------------------------------------
%%-********* END *************
%%--------------------------------------------------------------------------------------------------------

%  %%--------------------------------------------------------------------------------------------------------
%  % ********** Renumbering of EDGES and TRIANGLES************
%  %%--------------------------------------------------------------------------------------------------------
%  disp('Mesh re-numbering initiated...')
%  
%  for i = 1:ndom 
%     
%      node_renumET = [];
%      tempReN1  = ['nodes00',num2str(i)];
%      tempRead1 = dlmread([tempReN1 '.dat']);
%         
%     for j = 1:npi(i)  %nPnt
%      tempEg1 = tempRead1(j);
%      node_renumET(tempEg1) = j;
%     end
%     
%      %%%--------------------------------------------------------------------
%      %%% Exceptional handeling 
%      tempTest2  = ['edges00',num2str(i)];
%      test2 = load([tempTest2 '.dat']);
%      test3 = isempty(test2);
%      
%      if test3 == 1
%          disp('floating subdomain');
%          NewEdges = [];
%      %%%--------------------------------------------------------------------    
%      else        
%          tempReEg2  = ['edges00',num2str(i)];
%          tempRead2 = dlmread([tempReEg2 '.dat']);
%          
%          NewEdges = zeros(nei(i),7);
%          for j = 1:nei(i)  %nEdg
%              tempE2 = tempRead2(j);
%              tempE3 = tempRead2(j:j,2);
%              node_renum1 =  node_renumET(tempE2);
%              node_renum2 =  node_renumET(tempE3);
%              NewEdges(j,:) = [node_renum1 node_renum2 tempRead2(j:j,3) tempRead2(j:j,4) tempRead2(j:j,5) tempRead2(j:j,6) tempRead2(j:j,7)];
%          end        
%      end
%      newNameE  = ['nedges00',num2str(i)];
%      dlmwrite([newNameE '.dat'], NewEdges);  
%          
%      tempReT2  = ['triangles00',num2str(i)];
%      tempRead3 = dlmread([tempReT2 '.dat']);
%      
%      NewTriangles = zeros(nti(i),4);
%     for j = 1:nti(i)  %nEdg
%      tempT2 = tempRead3(j);
%      tempT3 = tempRead3(j:j,2);
%      tempT4 = tempRead3(j:j,3);
%      node_renumT1 =  node_renumET(tempT2);
%      node_renumT2 =  node_renumET(tempT3);
%      node_renumT3 =  node_renumET(tempT4);
%      NewTriangles(j,:) = [node_renumT1  node_renumT2  node_renumT3 tempRead3(j:j,4)];
%     end
%      newNameT  = ['ntriangles00',num2str(i)];
%      dlmwrite([newNameT '.dat'], NewTriangles);      
%  end
%  
%  %%--------------------------------------------------------------------------------------------------------
%  % ********** Renumbering of Boundary Nodes ************
%  %%--------------------------------------------------------------------------------------------------------
%  node_renum = [];
%  nBG = length(GlobalBN);
%  for j = 1:nBG
%      tempBG = GlobalBN(j);
%      node_renum(tempBG) = j;
%  end
%  
%  for i = 1:ndom
%      NewNodes = zeros(nbi(i),1);
%      tempBN1  = ['bnodes00',num2str(i)];
%      tempReadB1 = dlmread([tempBN1 '.dat']);
%     for j = 1:nbi(i)  %nPar
%      tempB1 = tempReadB1(j);
%      NewNodes(j,:) = node_renum(tempB1);
%     end   
%      newNameLB  = ['nbnodes00',num2str(i)];
%      dlmwrite([newNameLB '.dat'], NewNodes);  
%  end

%%--------------------------------------------------------------------------------------------------------
%%-********* END *************
%%--------------------------------------------------------------------------------------------------------

% % %%%%%%%%%%%%%%%%%%%% Old Method %%%%%%%%%%%%%%%%
% % %     node = [];
% % %     for l = 1:T
% % %         m = triangles(l:l,end:end);
% % %         node1 =[];
% % %         if m == k
% % %             m1 = k;
% % %             node1 = [node1; triangles(l,1:3)];
% % %         end
% % %         node = [node node1];
% % %         nodes = unique(node)';
% % %     end
% % %     
% % %     pointss = [];
% % %     for p = 1:length(nodes)
% % %         P = nodes(p);
% % %         points1 = points(P,:);
% % %         pointss = [pointss; points1];
% % %     end
% % %     
% % %     nameN  = ['nodes00',num2str(m1)];
% % %     dlmwrite([nameN '.dat'], nodes);
% % %     
% % %     bnode1 = dlmread([nameN '.dat']);
% % %     bnodes = intersect(GlobalBN, bnode1);
% % %     nameLB  = ['bnodes00',num2str(m1)];
% % %     dlmwrite([nameLB '.dat'], bnodes);
% % %     
% % %     nameP  = ['points00',num2str(m1)];
% % %     dlmwrite([nameP '.dat'], pointss, '\t');
% % %     
% % %     nPnt = length(pointss);
% % %     nPar = length(bnodes);
% % %     meshdim = [nPnt nEdg nTra nPar];
% % %     nameMeshDim  = ['meshdim00',num2str(m1)];
% % %     dlmwrite([nameMeshDim '.dat'], meshdim, '\t');
% % 
% % % %%******************************************************************************************
% % % %%******************************************************************************************

