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
clearvars
disp('Local 3D-Mesh decomposition initiated...')

%%%%%%% LOADING ALL GLOBAL DATA %%%%%%%
num_partition = dlmread('../num_partition.dat', '');
if num_partition == 1
   disp('No Mesh Partitions: DDM fails, need atleast two partitions'); exit;   
end

points = dlmread('../points.dat','');
edges = dlmread('../edges.dat','');
triangles = dlmread('../triangles.dat','');
tetrahedrons = dlmread('../tetrahedrons.dat','');
mtetrahedrons = dlmread('../mtetrahedrons.dat','');
GlobalBN = dlmread('../boundary_nodes.dat', '');
meshDim = dlmread('../meshdim.dat', '');

%%% Floatin subdomains if triangles doen't exit
flotDom = sort(unique(triangles(:,4)));

%%% Global mesh dimensions:
np = meshDim(1);  %% potins
ne = meshDim(2);  %% edges
nt = meshDim(3);  %% triangles 
nh = meshDim(4);  %% tetrahedrons
nb = meshDim(5);  %% interface nodes 
ndom = meshDim(6);   %% number of subdomains
nb_sus = meshDim(7); %% SharedEdges

%%% Local mesh dimensions:
npi(ndom) = 0; %Points
nei(ndom) = 0; %edges
nti(ndom) = 0; %triangles
nhi(ndom) = 0; %tetra
nbi(ndom) = 0; %interface

%%%%%%%  BOUNDARY NODES %%%%%%%
for i = 1:nb
    bnodes(i) = GlobalBN(i);
end

nzs = 4;  % number of zeros to pad
tetrahedrons_on_boundary = dlmread('../tetrahedrons_on_boundary.dat', '');
for i = 1:nb_sus
cur_node = tetrahedrons_on_boundary(i,1);
cur_dom = tetrahedrons_on_boundary(i,2);
    for j = 1:nb
        if cur_node == bnodes(j)    
        %tempBnodes  = ['tempBnodes00',num2str(cur_dom)];
        tempBnodes  = ['tempBnodes',pad(num2str(cur_dom),nzs,'left','0')];
        dlmwrite([tempBnodes '.dat'], cur_node, '-append');
        nbi(cur_dom) = nbi(cur_dom)+1;
        end
    end
end

for j = 1:ndom
    %tempNodes  = ['tempBnodes00',num2str(j)];
    tempNodes  = ['tempBnodes',pad(num2str(j),nzs,'left','0')];
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
           nodes(nbicur) = cur_node;
        end
    end

  nbi(j) = nbicur;
 
  %bnodes2  = ['../bnodes00',num2str(j)];
  bnodes2  = ['../bnodes',pad(num2str(j),nzs,'left','0')];
  dlmwrite([bnodes2 '.dat'], nodes', '\t');
end
!rm tempBnodes*.dat

%%% Edge-boundary nodes (type of corner nodes)
edges_on_interface = [];
for j = 1:length(edges)
    ii = edges(j:j,1);
    jj = edges(j:j,2);
    for i = 1:length(GlobalBN)
        nn = GlobalBN(i);
        if ((ii == nn) || (jj == nn))
            inteface_edge = GlobalBN(i);
            res = ismember(inteface_edge,edges_on_interface);
            if res == 0
                edges_on_interface = [edges_on_interface; inteface_edge];
            end
        end
    end
end
dlmwrite('../edges_boundary_nodes.dat', edges_on_interface, '\t')

%%% Triangle-boundary nodes (type of corner nodes)
triangle_on_interface = [];
for j = 1:length(triangles)
    ii = triangles(j:j,1);
    jj = triangles(j:j,2);
    kk = triangles(j:j,3);
    for i = 1:length(GlobalBN)
        nn = GlobalBN(i);
        if ((ii == nn) || (jj == nn) || (kk == nn))
            inteface_triangle = GlobalBN(i);
            res = ismember(inteface_triangle,triangle_on_interface);
            if res == 0
                triangle_on_interface = [triangle_on_interface; inteface_triangle];
            end
        end
    end
end
dlmwrite('../triangles_boundary_nodes.dat', triangle_on_interface, '\t')

    
%%--------------------------------------------------------------------------------------------------------
% ****** TO SEPARATE 'NODES','EDGES','TRIANGLES'&'BOUNDARY_NODES'*******
% ************ FOR EACH PARTITION (LOCAL DATA) **************
%%--------------------------------------------------------------------------------------------------------

%%%%%%% EDGES & TRIANGLES %%%%%%%
FSD = [];    %  floating sub-domain
for k = 1:ndom
    
%%% Need toccheck the edges (physical) in .geo file to extract here
%    %%%--------------------------------------------------------------------
%    %%% Exception handling 
%    test1 = unique(edges(:,5));
%     if k ~= test1
%        %disp('floating domain'); %disp(k); 
%        domainNumber = k;
%        FSD = [FSD; k];
%        aaa = [];
%        floatingEdge  = ['../edges00',num2str(domainNumber)];
%        dlmwrite([floatingEdge '.dat'], aaa );  
%     end
%     %%%--------------------------------------------------------------------
%     
%     for i = 1:ne
%         j = edges(i:i,5:5);             
%         if j == k
%             tempEdg  = ['../edges00',num2str(j)];
%             dlmwrite([tempEdg '.dat'], edges(i,1:end), '-append');
%             nei(j) = nei(j) + 1;
%         % else
%         % disp('file not found')
%         end       
%     end
    %%Temporery empty edges (we don't need it FEniCS)
    %tempEmptyEdge = ['../edges00',num2str(k)];
    tempEmptyEdge = ['../edges', pad(num2str(k),nzs,'left','0')];
    tempEmptyEdg = [];
    dlmwrite([tempEmptyEdge '.dat'], tempEmptyEdg, '\t');
    
    for l = 1:nt
        m = triangles(l:l,end:end);
        ifsubexist = 0;
        if m == k
            ifsubexist = 1;
            %nameT  = ['../triangles00',num2str(m)];
            nameT = ['../triangles', pad(num2str(k),nzs,'left','0')];
            dlmwrite([nameT '.dat'], triangles(l,1:end), '-append');
            nti(m) = nti(m) + 1;
        end    
    end      
    
    
    for l = 1:nh
        h = mtetrahedrons(l:l,end:end);
        if h == k
            %nameH  = ['../tetrahedrans00',num2str(h)];
            nameH = ['../tetrahedrans', pad(num2str(h),nzs,'left','0')];
            dlmwrite([nameH '.dat'], mtetrahedrons(l,1:end-1), '-append');
            %nameH  = ['../mtetrahedrans00',num2str(h)];
            nameH = ['../mtetrahedrans', pad(num2str(h),nzs,'left','0')];
            dlmwrite([nameH '.dat'], mtetrahedrons(l,1:end), '-append');
            nhi(h) = nhi(h) + 1;
        end
    end  
    
    res = ismember(k,flotDom);
    if res == 0
        dispout = ['Floating Subdomain', int2str(k)];
        disp(dispout)
        tempEmptytriangles = ['../triangles', pad(num2str(k),nzs,'left','0')];
        Emptytriangles = [];
        dlmwrite([tempEmptytriangles '.dat'], Emptytriangles, '\t');
    end
    
    
end

%%%%%%% POINTS & NODES %%%%%%%
for i = 1:ndom
    %tempNodes2  = ['../bnodes00',num2str(i)];
    tempNodes2  = ['../bnodes' , pad(num2str(i),nzs,'left','0')];
    nodeTemp1 = dlmread([tempNodes2 '.dat']);
   for j = 1:nbi(i)
       nodes(j) = nodeTemp1(j);
   end

    %tempTriangles2  = ['../triangles00',num2str(i)];
    %triangle2 = dlmread([tempTriangles2 '.dat']);
    
    %tempTetrahedrons  = ['../mtetrahedrans00',num2str(i)];
    tempTetrahedrons =  ['../mtetrahedrans', pad(num2str(i),nzs,'left','0')];
    tetrahedrons2 = dlmread([tempTetrahedrons '.dat']);
   
    npicur = nbi(i);
    for j = 1:nhi(i)
        trianglei = tetrahedrons2(j,1:4);
        for l = 1:4
            found = 0;
            for k = 1:npicur
                if nodes(k) == trianglei(l)
                    found = 1;
                end
            end
            
            if found == 0
                npicur = npicur+1;
                nodes(npicur) = trianglei(l);
                %nodes4  = ['../nodes00',num2str(i)];
                nodes4 = ['../nodes', pad(num2str(i),nzs,'left','0')];
                dlmwrite([nodes4 '.dat'], trianglei(l)', '-append');
            end
        end
    end
    
    for j = 1:nbi(i)   
    %nodes4  = ['../nodes00',num2str(i)];
    nodes4 = ['../nodes', pad(num2str(i),nzs,'left','0')];
    dlmwrite([nodes4 '.dat'], nodes(j)', '-append');
    end
   npi(i) = npicur; 
end
                
for i=1:ndom
    %tempPoints  = ['../nodes00',num2str(i)];
    tempPoints  = ['../nodes', pad(num2str(i),nzs,'left','0')];
    pointsTemp1 = dlmread([tempPoints '.dat']);
    
    for j=1:npi(i)
        pointsTemp2 = points(pointsTemp1(j),:);
        %nameP  = ['../points00',num2str(i)];
        nameP  = ['../points', pad(num2str(i),nzs,'left','0')];
        dlmwrite([nameP '.dat'], pointsTemp2, '-append');
    end
end

%%%%%%% LOCAL MESH DIMENSION %%%%%%%

for i=1:ndom
    meshdimTemp = [npi(i), nei(i), nti(i), nhi(i) nbi(i)];
    %nameMeshDim  = ['../meshdim00',num2str(i)];
    nameMeshDim  = ['../meshdim', pad(num2str(i),nzs,'left','0')];
    dlmwrite([nameMeshDim '.dat'], meshdimTemp, '\t');
end

dlmwrite('../floatingDomians.dat', FSD, '\t');
nFSD = length(FSD);  % number of floating sub-domains

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