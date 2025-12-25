%% EXTRACT THE CORNER NODES (a nodes belongs to 3 or more sub domain)
% Separate the Boundary(interface) nodes into corner (between 3 or more subdomains) & remaining nodes (between two-subdomains) 
% previously preprocessor3:  (nodes belongs to more than 2 sub domain) : Ajit/June/2013
% NOTE: Need to run globalDecomposer.m & localDecomposer.m before executing this script

%%-    Date:      Author:   Modifications/Comments:     
%   Jun/06/2013     AD      Original or Primary
%   Sup/08/2013     AD      Modified for corner/remaining nodes for each sub-domain & also to global corner nodes
%   Oct/02/2013     AD      Exceptional Handeling for reading empty edge00*.dat file
%   Oct/09/2013     AD      Position of Cnodes & Rnodes changed according to Bnodes position
%   Oct/23/2013     AD      Globalcorner_nodes.dat-->corner_nodes.dat:: corner_nodes.dat-->corner_nodesWithoutEdges.dat:: npoints, nedges, ntriangles-->rpoints, redges, rtriangles
%   Mar/19/2013     AD      Global remainig nodes extraction: remaining_nodes.dat  
%   Oct/23/2014     AD      Optimized


%%--------------------------------------------------------------------------------------------------------
% clear all; clc;
clearvars;
disp("Initiating corner/remaining node extraction")
load ../data/boundary_nodes.dat
ndom = dlmread('../data/num_partition.dat', '');
interface_edges = dlmread('../data/edges_boundary_nodes.dat','');
n = length(boundary_nodes);

%%--------------------------------------------------------------------------------------------------------
%%- Initiating matrix to assign count (to check in how many subdomain exist in each sub-domain.
aa = zeros(1,n);

for i = 1:ndom
    
    %nameN  = ['../data/bnodes00',num2str(i)];
    nameN  = ['../data/bnodes',pad(num2str(i),3,'left','0')];
    bnode1 = dlmread([nameN '.dat']);
    
    nb = length(bnode1);
    for j = 1:n
        for k = 1:nb
            if boundary_nodes(j) == bnode1(k)
                count = 1;
                aa(j) = aa(j) + count;
            end
        end        
    end
end

%%--------------------------------------------------------------------------------------------------------
%%- To Extracting the corner nodes: Node exist in more than two sub-domain
corner_nodes = [];
for i = 1:n
    cn = [];
    if aa(i) > 2
        cn = boundary_nodes(i);
    end
    corner_nodes = [corner_nodes; cn];
end
dlmwrite('../data/corner_nodesWithoutEdges.dat', corner_nodes, '\t')

%%------------------------------------------------------------------------------------
%%- Corner nodes(+plus the nodes at the ends of interface edges) for each sub-domain
GCnodes = [];
for i = 1:ndom
    
    %nameN2  = ['../data/bnodes00',num2str(i)];
    nameN2  = ['../data/bnodes',pad(num2str(i),3,'left','0')];
    bnode2 = dlmread([nameN2 '.dat']);   
    Cnodes = [];
    Indices = [];
    Rnodes = bnode2;
    
    %%- To Extract the nodes at the ends of interface edges
    nb = length(bnode2); ne = length(interface_edges); %nc = length(corner_nodes);
    for j = 1:nb
        for k = 1:ne
            if bnode2(j) == interface_edges(k)
                Cnodes2 = bnode2(j);
                Cnodes = [Cnodes; Cnodes2];
                GCnodes = [GCnodes; Cnodes2];
                Indices = [Indices; j];
            end
        end
    end
    
%     %%%----------------------------------------------------------------
%     %%% Exception Handling
%     tempTest2  = ['../data/edges00',num2str(i)];
%     test2 = load([tempTest2 '.dat']);
%     test3 = isempty(test2);
%     if test3 == 1
%         % disp('floating sub-domain')
%     else
%         %%%----------------------------------------------------------------
%         nameEd  = ['../data/edges00',num2str(i)];
%         edge1 = dlmread([nameEd '.dat']);
%         edge2 = [];
%         edge2 = [edge2; edge1(1:end,1)];
%         edge2 = [edge2; edge1(1:end,2)];
%         edge2 = unique(edge2);
%         
%         %%- To Extract the nodes at the ends of interface edges
%         nb = length(bnode2); ne = length(edge2); nc = length(corner_nodes);
%         for j = 1:nb
%             for k = 1:ne
%                 if bnode2(j) == edge2(k)
%                     Cnodes2 = bnode2(j);
%                     Cnodes = [Cnodes; Cnodes2];
%                     GCnodes = [GCnodes; Cnodes2];
%                     Indices = [Indices; j];
%                 end
%             end
%         end
%    end
    
    %%- To Extract the nodes shared between more than two sub-domains
    nb = length(bnode2); nc = length(corner_nodes); Indices2 = [];
    for j = 1:nb
        for l = 1:nc
            if bnode2(j) == corner_nodes(l)
                Cnodes3 = bnode2(j);
                Cnodes = [Cnodes; Cnodes3];
                GCnodes = [GCnodes; Cnodes3];
                Indices2 = [Indices2; j];
                Indices = [Indices; j];
            end
        end
    end
    
    %%- Corner nodes for each sub-domain (both, corner & edge)
    CnodesTWR = [];
    CnodesUni = unique(Cnodes);
    for k = 1:length(bnode2)
        for j = 1:length(CnodesUni)
            if bnode2(k) == CnodesUni(j)
                CnodesTWR = [CnodesTWR; CnodesUni(j)];
            end
        end
    end
    
    %CnodesN  = ['../data/cnodes00',num2str(i)];
    CnodesN  = ['../data/cnodes',pad(num2str(i),3,'left','0')];
    dlmwrite([CnodesN '.dat'], CnodesTWR, '\t');
    
    %%- Remaining Boundary nodes
%     %%%----------------------------------------------------------------
%     %%% Exception Handling
%     tempTest3  = ['../data/edges00',num2str(i)];
%     test4 = load([tempTest3 '.dat']);
%     test5 = isempty(test4);
%     if test5 == 1
%         Rnodes(Indices2) = [];
%         RnodesN  = ['../data/rnodes00',num2str(i)];
%         dlmwrite([RnodesN '.dat'], Rnodes, '\t');
%     else
%        %%%----------------------------------------------------------------
%         
%         Rnodes(Indices) = [];
%         RnodesTWR = [];
%         for k = 1:length(bnode2)
%             for j = 1:length(Rnodes)
%                 if bnode2(k) == Rnodes(j)
%                     RnodesTWR = [RnodesTWR; Rnodes(j)];
%                 end
%             end
%         end
%         
%         RnodesN  = ['../data/rnodes00',num2str(i)];
%         dlmwrite([RnodesN '.dat'], RnodesTWR, '\t');
%     end

    nb = length(bnode2); RnodesTWR=[];
    for k = 1:nb
        res = ismember(bnode2(k),CnodesTWR);
        if res == 0
            RnodesTWR = [RnodesTWR; bnode2(k)];
        end
    end
    
    %RnodesN  = ['../data/rnodes00',num2str(i)];
    RnodesN  = ['../data/rnodes',pad(num2str(i),3,'left','0')];
    dlmwrite([RnodesN '.dat'], RnodesTWR, '\t');

    
    CnodesNum = unique(Cnodes); numberCN = length(CnodesNum);
    RnodesNum = unique(RnodesTWR); numberRN = length(RnodesNum);
    numberCRN = [numberCN, numberRN];
    nameRCN  = ['../data/dimcrn',pad(num2str(i),3,'left','0')];
    dlmwrite([nameRCN '.dat'], numberCRN, '\t');
end

%%- Global, corner nodes + the nodes at the ends of interface edges
GCnodes = unique(GCnodes);
GlobalCnodes1 = vertcat(corner_nodes, GCnodes);
GlobalCnodes1 = unique(GlobalCnodes1);

GlobalCnodes = [];
for i = 1:length(boundary_nodes)
    for j = 1:length(GlobalCnodes1)
        if boundary_nodes(i) == GlobalCnodes1(j)
            GlobalCnodes = [GlobalCnodes; GlobalCnodes1(j)];
        end
    end
end
%GlobalCnodes;
dlmwrite('../data/corner_nodes.dat', GlobalCnodes, '\t')


%%- To Extract the Global Remaining Nodes 
GRnodes = setdiff(boundary_nodes, GlobalCnodes);

GlobalRnodes = [];
for i = 1:length(boundary_nodes)
    for j = 1:length(GRnodes)
        if boundary_nodes(i) == GRnodes(j)
            GlobalRnodes = [GlobalRnodes; GRnodes(j)];
        end
    end
end
%GlobalRnodes;
dlmwrite('../data/remaining_nodes.dat', GlobalRnodes, '\t')

%%--------------------------------------------------------------------------------------------------------
%%-********* END *************
%%--------------------------------------------------------------------------------------------------------
% corner_nodes = [];
% for i = 1:length(boundary_nodes)
%     for j = 1:length(corner_nodes2)
%     if boundary_nodes(i) == corner_nodes2(j)
%        corner_nodes = [corner_nodes; corner_nodes2(j)];
%     end
%     end
% end