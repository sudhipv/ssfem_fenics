%
% Copyright (C) 2017 Ajit Desai, Ph.D. Candidate ajit.ndesai@gmail.com
%
% EXTRACT THE CORNER NODES (a nodes belongs to 3 or more sub domain + physical edges
% EXTRACT THE WIRE-BASKET NODES (Corner + Interface Edges)
% Also extracts remaining nodes (after corner or after wire-basket)
% NOTE: Need to run globalDecomposer3D.m & localDecomposer3D.m before executing this script
%%--------------------------------------------------------------------------------------------------------

clearvars;
disp("Initiating wire-basek (corener+edges) & remaining (faces) nodes extraction")
load ../boundary_nodes.dat
ndom = dlmread('../num_partition.dat', '');
interface_edges = dlmread('../edges_boundary_nodes.dat','');
interface_triangles = dlmread('../triangles_boundary_nodes.dat','');
n = length(boundary_nodes);

npad = 4; 

%%--------------------------------------------------------------------------------------------------------
decomposition = 1;   %% For wire basket select decomposition = 1
%%--------------------------------------------------------------------------------------------------------
%%--------------------------------------------------------------------------------------------------------
%%- Initiating matrix to assign count (to check in how many subdomain exist in each sub-domain.
aa = zeros(1,n);

for i = 1:ndom
    
    nameN  = ['../bnodes',pad(num2str(i),npad,'left','0')];
    bnodes = dlmread([nameN '.dat']);
    
    nb = length(bnodes);
    for j = 1:n
        for k = 1:nb
            if boundary_nodes(j) == bnodes(k)
                count = 1;
                aa(j) = aa(j) + count;
            end
        end        
    end
end

%%--------------------------------------------------------------------------------------------------------
if decomposition == 1 
    disp('Using Wire-Basket')
%%--------------------------------------------------------------------------------------------------------
    %%- To Extracting the wire basek (C+E) nodes: Node exist in more than two sub-domain
    wire_basket = [];
    remaining_nodes = [];
    for i = 1:n
        cn = [];
        rn = [];
        if aa(i) > 2
            cn = boundary_nodes(i);
        else 
            res = ismember(boundary_nodes(i),interface_triangles);
            if res == 1
                cn = boundary_nodes(i);
            else
                rn = boundary_nodes(i);
            end
        end
        wire_basket = [wire_basket; cn];
        remaining_nodes = [remaining_nodes; rn];
    end
    dlmwrite('../corner_nodes.dat', wire_basket, '\t')
    dlmwrite('../remaining_nodes.dat', remaining_nodes, '\t')

    %%--------------------------------------------------------------------------------------------------------
    %%- Subdomain level wire-basket and remainig nodes 

    for i = 1:ndom

        nameN  = ['../bnodes',pad(num2str(i),npad,'left','0')];
        bnodes = dlmread([nameN '.dat']);
        nb = length(bnodes);

        local_wb = [];
        local_rn = [];

        for j = 1:nb
            cn =[];
            rn =[];
            res = ismember(bnodes(j),wire_basket);
            if res == 1
                cn = bnodes(j);
            else
                rn = bnodes(j);
            end
            local_wb = [local_wb; cn];
            local_rn = [local_rn; rn];
        end        

        wbnodes  = ['../cnodes',pad(num2str(i),npad,'left','0')];
        dlmwrite([wbnodes '.dat'], local_wb, '\t');

        rnodes  = ['../rnodes',pad(num2str(i),npad,'left','0')];
        dlmwrite([rnodes '.dat'], local_rn, '\t');

        numberWB = length(unique(local_wb));
        numberRN = length(unique(local_rn));
        numberCRN = [numberWB, numberRN];
        nameRCN  = ['../dimcrn',pad(num2str(i),npad,'left','0')];
        dlmwrite([nameRCN '.dat'], numberCRN, '\t');

    end
    
%%--------------------------------------------------------------------------------------------------------    
else 
    disp('Using Only-Corener')
%%--------------------------------------------------------------------------------------------------------
    %%- To Extracting the corner nodes: Node exist in more than two sub-domain
    corner_nodes = [];
    remaining_nodes = [];
    for i = 1:n
        cn = [];
        rn = [];
        if aa(i) > 2
            cn = boundary_nodes(i);
        else 
            res = ismember(boundary_nodes(i),interface_edges);
            if res == 1
                cn = boundary_nodes(i);
            else
                rn = boundary_nodes(i);
            end
        end
        corner_nodes = [corner_nodes; cn];
        remaining_nodes = [remaining_nodes; rn];
    end
    dlmwrite('../corner_nodes.dat', corner_nodes, '\t')
    dlmwrite('../remaining_nodes.dat', remaining_nodes, '\t')

    %%--------------------------------------------------------------------------------------------------------
    %%- Subdomain level corner and remainig nodes 
    for i = 1:ndom

        nameN  = ['../bnodes',pad(num2str(i),npad,'left','0')];
        bnodes = dlmread([nameN '.dat']);
        nb = length(bnodes);

        local_cn = [];
        local_rn = [];

        for j = 1:nb
            cn =[];
            rn =[];
            res = ismember(bnodes(j),corner_nodes);
            if res == 1
                cn = bnodes(j);
            else
                rn = bnodes(j);
            end
            local_cn = [local_cn; cn];
            local_rn = [local_rn; rn];
        end        

        cnodes  = ['../cnodes',pad(num2str(i),npad,'left','0')];
        dlmwrite([cnodes '.dat'], local_cn, '\t');

        rnodes  = ['../rnodes',pad(num2str(i),npad,'left','0')];
        dlmwrite([rnodes '.dat'], local_rn, '\t');


        numberCN = length(unique(local_cn));
        numberRN = length(unique(local_rn));
        numberCRN = [numberCN, numberRN];
        nameRCN  = ['../dimcrn',pad(num2str(i),npad,'left','0')];
        dlmwrite([nameRCN '.dat'], numberCRN, '\t');

    end
end



%%--------------------------------------------------------------------------------------------------------
%%-********* END *************
%%--------------------------------------------------------------------------------------------------------