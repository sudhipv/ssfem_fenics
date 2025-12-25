
ndom = dlmread('../data/num_partition.dat', '');

for pn = 1:ndom 
    
    %%: Mesh-dimensions
    tempMd2  = ['../data/meshdim00',num2str(pn)]; % meshdim = [nPnt nEdg nTra nPar]
    meshdim12 = dlmread([tempMd2 '.dat']);
    
    nP = meshdim12(1);                    % number of points
    nE = meshdim12(2);                    % number of edges
    nT = meshdim12(3);                    % number of triangles
    nB = meshdim12(4);                    % number of global boundary nodes
    
    %%: FEniCS-Assembled Mat-Vecs
    tempAb  = ['../data/Ab',num2str(pn)];
    load(tempAb)
    
    [Aii, Aig, Agi, Agg, Fi, Fg] = oneLevelSchur(A, b, nP, nB);
    
    if (pn < 10)
        createD = horzcat('../data/FEniCS/Amats/subdom000',num2str(pn));
    else
        createD = horzcat('../data/FEniCS/Amats/subdom00',num2str(pn));
    end
    mkdir (createD)
    APii  = [createD,'/Aii'];
    dlmwrite([APii '.dat'], Aii, '\t');
    APgi  = [createD,'/Agi'];
    dlmwrite([APgi '.dat'], Agi, '\t');
    APgg  = [createD,'/Agg'];
    dlmwrite([APgg '.dat'], Agg, '\t');
    bFi  = [createD,'/Fi'];
    dlmwrite([bFi '.dat'], Fi, '\t');
    bFg = [createD,'/Fg'];
    dlmwrite([bFg '.dat'], Fg, '\t');
    
end