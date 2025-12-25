%% Domain Decomposition Method (DDM Sub-structuring Methods) 
%%- 3D Linear Elasticity Equation (Vector Valued PDE).
%
% Copyright (C) 2017 Ajit Desai, Ph.D. Candidate ajit.ndesai@gmail.com
%    
%%% For details of Sub-structuring-DDM (Also called, Direct-DDM) refer to Waad's Thesis
%%% REFER WAAD's THESIS: ALGORITHM-2(Page-114)

%%-    Date:       Author:   Comments/Modifications:     
%   July/21/2013     AD       Original
%   Apr/02/2014      AD       getPET, getMehsDim & oneLevelSchur "Function" added

clearvars; clc;
disp('Direct-DDM for 3D-Elasticity Initiated...')
tic
%%%%%%% Load the required mesh data %%%%%%%
GlobalBN = dlmread('../data/boundary_nodes.dat', '');    %% Global boundary nodes (Global interface)
ndom = dlmread('../data/num_partition.dat', '');         %% number of partitions

nVec = 3;                                        %% number of vector components  

nBG = length(GlobalBN);                          %% number of global boundary nodes
RSiR = zeros(nVec*nBG,nVec*nBG);                 %% Initiating the require matrices
RGi = zeros(nVec*nBG,1);

for i = 1:ndom                                   %% ndom = number of sudomains    
    pn = i;                                      %% partition number
    
    %[p, e, t, h] = getPETH(pn);                 %% Mesh data: P E T H
    [nP, nE, nT, nH, nB] = getMeshDim(pn);       %% Mesh dimensions
    nI = nP-nB;                                  %% Interior nodes
    
    %%%%%%% Restriction matrix (Scatter & gather operator)
    Rmat = zeros(nVec*nB,nVec*nBG);
    for j = 1:nVec     %% j=1:3 for 3D and j=1:2 for 2D
        Rmatc = getRM(GlobalBN, nBG, nB, pn);  
        Rmat( (j-1)*nB+1:j*nB, (j-1)*nBG+1:j*nBG ) = Rmatc(1:nB, 1:nBG);
    end       
    
    %%%%%%% Assemble matrix & Assemble vector %%%%%%%
    %%: FEniCS/dolfin-Assembled Mat-Vecs
    tempAb  = ['../outputs/Ab',num2str(pn)];
    load(tempAb)
    Amat = An;
    bvec = bn';
    
    %%: FEniCS/puffin-Assembled Mat-Vecs
    %Amat = AssembleMatrix(p, e, t, 'PoissonModi', [], 0);
    %bvec = AssembleVector(p, e, t, 'PoissonModi', [], 0);
    
    %%% All local stiffness matrices: notations are similar to Waad's Thesis:
    %%% (A_ii)^s =Aii, (A_itau)^s =Aig, (A_taui)^s =Agi,
    %%% (A_tautau)^s =Agg, (f_i)^s =Fi & (f_tau)^s =Fg
    [Aii, Aig, Agi, Agg, Fi, Fg] = oneLevelVecSchur(Amat, bvec, nP, nB, nVec);
    
    %%%%%%% DDM Direct Solver %%%%%%%
    % Thesis Form: S_s = (A_tautau)^s - (A_taui)^s * ((Aii)^s)^-1 * (A_itau)^s
    % out1 = ((Aii)^s)^-1 * (A_itau)^s :: Si = S_s
    out1 = (Aii)\Aig;
    Si   = Agg - (Agi*out1);
    
    %%% Thesis Form: (g_tau)^s = (f_tau)^s - (A_taui)^s * ((Aii)^s)^-1 * (f_i)^s
    % out2 = ((Aii)^s)^-1 * (f_i)^s :: Gi = (g_tau)^s
    out2 = (Aii)\Fi;
    Gi   = Fg - (Agi*out2);
     
    %%% Thesis Form: Rmat = R_s    % Out3 = (S_s)*(R_s)    % IMP: out3 = Si*Rmat;
    % Thesis Form: S = Sum_1_Ns{(R_s)'*(S_s)*(R_s)}
    % g_tau = Sum_1_Ns{(R_s)'*(g_tau)^s} :: RSiR =S, RGi =g_tau
    
    RSiRi = (Rmat)'*Si*Rmat;  %out3;
    RGii  = (Rmat)'*Gi;
    RSiR = RSiR + RSiRi;
    RGi = RGi + RGii;
end

%%% Thesis Form: S*(u_tau) = g_tau  :: u_tau =Ub_g
Ub_g = (RSiR)\RGi;                               %% Global boundary node solution vector
toc
dlmwrite('../data/dolfin_globalInterface_sol.dat', Ub_g, '\t')
niter = 0;
dlmwrite('../data/num_itterations.dat', niter, '\t');

disp('Next: run "PostAssemblingVec" ')

%%% Until here, script solves for global interface solution.
%%--------------------------------------------------------------------------------------------------------
%%************* END *****************
%%--------------------------------------------------------------------------------------------------------
