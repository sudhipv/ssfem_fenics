%% Domain Decomposition Method (DDM Sub-structuring Methods) for Poisson's equation.       Date/Author: July/2013/Ajit
%%% Adapted Original Puffin Solver & MK's FORTEAN-DDM code        
%%% For details of Sub-structuring-DDM (Also called, Direct-DDM) refer to Waad's Thesis
%%% REFER WAAD's THESIS: ALGORITHM-2(Page-114)

%%-    Date:       Author:   Comments/Modifications:     
%   July/21/2013     AD       Original
%   Apr/02/2014      AD       getPET, getMehsDim & oneLevelSchur "Function" added

clearvars; clc;
disp('Direct-DDM Solution Initiated...')
tic
%%%%%%% Load the required mesh data %%%%%%%
GlobalBN = dlmread('../data/boundary_nodes.dat', '');    %% Global boundary nodes (Global interface)
ndom = dlmread('../data/num_partition.dat', '');         %% number of partitions

nBG = length(GlobalBN);                          %% number of global boundary nodes
RSiR = zeros(nBG,nBG);                           %% Initiating the require matrices
RGi = zeros(nBG,1);

for i = 1:ndom                                   %% ndom = number of sudomains    
    pn = i;                                      %% partition number
    
    [p, e, t] = getPET(pn);                      %% Mesh data: P E T
    [nP, nE, nT, nB] = getMeshDim(pn);           %% Mesh dimensions
    
    %%%%%%% Restriction matrix (Scatter & gather operator)
    Rmat = getRM(GlobalBN, nBG, nB, pn);        
    
    %%%%%%% Assemble matrix & Assemble vector %%%%%%%
    Amat = AssembleMatrix(p, e, t, 'PoissonModi', [], 0);
    bvec = AssembleVector(p, e, t, 'PoissonModi', [], 0);
    
    %%% All local stiffness matrices: notations are similar to Waad's Thesis:
    %%% (A_ii)^s =Aii, (A_itau)^s =Aig, (A_taui)^s =Agi,
    %%% (A_tautau)^s =Agg, (f_i)^s =Fi & (f_tau)^s =Fg
    
    [Aii, Aig, Agi, Agg, Fi, Fg] = oneLevelSchur(Amat, bvec, nP, nB);
    
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
dlmwrite('../data/globalInterface_sol.dat', Ub_g, '\t')
niter = 0;
dlmwrite('../data/num_itterations.dat', niter, '\t');

%%% Until here, script solves for global interface solution.
%%--------------------------------------------------------------------------------------------------------
%%************* END *****************
%%--------------------------------------------------------------------------------------------------------
