%% PCGM-DDM (Preconditioned Conjugate Gradient Method) for Poisson's equation.       Date/Author: Aug/2013/Ajit
%%% Preconditioning of Schur Complement System
%%% Adapted Original Puffin Solver & MK's FORTRAN-PCGM-DDM code
%%% PCGM algorithm refer to Appendix-D.1 (page 174) of the MK's Thesis

%%-    Date:      Author:   Modification:
%   Aug/06/2013     AD      If-Else loop @line Number 130, in align with MK's thesis algorithm (previously it was aligned with MK's Fortran Code, another way of doing same)
%   Oct/05/2013     AD      Exceptional Handling with empty edge file is added
%   Oct/23/2013     AD      Globalcorner_nodes.dat --> corner_nodes.dat: npoints
%                            nedges, ntriangles --> rpoints, redges, rtriangles
%   Apr/02/2014     AD      getPET, getMehsDim & oneLevelSchur "Function" added

%%--------------------------------------------------------------------------------
%*********************************************************************************
%%%%%%% MESH DECOMPOSITION: PREPROCESSING: %%%%%%%
clear all; clc;
% Preprocessor
tic
%*********************************************************************************
%%%%%%% MAIN SOLVER: PROCESSING: %%%%%%%
clear all;
disp('Calculating The First Iteration Value (r_0)...')

%%%%%%% Load the required global data %%%%%%%
GlobalBN = dlmread('boundary_nodes.dat', '');    %% Global boundary nodes (Global interface)
ndom = dlmread('num_partition.dat', '');         %% number of partitions
nBG  = length(GlobalBN);                         %% number of global boundary nodes

%%%%%%% r_0 = G_tau - S*X_tau0: First itteration value %%%%%%%
%%% Since X_tau0 = 0, Initial Guess, r_0 = G_tau
rb_g = zeros(nBG,1);     %% RGi  = zeros(nBG,1);
for i = 1:ndom
    pn = i;                                      %% partition number
    
    [p, e, t] = getPET(pn);                      %% Mesh data: P E T
    [nP, nE, nT, nB] = getMeshDim(pn);           %% Mesh dimensions
    
%%% Assemble Matrix & Vector for each sub-domain
    Amat = AssembleMatrix(p, e, t, 'PoissonModi', [], 0);
    bvec = AssembleVector(p, e, t, 'PoissonModi', [], 0);
    
%%% Local stiffness matrices: oneLevel: notations are similar to Waad's Thesis:
    [Aii, Aig, Agi, Agg, Fi, Fg] = oneLevelSchur(Amat, bvec, nP, nB);
    
    out2 = (Aii)\Fi;
    Gi = Agi*out2;
    Gi = Fg - Gi;
    Rs = getRM(GlobalBN, nBG, nB, pn);                   %% getubg(Gi, m1, nB, nBG);
    RGi = Rs'*(Gi);
    rb_g = rb_g + RGi;
end

%%-------------------------------------------------------------------------
%%%%%%% PCGM Solver %%%%%%%
disp('PCGM Solution Initiated...')
Ub_g = zeros(nBG,1);       %% Zb_g = zeros(nBG,1);  % Qb_g = zeros(nBG,1);
tol = 10^(-5);             %%% Tolerance limit
maxiter = 100;             %%% Max Iterations

for j = 1:maxiter
    disp('Preconditioning of Ext. Schur Complement System')
    Zb_g = zeros(nBG,1);
    for i = 1:ndom
        pn = i;                                      %% partition number
        
        [p, e, t] = getPET(pn);                      %% Mesh data: P E T
        [nP, nE, nT, nB] = getMeshDim(pn);           %% Mesh dimensions
        
%%% Assemble Matrix & Vector for each sub-domain
        Amat = AssembleMatrix(p, e, t, 'PoissonModi', [], 0);
        bvec = AssembleVector(p, e, t, 'PoissonModi', [], 0);
        
        Agg  = Amat((nP-nB+1):nP,(nP-nB+1):nP);
        
        rb = getub(rb_g, pn, nB, nBG);    %% Zb = rb;
        Zb = (Agg)\rb;
        RZb = getubg(Zb, pn, nB, nBG);
        Zb_g = Zb_g + RZb;
    end
%%% Preconditioning of Extended Schur Complement System Ends Here.
    
    rho_next = dot(rb_g,Zb_g);
    if j == 1
        Pb_g = Zb_g;
    else
        beta = rho_next/rho_curr;
        Pb_g = Zb_g + (beta*Pb_g);
    end
    
%%-------------------------------------------------------------------------
    disp('Matrix-Vector Product for the Ext. Schur Complement System')
    Qb_g = zeros(nBG,1);
    for i = 1:ndom   % k --> np
        pn = i;                                      %% partition number
        
        [p, e, t] = getPET(pn);                      %% Mesh data: P E T
        [nP, nE, nT, nB] = getMeshDim(pn);           %% Mesh dimensions
        
%%% Assemble Matrix & Vector for each sub-domain
        Amat = AssembleMatrix(p, e, t, 'PoissonModi', [], 0);
        bvec = AssembleVector(p, e, t, 'PoissonModi', [], 0);
        
%%% Local stiffness matrices: oneLevel: notations are similar to Waad's Thesis:        
        [Aii, Aig, Agi, Agg] = oneLevelSchur(Amat, bvec, nP, nB);
        
        Pb   = getub(Pb_g, pn, nB, nBG);
        out2 = (Aig*Pb);
        Ui = (Aii)\out2;
        
        out2 = (Agi*Ui);
        Qb = (Agg*Pb);
        Qb = Qb - out2;
        
        RQb = getubg(Qb, pn, nB, nBG);
        Qb_g = Qb_g + RQb;
    end
%%% Matrix-Vector Product for the Ext. Schur Complement System Ends Here.
%%-------------------------------------------------------------------------
    
    rho_curr = rho_next;
    alpha = rho_curr/(dot(Qb_g,Pb_g));
    err = (alpha*alpha*(dot(Pb_g,Pb_g)))/(dot(Ub_g,Ub_g));
    Ub_g = Ub_g + (alpha*Pb_g);
    rb_g = rb_g - (alpha*Qb_g);
    
    disp('iteration'); disp(j); disp('relative error of '); disp(err);
    if err < tol
        break
    end
end
toc
dlmwrite('globalInterface_sol.dat', Ub_g, '\t');
dlmwrite('num_itterations.dat', j, '\t');

%t2 = toc;
%fprintf('Elapsed time is %0.2f seconds.\n', t2);

%%% Until here, the script solves for the global interface solution.
%*********************************************************************************

%%%%%% POST PROCESSING: RE-ASSEMBLING THE SOLUTION VECTOR: %%%%%%%
clear all;
% PostAssembling
%%--------------------------------------------------------------------------------------------------------
%%************* END *****************
%%--------------------------------------------------------------------------------------------------------

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%% Direct Solver %%%%%%%e
% % Thesis Form: S_s = (A_tautau)^s - (A_taui)^s * ((Aii)^s)^-1 * (A_itau)^s
% % Fortran code: call solve((np-nb)*npceout,nb*npceout,Aii,out1,Aig)
% % out1 = ((Aii)^s)^-1 * (A_itau)^s :: Si = S_s
%     out1 = inv(Aii)*Aig;
%     Si   = Agg - (Agi*out1);
%
% % Thesis Form: (g_tau)^s = (f_tau)^s - (A_taui)^s * ((Aii)^s)^-1 * (f_i)^s
% % Fortran code: call solve((np-nb)*npceout,1,Aii,out2,Fi)
% % out2 = ((Aii)^s)^-1 * (f_i)^s :: Gi = (g_tau)^s
%     out2 = inv(Aii)*Fi;
%     Gi   = Fg - (Agi*out2);
%
% % Thesis Form: Rmat = R_s
% % Out3 = (S_s)*(R_s)
%     %IMP: out3 = Si*Rmat;
%
% % Thesis Form: S = Sum_1_Ns{(R_s)'*(S_s)*(R_s)}
% % g_tau = Sum_1_Ns{(R_s)'*(g_tau)^s} :: RSiR =S, RGi =g_tau
% %%- FIXME: It should be summation, Check how it is done in fortran code
%
%     RSiRi = (Rmat)'*Si*Rmat;  %out3;
%     RGii  = (Rmat)'*Gi;
%     RSiR = RSiR + RSiRi;
%     RGi = RGi + RGii;
% end
%
% % Thesis Form: S*(u_tau) = g_tau  :: u_tau =Ub_g
% % Fortran Code: call solve(nbg*npceout,1,S,Ub_g,calG)
%     Ub_g = inv(RSiR)*RGi;
%     dlmwrite('globalInterface_sol.dat', Ub_g, '\t')
%
% % Until here, script solves for global interface solution.

% %%%%%%% RESTRICTION MATRIX (SCATTER & GATHER OPERATOR) %%%%%%%
%     Rmat = zeros(length(bnode12), length(GlobalBN));
%
%     for i = 1:length(bnode12)
%         Ri = bnode12(i);
%         for j = 1:length(GlobalBN)
%             Rj = GlobalBN(j);
%             if Ri == Rj
%             Rmat(i,j) = 1;
%             else
%             Rmat(i,j) = 0;
%             end
%         end
%     end

% %%%%%%% OLD TRY TO FIND LOCAL & GLOBAL SOLUTION %%%%%%%
%  for i = 1:ndom
%      pn = i;
%  Ub = getub(Ub_g, pn, nB, nBG);
%  out2 = (Aig*Ub);
%
%  Ui = Fi-out2;
%  Ui = inv(Aii)*Ui;
%
%  tempReadUi  = ['Ui00',num2str(pn)];
%  dlmwrite([tempReadUi '.dat'], Ui, '\t');
%
%  tempReadUb_g  = ['Ub_g00',num2str(pn)];
%  dlmwrite([tempReadUb_g '.dat'], Ub_g, '\t');
%
%  Ub_gFinal = Ub_gFinal + Ub_g;
% % end
%%--------------------------------------------------------------------------------------------------------
%%************* END *****************
%%--------------------------------------------------------------------------------------------------------
