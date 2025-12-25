%% PCGM-DDM (Preconditioned Conjugate Gradient Method) for Poisson's equation.       Date/Author: Aug/2013/Ajit
%%% Two-Level Neumann-Neumann Preconditioning of Schur Complement System
%%% For detail refer Algorithm-12 of WS_Thesis (Pg-146)
%%% Adapted Original Puffin Solver & MK's FORTRAN-PCGM-DDM code
%%% For details for the PCGM algorithm refer to Appendix-D or the MK's Thesis

%%-    Date:      Author:   Modification:
%   Aug/06/2013     AD      If-Else loop @line Number 130, in align with MK's thesis algorithm (previously it was aligned with MK's Fortran Code, another way of doing same)
%   Sup/09/2013     AD      Aic, Air, Arr, Acc etc & GlobalCornerNodes, nC, nR etc are added in PNNCGM solver to incorporate TwoNN-PCGM solver, As per waad's thesis Pg.(170).
%   Sup/13/2013     AD      PCGM with lumped preconditioner is added for step-9 of Algorithm 12 or WS_Thesis
%   Oct/02/2013     AD      Fcc*Zc = dc solved directly instead of lumped PCGM solver/ getRM is used to extract Rmat instead of getub & getubg
%   Oct/05/2013     AD      Exception Handling with the empty edge file is added
%   Oct/23/2013     AD      Globalcorner_nodes.dat --> corner_nodes.dat: npoints, nedges, ntriangles --> rpoints, redges, rtriangles (All Corrected)
%   Feb/19/2013     AD      Added Lumped-PCGM solver for Corner-Schur-Compliment Problem
%   Apr/02/2014     AD      getPET, getMehsDim, oneLevelSchur & twoLevelSchur "Function" added

%%--------------------------------------------------------------------------------
%*********************************************************************************
%%%%%%% MESH DECOMPOSITION: PREPROCESSING: %%%%%%%
clear all; clc;
% twoLevelPreprocessor

%*********************************************************************************
%%%%%%% MAIN SOLVER: PROCESSING: %%%%%%%
clear all;
disp('Calculating The First Iteration Value (r_0)...')

%%%%%%% Load the required global data %%%%%%%
GlobalBN = dlmread('../boundary_nodes.dat', '');        %% Global boundary nodes (Global interface)
GlobalCN = dlmread('../corner_nodes.dat', '');          %% Global corner nodes
ndom = dlmread('../num_partition.dat', '');             %% number of partitions
nBG  = length(GlobalBN);
nBGc = length(GlobalCN);

%%%%%%% r_0 = G_tau - S*X_tau0 %%%%%%%
%%% Since X_tau0 = 0, Initial Guess, r_0 = G_tau
rb_g = zeros(nBG,1);   % RGi  = zeros(nBG,1);

for i = 1:ndom
    pn = i;                                      %% partition number
    
    [p, e, t] = getPET(pn);                      %% Mesh data: P E T
    [nP, nE, nT, nB] = getMeshDim(pn);           %% Mesh dimensions
    
    %%%%%%% Assemble matrix %%%%%%%
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

%*********************************************************************************
%%%%%%% PCGM Solver %%%%%%%
disp('PTwo-NNCGM Solution Initiated...')
Ub_g = zeros(nBG,1);      % Zb_g = zeros(nBG,1); % Qb_g = zeros(nBG,1);

tol = 10^(-5);                             %% Tolerance limit
maxiter = 100;                             %% Max Iteration
for j = 1:maxiter
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Two-Level NN Preconditioning of Ext. Schur Complement System')
    dc  = zeros(nBGc,1);
    Fcc = zeros(nBGc,nBGc);
    Mpre = zeros(nBGc,nBGc);
    %%% Zb_g = zeros(nBG,1); %%% RZb  = zeros(nBG,1);
    for i = 1:ndom
        pn = i;                                      %% partition number
        
        [p, e, t] = getPET(pn);                      %% Mesh data: P E T
        [nP, nE, nT, nB] = getMeshDim(pn);           %% Mesh dimensions
        
        [nC, nR] = getCorRemDim(pn);                 %% Mesh dim: corner & remaining
        
%%% Assemble Matrix & Vector for each sub-domain
        Amat = AssembleMatrix(p, e, t, 'PoissonModi', [], 0);
        bvec = AssembleVector(p, e, t, 'PoissonModi', [], 0);
        
%%% Local stiffness matrices: oneLevel: notations are similar to Waad's Thesis:
        [Aii] = oneLevelSchur(Amat, bvec, nP, nB);
        [Air, Aic, Ari, Aci, Arr, Acc, Arc, Acr] = twoLevelSchur(Amat, nP, nB, nC);
        
%%% Preconditioner Steps Initiated
        %%% The First stage of Preconditioning (Step 1 to 8) of Algorithm-12 (WS)
        Ds = getBDM(pn, nBG);
        Rs = getRM(GlobalBN, nBG, nB, pn);  %%% rb = getub(rb_g, m1, nB, nBG);
        rb = Rs*rb_g;                          %% Zb = rb;
        rbs = Ds*rb;
        
        RrSmat = getRrS(pn, nB);
        FrS = RrSmat*rbs;
        RcSmat = getRcS(pn, nB);
        FcS = RcSmat*rbs;
        
        out1 = (Aii)\Air;                      %% Sab = Aab - (Aai *inv(Aii)*Aib)
        Srr = Arr - (Ari*out1);                %% Srr = Arr - (Ari *inv(Aii)*Air)
        v1 = (Srr)\FrS;
        
        Scr = Acr - (Aci*out1);                %% Scr = Acr - (Aci *inv(Aii)*Air)
        dcS = FcS - (Scr*v1);
        BcSmat = getBcS(pn, nBGc);
        DcBc = (BcSmat')*dcS;
        dc = dc + DcBc;
        
        outC = (Aii)\Aic;
        Scc  = Acc - (Aci*outC);               %% Scc = Acc - (Aci *inv(Aii)*Aic)
        %%% outD = (Aii)\Aic;
        Src   = Arc - (Ari*outC);              %% Src = Arc - (Ari *inv(Aii)*Aic)
        
        temp1 = Srr\Src;
        temp2 = Scc - (Scr*temp1);
        temp3 = temp2*BcSmat;
        BcSc = (BcSmat'*temp3);
        Fcc = Fcc + BcSc;
        
%%Preconditoner for PCGM for next step: Corner-Nodes-Schur-Complement
        %%Mpre = BcSmat'*Acc*BcSmat: Reffer Waad's Thesis(Pg:147)
        Mpre1 = Acc*(BcSmat);
        Mpre2 = BcSmat'*(Mpre1);
        Mpre  = Mpre + Mpre2;
    end
    
%%%----------------------------------------------------------------------------------------
%%% PCGM with Lumped solver is called here to solve Fcc*Zc = dc
%%% The Second stage of Preconditioning (Step 9) of Algorithm-12 (WS)
    
%%% Method #1: directly solving for Fcc*Zc = dc, i.e. Zc = inv(Fcc)*dc
    % disp('Solving Corner-Nodes-Schur-Complement directly...')
    % Zc = Fcc\dc;
    
%%% Method #2 Using PCGM with Lumped Preconditioner
    disp('Solving Corner-Nodes-Schur-Complement using PCGM...')
    tol2 = 1e-4;
    maxit = 100;
    PreCon = Mpre;
    Zc = pcg(Fcc,dc,tol2,maxit,PreCon);
    
%%%----------------------------------------------------------------------------------------
%%% The Third stage of Preconditioning (Step 10 to 16) of Algorithm-12 (WS)
%%% clear rb rbs Rs Ds FrS RrSmat RcSmat Amat bvec Fcc dc
    Zb_g = zeros(nBG,1);
    for i = 1:ndom
        pn = i;                                      %% partition number
        
        [p, e, t] = getPET(pn);                      %% Mesh data: P E T
        [nP, nE, nT, nB] = getMeshDim(pn);           %% Mesh dimensions
        
        [nC, nR] = getCorRemDim(pn); 
        
        %%% Assemble Matrix & Vector for each sub-domain
        Amat = AssembleMatrix(p, e, t, 'PoissonModi', [], 0);
        bvec = AssembleVector(p, e, t, 'PoissonModi', [], 0);
        
%%% Local stiffness matrices: oneLevel: notations are similar to Waad's Thesis:
        [Aii] = oneLevelSchur(Amat, bvec, nP, nB);
        [Air, Aic, Ari, Aci, Arr, Acc, Arc, Acr] = twoLevelSchur(Amat, nP, nB, nC);
        
%%% The Third stage of Preconditioning (Step 10 to 16) of Algorithm-12 (WS)
        BcSmat = getBcS(pn, nBGc);
        ZcS = (BcSmat)*Zc;
        
        Rs = getRM(GlobalBN, nBG, nB, pn);                   %% rb = getub(rb_g, m1, nB, nBG);
        rb = Rs*rb_g;
        Ds = getBDM(pn, nBG);                  %% Zb = rb;
        rbs = Ds*rb;
        
        RrSmat = getRrS(pn, nB);
        FrS = RrSmat*rbs;
        RcSmat = getRcS(pn, nB);
        %%% FcS = RcSmat*rbs;
        
        out12 = Aii\Aic;
        Src = Arc - (Ari*out12);               %% Src = Arc - (Ari *inv(Aii)*Aic)
        v2 = FrS - (Src*ZcS);
        
        out13 = Aii\Air;
        Srr = Arr - (Ari*out13);               %% Srr = Arr - (Ari *inv(Aii)*Air)
        ZrS = Srr\v2;
        
        Zs = (RrSmat'*ZrS) + (RcSmat'*ZcS);
        ZsD = Ds*Zs;
        RZb = (Rs'*ZsD);                       %% RZb = getubg(ZsD, m1, nB, nBG);
        Zb_g = Zb_g + RZb;
    end
%%% Preconditioning of Extended Schur Complement System Ends Here.
%%%--------------------------------------------------------------------
    rho_next = dot(rb_g,Zb_g);
    if j == 1
        Pb_g = Zb_g;
    else
        beta = rho_next/rho_curr;
        Pb_g = Zb_g + (beta*Pb_g);
    end
    
%%%--------------------------------------------------------------------
    disp('Matrix-Vector Product for the Ext. Schur Complement System')
    Qb_g = zeros(nBG,1);
    for i = 1:ndom   % k --> np
        pn = i;                                      %% partition number
        
        [p, e, t] = getPET(pn);                      %% Mesh data: P E T
        [nP, nE, nT, nB] = getMeshDim(pn);           %% Mesh dimensions
        
        %%% Assemble Materix & Vector for each subdomain
        Amat = AssembleMatrix(p, e, t, 'PoissonModi', [], 0);
        bvec = AssembleVector(p, e, t, 'PoissonModi', [], 0);
        
%%% Local stiffness matrices: oneLevel: notations are similar to Waad's Thesis:        
        [Aii, Aig, Agi, Agg] = oneLevelSchur(Amat, bvec, nP, nB);
        
        Pb   = getub(Pb_g, pn, nB, nBG);
        out22 = (Aig*Pb);
        Ui = out22;
        Ui = (Aii)\Ui;
        
        out22 = (Agi*Ui);
        Qb = (Agg*Pb);
        Qb = Qb - out22;
        
        RQb = getubg(Qb, pn, nB, nBG);
        Qb_g = Qb_g + RQb;
    end
%%% Matrix-Vector Product for the Ext. Schur Complement System Ends Here.
%%%------------------------------------------------------------------------
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

dlmwrite('globalInterface_sol.dat', Ub_g, '\t');
dlmwrite('num_itterations.dat', j, '\t');
%%% Until here, the script solves for the global interface solution.

%*********************************************************************************
%%%%%% POST PROCESSING: RE-ASSEMBLING THE SOLUTION VECTOR: %%%%%%%
%clear all;
%PostAssembling

%%--------------------------------------------------------------------------------------------------------
%%************* END *****************
%%--------------------------------------------------------------------------------------------------------


%%% Useful Extra Stuff:
% Sab = Aab - (Aai *inv(Aii)*Aib)
% outA = out1; %%(Aii)\Air;
% Srr   = Arr - (Ari*outA); %% Srr = Arr - (Ari *inv(Aii)*Air)
% outB = (Aii)\Air;
% Scr   = Acr - (Aci*outA); %% Scr = Acr - (Aci *inv(Aii)*Air)


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
%      m1 = i;
%  Ub = getub(Ub_g, m1, nB, nBG);
%  out2 = (Aig*Ub);
%
%  Ui = Fi-out2;
%  Ui = inv(Aii)*Ui;
%
%  tempReadUi  = ['Ui00',num2str(m1)];
%  dlmwrite([tempReadUi '.dat'], Ui, '\t');
%
%  tempReadUb_g  = ['Ub_g00',num2str(m1)];
%  dlmwrite([tempReadUb_g '.dat'], Ub_g, '\t');
%
%  Ub_gFinal = Ub_gFinal + Ub_g;
% % end
%%--------------------------------------------------------------------------------------------------------
%%************* END *****************
%%--------------------------------------------------------------------------------------------------------
