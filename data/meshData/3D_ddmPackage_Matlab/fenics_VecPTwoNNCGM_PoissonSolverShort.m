%% NNC/BDDC PCGM-DDM Two-Level NNC Preconditined CGM for       
%%- 3D Linear Elasticity Equation (Vector Valued PDE).  
%
% Copyright (C) 2017 Ajit Desai, Ph.D. Candidate ajit.ndesai@gmail.com
%
%%% Two-Level Neumann-Neumann Preconditioning of Schur Complement System
%%% For detail refer Algorithm-12 of WS_Thesis (Pg-146)
%%% For details for the PCGM algorithm refer to Appendix-D or the MK's Thesis

%%--------------------------------------------------------------------------------
%*********************************************************************************
%%%%%%% MESH DECOMPOSITION: PREPROCESSING: %%%%%%%
clearvars;
% twoLevelPreprocessor

%*********************************************************************************
%%%%%%% MAIN SOLVER: PROCESSING: %%%%%%%
disp('Calculating The First Iteration Value (r_0)...')

%%%%%%% Load the required global data %%%%%%%
GlobalBN = dlmread('../boundary_nodes.dat', '');        %% Global boundary nodes (Global interface)
GlobalCN = dlmread('../corner_nodes.dat', '');          %% Global corner nodes
ndom = dlmread('../num_partition.dat', '');             %% number of partitions
nBG  = length(GlobalBN);
nBGc = length(GlobalCN);

nVec = 3;                                                %% number of vector components  

%%%%%%% r_0 = G_tau - S*X_tau0 %%%%%%%
%%% Since X_tau0 = 0, Initial Guess, r_0 = G_tau
rb_g = zeros(nVec*nBG,1);   

for i = 1:ndom
    pn = i;                                      %% partition number
    
    [nP, nE, nT, nH, nB] = getMeshDim(pn);       %% Mesh dimensions
    
    %%%%%%% Restriction matrix (Scatter & gather operator)
    Rs = zeros(nVec*nB,nVec*nBG);
    for j = 1:nVec     %% j=1:3 for 3D and j=1:2 for 2D
        Rmat = getRM(GlobalBN, nBG, nB, pn);  
        Rs( (j-1)*nB+1:j*nB, (j-1)*nBG+1:j*nBG ) = Rmat(1:nB, 1:nBG);
    end       
   
%     Dmat = getBDM(pn, nBG);
    %%%%%%% Assemble matrix FEniCS/dolfin-Assembled Mat-Vecs %%%%%%%%
% % %     tempAb  = ['../../outputs/Ab',num2str(pn)];
% % %     load(tempAb)
% % %     Amat = An;
% % %     bvec = bn';
% % %     
% % %     %%% Local stiffness matrices: oneLevel: notations are similar to Waad's Thesis:
% % %     [Aii, Aig, Agi, Agg, Fi, Fg] = oneLevelVecSchur(Amat, bvec, nP, nB, nVec);
    
    
    % Code by sudhi -Loading all the deterministic 3D decomposed matrices
    subpath = ['../../../external/dolfin/data/Amats/subdom',num2str(pn)];
    filepath = strcat(subpath,'/ADgg1.dat');
    Agg = load(filepath);
    filepath = strcat(subpath,'/ADii1.dat');
    Aii = load(filepath);
    filepath = strcat(subpath,'/ADig1.dat');
    Aig = load(filepath);
    filepath = strcat(subpath,'/bg1.dat');
    Fg = load(filepath);
    filepath = strcat(subpath,'/bi1.dat');
    Fi = load(filepath);
    
    
    Agi = Aig';
    
    
    
    out2 = (Aii)\Fi;
    Gi = Agi*out2;
    Gi = Fg - Gi;
    RGi = Rs'*(Gi);
    rb_g = rb_g + RGi;
end

%*********************************************************************************
%%%%%%% PCGM Solver %%%%%%%
disp('TwoLevel-NNC-PCGM Solution Initiated...')
Ub_g = zeros(nVec*nBG,1);      

tol = 10^(-8);                             %% Tolerance limit
maxiter = 100;                             %% Max Iteration
for j = 1:maxiter
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if j==0
       disp('Two-Level NNC/BDDC Preconditioning of Ext. Schur Complement System')
    end
    dc  = zeros(nVec*nBGc,1);
    Fcc = zeros(nVec*nBGc,nVec*nBGc);
    Mpre = zeros(nVec*nBGc,nVec*nBGc);
    %%% Zb_g = zeros(nBG,1); %%% RZb  = zeros(nBG,1);
    for i = 1:ndom
        pn = i;                                      %% partition number

        [nP, nE, nT, nH, nB] = getMeshDim(pn);       %% Mesh dimensions
        [nC, nR] = getCorRemDim(pn);                 %% Mesh dim: corner & remaining
        
        
        
%%% Assemble Matrix & Vector for each sub-domain
        %%: FEniCS/dolfin-Assembled Mat-Vecs
% % %         tempAb  = ['../outputs/Ab',num2str(pn)];
% % %         load(tempAb)
% % %         Amat = An;
% % %         bvec = bn';
% % %         
% % % %%% Local stiffness matrices: oneLevel: notations are similar to Waad's Thesis:
% % %         [Aii] = oneLevelVecSchur(Amat, bvec, nP, nB, nVec);
% % %         [Air, Aic, Ari, Aci, Arr, Acc, Arc, Acr] = twoLevelVecSchur(Amat, nP, nB, nC, nVec);
        
        % Code by sudhi -Loading all the deterministic 3D decomposed matrices
        subpath = ['../../../external/dolfin/data/Amats/subdom',num2str(pn)];
        filepath = strcat(subpath,'/ADii1.dat');
        Aii = load(filepath);
        filepath = strcat(subpath,'/ADir1.dat');
        Air = load(filepath);
        filepath = strcat(subpath,'/ADic1.dat');
        Aic = load(filepath);
        filepath = strcat(subpath,'/ADcc1.dat');
        Acc = load(filepath);
        filepath = strcat(subpath,'/ADrr1.dat');
        Arr = load(filepath);
        filepath = strcat(subpath,'/ADrc1.dat');
        Arc = load(filepath);
        
        Ari = Air';
        Aci = Aic';
        Acr = Arc';
        
        
        
        
%%% Preconditioner Steps Initiated
        %%% The First stage of Preconditioning (Step 1 to 8) of Algorithm-12 (WS)
        
        %%: Get all restriction operators
        %Ds = getBDM(pn, nBG);
        Ds = zeros(nVec*nB, nVec*nB);
        for k = 1:nVec     %% k=1:3 for 3D and k=1:2 for 2D
            Dmat = getBDM(pn, nBG);  
            Ds( (k-1)*nB+1:k*nB, (k-1)*nB+1:k*nB ) = Dmat(1:nB, 1:nB);
        end
        
        %Rs = getRM(GlobalBN, nBG, nB, pn);     %% rb = getub(rb_g, m1, nB, nBG);
        Rs = zeros(nVec*nB,nVec*nBG);
        for k = 1:nVec     %% k=1:3 for 3D and k=1:2 for 2D
            Rmat = getRM(GlobalBN, nBG, nB, pn);  
            Rs( (k-1)*nB+1:k*nB, (k-1)*nBG+1:k*nBG ) = Rmat(1:nB, 1:nBG);
        end
        
        %RrSmat = getRrS(pn, nB);
        RrSmat = zeros(nVec*nR, nVec*nB);
        for k = 1:nVec     %% k=1:3 for 3D and k=1:2 for 2D
            RrS = getRrS(pn, nB);
            RrSmat( (k-1)*nR+1:k*nR, (k-1)*nB+1:k*nB ) = RrS(1:nR, 1:nB);
        end
        
        %RcSmat = getRcS(pn, nB);
        RcSmat = zeros(nVec*nC, nVec*nB);
        for k = 1:nVec     %% k=1:3 for 3D and k=1:2 for 2D
            RcS = getRcS(pn, nB);  
            RcSmat( (k-1)*nC+1:k*nC, (k-1)*nB+1:k*nB ) = RcS(1:nC, 1:nB);
        end
        
        %BcSmat = getBcS(pn, nBGc);
        BcSmat = zeros(nVec*nC, nVec*nBGc);
        for k = 1:nVec     %% k=1:3 for 3D and k=1:2 for 2D
            BcS = getBcS(pn, nBGc);
            BcSmat( (k-1)*nC+1:k*nC, (k-1)*nBGc+1:k*nBGc ) = BcS(1:nC, 1:nBGc);
        end
        
        rb = Rs*rb_g;                          %% Zb = rb;
        rbs = Ds*rb;
        
        FrS = RrSmat*rbs;
        FcS = RcSmat*rbs;
        
        out1 = (Aii)\Air;                      %% Sab = Aab - (Aai *inv(Aii)*Aib)
        Srr = Arr - (Ari*out1);                %% Srr = Arr - (Ari *inv(Aii)*Air)
        v1 = (Srr)\FrS;
        
        Scr = Acr - (Aci*out1);                %% Scr = Acr - (Aci *inv(Aii)*Air)
        dcS = FcS - (Scr*v1);
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
    
%%%-Method #1: directly solving for Fcc*Zc = dc, i.e. Zc = inv(Fcc)*dc
    disp('Solving Corner-Nodes-Schur-Complement directly...')
    Zc = Fcc\dc;
    
% %%%-Method #2 Using PCGM with Lumped Preconditioner
%     if j==0
%         disp('Solving Corner-Nodes-Schur-Complement using Matlab-PCGM...')
%     end
%     tol2 = 1e-4;
%     maxit = 100;
%     PreCon = Mpre;
%     Zc = pcg(Fcc,dc,tol2,maxit,PreCon);
    
%%%----------------------------------------------------------------------------------------
%%% The Third stage of Preconditioning (Step 10 to 16) of Algorithm-12 (WS)
%%% clear rb rbs Rs Ds FrS RrSmat RcSmat Amat bvec Fcc dc
    Zb_g = zeros(nVec*nBG,1);
    for i = 1:ndom
        pn = i;                                      %% partition number
        
        %[p, e, t, h] = getPETH(pn);                 %% Mesh data: P E T H
        [nP, nE, nT, nH, nB] = getMeshDim(pn);       %% Mesh dimensions
        [nC, nR] = getCorRemDim(pn); 
        
        %%% FEniCS/dolfin-Assembled Mat-Vecs for each sub-domain %%%%
% % % %         tempAb  = ['../outputs/Ab',num2str(pn)];
% % % %         load(tempAb)
% % % %         Amat = An;
% % % %         bvec = bn';
% % % %         
% % % % %%% Local stiffness matrices: oneLevel: notations are similar to Waad's Thesis:
% % % %         [Aii] = oneLevelVecSchur(Amat, bvec, nP, nB, nVec);
% % % %         [Air, Aic, Ari, Aci, Arr, Acc, Arc, Acr] = twoLevelVecSchur(Amat, nP, nB, nC, nVec);
        
        
        
        
        % Code by sudhi -Loading all the deterministic 3D decomposed matrices
        subpath = ['../../../external/dolfin/data/Amats/subdom',num2str(pn)];
        filepath = strcat(subpath,'/ADii1.dat');
        Aii = load(filepath);
        filepath = strcat(subpath,'/ADir1.dat');
        Air = load(filepath);
        filepath = strcat(subpath,'/ADic1.dat');
        Aic = load(filepath);
        filepath = strcat(subpath,'/ADcc1.dat');
        Acc = load(filepath);
        filepath = strcat(subpath,'/ADrr1.dat');
        Arr = load(filepath);
        filepath = strcat(subpath,'/ADrc1.dat');
        Arc = load(filepath);
        
        Ari = Air';
        Aci = Aic';
        Acr = Arc';
        
        
%%% The Third stage of Preconditioning (Step 10 to 16) of Algorithm-12 (WS)
%%: Get all restriction operators
        %Ds = getBDM(pn, nBG);
        Ds = zeros(nVec*nB, nVec*nB);
        for k = 1:nVec     %% k=1:3 for 3D and k=1:2 for 2D
            Dmat = getBDM(pn, nBG);  
            Ds( (k-1)*nB+1:k*nB, (k-1)*nB+1:k*nB ) = Dmat(1:nB, 1:nB);
        end
        
        %Rs = getRM(GlobalBN, nBG, nB, pn);     %% rb = getub(rb_g, m1, nB, nBG);
        Rs = zeros(nVec*nB,nVec*nBG);
        for k = 1:nVec     %% k=1:3 for 3D and k=1:2 for 2D
            Rmat = getRM(GlobalBN, nBG, nB, pn);  
            Rs( (k-1)*nB+1:k*nB, (k-1)*nBG+1:k*nBG ) = Rmat(1:nB, 1:nBG);
        end
        
        %RrSmat = getRrS(pn, nB);
        RrSmat = zeros(nVec*nR, nVec*nB);
        for k = 1:nVec     %% k=1:3 for 3D and k=1:2 for 2D
            RrS = getRrS(pn, nB);
            RrSmat( (k-1)*nR+1:k*nR, (k-1)*nB+1:k*nB ) = RrS(1:nR, 1:nB);
        end
        
        %RcSmat = getRcS(pn, nB);
        RcSmat = zeros(nVec*nC, nVec*nB);
        for k = 1:nVec     %% k=1:3 for 3D and k=1:2 for 2D
            RcS = getRcS(pn, nB);  
            RcSmat( (k-1)*nC+1:k*nC, (k-1)*nB+1:k*nB ) = RcS(1:nC, 1:nB);
        end
        
        %BcSmat = getBcS(pn, nBGc);
        BcSmat = zeros(nVec*nC, nVec*nBGc);
        for k = 1:nVec     %% k=1:3 for 3D and k=1:2 for 2D
            BcS = getBcS(pn, nBGc);
            BcSmat( (k-1)*nC+1:k*nC, (k-1)*nBGc+1:k*nBGc ) = BcS(1:nC, 1:nBGc);
        end

        ZcS = (BcSmat)*Zc;
        rb = Rs*rb_g;
        rbs = Ds*rb;
        
        FrS = RrSmat*rbs;
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
    if j==0
    disp('Matrix-Vector Product for the Ext. Schur Complement System')
    end
    Qb_g = zeros(nVec*nBG,1);
    for i = 1:ndom   % k --> np
        pn = i;                                      %% partition number
        
        %[p, e, t] = getPET(pn);                      %% Mesh data: P E T H
        [nP, nE, nT, nH, nB] = getMeshDim(pn);        %% Mesh dimensions
        
        %%% FEniCS/dolfin Assemble Materix & Vector for each subdomain
% % % % % %         tempAb  = ['../outputs/Ab',num2str(pn)];
% % % % % %         load(tempAb)
% % % % % %         Amat = An;
% % % % % %         bvec = bn';
% % % % % %         
% % % % % % %%% Local stiffness matrices: oneLevel: notations are similar to Waad's Thesis:        
% % % % % %         [Aii, Aig, Agi, Agg] = oneLevelVecSchur(Amat, bvec, nP, nB, nVec);
% % % % % %         

% Code by sudhi -Loading all the deterministic 3D decomposed matrices
    subpath = ['../../../external/dolfin/data/Amats/subdom',num2str(pn)];
    filepath = strcat(subpath,'/ADgg1.dat');
    Agg = load(filepath);
    filepath = strcat(subpath,'/ADii1.dat');
    Aii = load(filepath);
    filepath = strcat(subpath,'/ADig1.dat');
    Aig = load(filepath);
    filepath = strcat(subpath,'/bg1.dat');
    Fg = load(filepath);
    filepath = strcat(subpath,'/bi1.dat');
    Fi = load(filepath);
    
    
    Agi = Aig';


        %Rs = getRM(GlobalBN, nBG, nB, pn);     %% rb = getub(rb_g, m1, nB, nBG);
        Rs = zeros(nVec*nB,nVec*nBG);
        for k = 1:nVec     %% k=1:3 for 3D and k=1:2 for 2D
            Rmat = getRM(GlobalBN, nBG, nB, pn);  
            Rs( (k-1)*nB+1:k*nB, (k-1)*nBG+1:k*nBG ) = Rmat(1:nB, 1:nBG);
        end
        
        Pb = Rs*Pb_g;       %Pb   = getub(Pb_g, pn, nB, nBG);
        out22 = (Aig*Pb);
        Ui = out22;
        Ui = (Aii)\Ui;
        
        out22 = (Agi*Ui);
        Qb = (Agg*Pb);
        Qb = Qb - out22;
        
        RQb = Rs'*Qb;        %RQb = getubg(Qb, pn, nB, nBG);
        Qb_g = Qb_g + RQb;
    end
%%% Matrix-Vector Product for the Ext. Schur Complement System Ends Here.
%%%------------------------------------------------------------------------
    rho_curr = rho_next;
    alpha = rho_curr/(dot(Qb_g,Pb_g));
    err = (alpha*alpha*(dot(Pb_g,Pb_g)))/(dot(Ub_g,Ub_g));
    Ub_g = Ub_g + (alpha*Pb_g);
    rb_g = rb_g - (alpha*Qb_g);
    
    dispout = ['Global-Interface PCGM with NNC/BDDC Preconditioner Iteration #',pad(num2str(j),2,'left','0'), ', Relative error of: ', num2str(err)];
    disp(dispout);
    if err < tol
        break
    end
end

solout = '../../../data/solution/global_interface_matlab.dat';
dlmwrite(solout, Ub_g,'\t');
dlmwrite('../../../data/solution/num_iterations.dat', j,'\t');

dispout = ['Interface Solution: ', solout];
disp(dispout)

%%% Until here, the script solves for the global interface solution.

disp('Next: run "PostAssemblingVec" ')

%*********************************************************************************
%%%%%% POST PROCESSING: RE-ASSEMBLING THE SOLUTION VECTOR: %%%%%%%
%clear all;
%PostAssembling

%%--------------------------------------------------------------------------------------------------------
%%************* END *****************
%%--------------------------------------------------------------------------------------------------------