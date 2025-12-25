%% PCGM-DDM (Preconditioned Conjugate Gradient Method) for 
%%- 3D Linear Elasticity Equation (Vector Valued PDE).  
%
% Copyright (C) 2017 Ajit Desai, Ph.D. Candidate ajit.ndesai@gmail.com
%
%%% One Level Preconditioning of Schur Complement System
%%% PCGM algorithm refer to Appendix-D.1 (page 174) of the MK's Thesis

%%--------------------------------------------------------------------------------
%*********************************************************************************
clearvars; clc;
disp('One-Level PCGM-DDM for 3D-Elasticity Initiated...')
tic

%*********************************************************************************
%%%%%%% MAIN SOLVER: PROCESSING: %%%%%%%
disp('Calculating The First Iteration Value (r_0)...')

%%%%%%% Load the required global data %%%%%%%
GlobalBN = dlmread('../boundary_nodes.dat', '');    %% Global boundary nodes (Global interface)
ndom = dlmread('../num_partition.dat', '');         %% number of partitions
nBG  = length(GlobalBN);                                 %% number of global boundary nodes

nVec = 3;                                                %% number of vector components  

%%%%%%% r_0 = G_tau - S*X_tau0: First itteration value %%%%%%%
%%% Since X_tau0 = 0, Initial Guess, r_0 = G_tau
rb_g = zeros(nVec*nBG,1);     %% RGi  = zeros(nBG,1);
for i = 1:ndom
    pn = i;                                      %% partition number
    
    [nP, nE, nT, nH, nB] = getMeshDim(pn);       %% Mesh dimensions
        
    %%%%%%% Restriction matrix (Scatter & gather operator)
    Rs = zeros(nVec*nB,nVec*nBG);
    for j = 1:nVec     %% j=1:3 for 3D and j=1:2 for 2D
        Rmat = getRM(GlobalBN, nBG, nB, pn);  
        Rs( (j-1)*nB+1:j*nB, (j-1)*nBG+1:j*nBG ) = Rmat(1:nB, 1:nBG);
    end       
    
    %%: FEniCS/dolfin-Assembled Mat-Vecs
    tempAb  = ['../outputs/Ab',num2str(pn)];
    load(tempAb)
    Amat = An;
    bvec = bn';
    
    %%: Local stiffness matrices: oneLevel: notations are similar to Waad's Thesis:
    [Aii, Aig, Agi, Agg, Fi, Fg] = oneLevelVecSchur(Amat, bvec, nP, nB, nVec);
    
    out2 = (Aii)\Fi;
    Gi = Agi*out2;
    Gi = Fg - Gi;
    RGi = Rs'*(Gi);
    rb_g = rb_g + RGi;
end

%%-------------------------------------------------------------------------
%%%%%%% PCGM Solver %%%%%%%
disp('PCGM Solution Initiated...')
Ub_g = zeros(nVec*nBG,1);       %% Zb_g = zeros(nBG,1);  % Qb_g = zeros(nBG,1);
tol = 10^(-9);             %%% Tolerance limit
maxiter = 200;             %%% Max Iterations

for j = 1:maxiter
    if j==0
        disp('Preconditioning of Ext. Schur Complement System')
    end
    Zb_g = zeros(nVec*nBG,1);
    for i = 1:ndom
        pn = i;                                      %% partition number
        
        [nP, nE, nT, nH, nB] = getMeshDim(pn);       %% Mesh dimensions

        %%: FEniCS/dolfin-Assembled Mat-Vecs
        tempAb  = ['../outputs/Ab',num2str(pn)];
        load(tempAb)
        Amat = An;
        bvec = bn';
        
        [Aii, Aig, Agi, Agg, Fi, Fg] = oneLevelVecSchur(Amat, bvec, nP, nB, nVec);
        
        Rs = zeros(nVec*nB,nVec*nBG);
        for k = 1:nVec     %% k=1:3 for 3D and k=1:2 for 2D
            Rmat = getRM(GlobalBN, nBG, nB, pn);  
            Rs( (k-1)*nB+1:k*nB, (k-1)*nBG+1:k*nBG ) = Rmat(1:nB, 1:nBG);
        end     
        
        rb = Rs*rb_g;      %% rb = getub(rb_g, pn, nB, nBG);
        
        Zb = (Agg)\rb;
        RZb = Rs'*Zb;      %% RZb = getubg(Zb, pn, nB, nBG);

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
    if j==0
    disp('Matrix-Vector Product for the Ext. Schur Complement System')
    end
    Qb_g = zeros(nVec*nBG,1);
    for i = 1:ndom   % k --> np
        pn = i;                                      %% partition number
        
        [nP, nE, nT, nH, nB] = getMeshDim(pn);       %% Mesh dimensions

        %%: FEniCS/dolfin-Assembled Mat-Vecs
        tempAb  = ['../outputs/Ab',num2str(pn)];
        load(tempAb)
        Amat = An;
        bvec = bn';
        
        [Aii, Aig, Agi, Agg, Fi, Fg] = oneLevelVecSchur(Amat, bvec, nP, nB, nVec);
        
        Rs = zeros(nVec*nB,nVec*nBG);
        for k = 1:nVec     %% k=1:3 for 3D and k=1:2 for 2D
            Rmat = getRM(GlobalBN, nBG, nB, pn);  
            Rs( (k-1)*nB+1:k*nB, (k-1)*nBG+1:k*nBG ) = Rmat(1:nB, 1:nBG);
        end    
        
        Pb = Rs*Pb_g;       
        out2 = (Aig*Pb);
        Ui = (Aii)\out2;
        
        out2 = (Agi*Ui);
        Qb = (Agg*Pb);
        Qb = Qb - out2;
        
        RQb = Rs'*Qb;       
        Qb_g = Qb_g + RQb;
    end
%%% Matrix-Vector Product for the Ext. Schur Complement System Ends Here.
%%-------------------------------------------------------------------------
    
    rho_curr = rho_next;
    alpha = rho_curr/(dot(Qb_g,Pb_g));
    err = (alpha*alpha*(dot(Pb_g,Pb_g)))/(dot(Ub_g,Ub_g));
    Ub_g = Ub_g + (alpha*Pb_g);
    rb_g = rb_g - (alpha*Qb_g);
    
    dispout = ['PCGM iteration: ',pad(num2str(j),2,'left','0'), ', Relative error of: ', num2str(err)];
    disp(dispout);
    if err < tol
        break
    end
end
toc

solout = '../data/dolfin_globalInterface_sol.dat';
dlmwrite(solout, Ub_g, '\t');
dlmwrite('../data/num_itterations.dat', j, '\t');

dispout = ['Interface Solution: ', solout];
disp(dispout)

disp('Next: run "PostAssemblingVec" ')

%t2 = toc;
%fprintf('Elapsed time is %0.2f seconds.\n', t2);

%%% Until here, the script solves for the global interface solution.
%*********************************************************************************

%%%%%% POST PROCESSING: RE-ASSEMBLING THE SOLUTION VECTOR: %%%%%%%
% clearvars;
% PostAssembling
%%--------------------------------------------------------------------------------------------------------
%%************* END *****************
%%--------------------------------------------------------------------------------------------------------