%% Finding the local solutions Boundary and Interior solutions :  Date/Author July/2013/Ajit 
%%% Refer Waad's/MK's Thesis for the analytical formulation

%%-----------------------------------------------------------------------------------------------
clearvars;
disp('For Dolfin-Assembled Mat-Vec: UVW-Vec-Post-Assembling Initiated...')
Ub_g = dlmread('../../../data/solution/global_interface_matlab.dat', ''); %% dolfin-Global interface solution vector

%%%%%%% Loading the requied data %%%%%%%
npoints = dlmread('../points.dat','');   %dlmread('points.dat','');
GlobalBN = dlmread('../boundary_nodes.dat', '');
ndom = dlmread('../num_partition.dat', '');
nitter = dlmread('../../../data/solution/num_iterations.dat', '');
% if nitter ~= 0 
% disp('number of itterations'); disp(nitter)
% end

%%%%%%% Initiation of solution vector %%%%%%%
nBG = length(GlobalBN);                          %% Number of interface nodes
nPG = length((npoints));                         %% Number of nodes
nVec = 3;
U = zeros(nVec*nPG,1);


for i = 1:ndom   %% k --> np
    pn = i;                                      %% partition number
    
    %[p, e, t] = getPET(pn);                     %% Mesh data: P E T
    [nP, nE, nT, nH, nB] = getMeshDim(pn);       %% Mesh dimensions
    nI = nP-nB;                                  %% Interior nodes

    %%%%%%% Restriction matrix (Scatter & gather operator)
    %%%%%%% Restriction matrix (Scatter & gather operator)
    Rmat = zeros(nVec*nB,nVec*nBG);
    for j = 1:nVec     %% j=1:3 for 3D and j=1:2 for 2D
        Rmatc = getRM(GlobalBN, nBG, nB, pn);  
        Rmat( (j-1)*nB+1:j*nB, (j-1)*nBG+1:j*nBG ) = Rmatc(1:nB, 1:nBG);
    end      
    %%%%%%% Assemble matrix & Assemble vector %%%%%%%
    %%: FEniCS/dolfin-Assembled Mat-Vecs
% % % %     tempAb  = ['../outputs/Ab',num2str(pn)];
% % % %     load(tempAb)
% % % %     Amat = An;
% % % %     bvec = bn';
% % % %     
% % % %     [Aii, Aig, Agi, Agg, Fi, Fg] = oneLevelVecSchur(Amat, bvec, nP, nB, nVec);


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
    
%%%%%%% Local interface solution Vector %%%%%%%
% Thesis Form: (u_tau)^s = R_s * (u_tau) :: (u_tau)^s = Ub
    %Ub = getub(Ub_g, pn, nB, nBG);   
    Ub = (Rmat*Ub_g);                      % For Direct DDM

%%%%%%% Local Interior Solution Vector: %%%%%%%
% Thesis Form: (A_ii)^s *(u_i)^s = (f_i)^s - {(A_taui)^s * (R_s) * (u_tau)}
% (u_i)^s = UI :: Local solution for each sub-domain
% Fortran code: solve((np-nb)*npceout,1,Aii,Ui,Fi-matmul(Aig,Ub))
    Ui   = (Aii)\(Fi - (Aig*Ub));   %%%out2 = (Aig*Ub); %%Ui = Fi-out2; %%Ui = inv(Aii)*Ui;
    
%%%%%%% Re-numbering to previously original nodes & %%%%%%% 
%- Constructing final solution vector:
    npi = length(Ui);                            %% npi = number of points in each sub-domain
    %%% NOTE: Here we need to use re-arranged nodes (Nnondes0*) 
    %%% because we use nbnodes while solving
    %tempNi2  = ['../data/Nnodes00',num2str(pn)];         %%['nodes00',num2str(m1)];
    tempNi2  = ['../nodes',pad(num2str(pn),4,'left','0')];
    nodei = dlmread([tempNi2 '.dat']);
    
    %%- Scalar Assembly 
    %for j = 1:npi
    %    U(nodei(j)) = Ui(j);
    %end
    
    %%- 3D Vec Assembly 
    for k = 1:npi
        if k <= nI
            U(nodei(k)) = Ui(k);
        elseif k <= 2*nI
            U(nPG+nodei(k-nI)) = Ui(k);
        else
            U((2*nPG)+nodei(k-(2*nI))) = Ui(k);
        end
    end     
end

%%- Scalar Assembly 
%  nBG = length(Ub_g);
%  for j = 1:nBG
%    tempGB = GlobalBN(j);
%    U(tempGB) = Ub_g(j);
%  end
 
%%- 3D Vec Assembly 
 nBGc = length(Ub_g);
 for j = 1:nBGc
     if j <= nBG 
        tempGB = GlobalBN(j);
        U(tempGB) = Ub_g(j);
     elseif j <= 2*nBG
        tempGB = nPG+GlobalBN(j-nBG);
        U(tempGB) = Ub_g(j);
     else
        tempGB = (2*nPG)+GlobalBN(j-(2*nBG));
        U(tempGB) = Ub_g(j);
     end
 end
 
 strout = '../../solution/final_solutionMatlab.dat';
 dlmwrite(strout, U, '\t');
 dispout = ['Final Solution: ', strout];
 disp(dispout)
  
 disp('Next: run "myWriteVecVtk" ')

%%--------------------------------------------------------------------------------------------------------
%%-********* PLOTTING SOLUTION VECTOR *************
%%--------------------------------------------------------------------------------------------------------
% %%%%%%% DDM Solution %%%%%%%
% P = npoints';
% T = triangles';
% if dolfin == 1
%     figure(2); clf
% else
%     figure(1); clf
% end
% pdesurf(P,T,U)
% shading faceted
% title('DDM Computed solution');
% axis tight %% axis([0 2 0 1 -1.2 1.2])   %% clc;
% 
% % %%%%%%% Computing & Plotting exact solution %%%%%%%
% % uE = zeros(size(U));
% % for i = 1:size(P,2)
% %   xE = P(:,i);
% %   uE(i) = sin(pi*xE(1)) * sin(2*pi*xE(2)); 
% %   %uE(i) = 1 + xE(1)^2 + (2*xE(2)^2);
% % end
% % 
% % figure(3); clf
% % pdesurf(P,T,uE);
% % shading faceted
% % title('Exact solution')
% % axis tight %% axis([0 2 0 1 -1.2 1.2])
% 
% %%%%%%% Compute the error %%%%%%%
% %error = U - u;
% %enorm = max(abs(error));
% %disp(['Maximum norm error: ' num2str(enorm)]);

%%--------------------------------------------------------------------------------------------------------
%%-********* END *************
%%--------------------------------------------------------------------------------------------------------

