%% Intrusive PCE Code to solve Ax = b
%%: Inputs are charectorized by Log-Normal RVs
%%: PCE-LNRV :: K = K_i Psi_i  :: Kml*[K0*1 + K1*x + K2*(x^2-1)]
%%: Kmg, Ksg :: Mean and SD of Underlaying Gaussian RVs
%%: A = deterministic finite element matrix : Puffin Assembly
%%: b = deterministic source vector         : Puffin Assembly
%%: As = Stochastic finite element matrix
%%: Fs = Stochastic source vector
%%: u = Sum(ui*Pis_i)
%%: Cijk = <Psi_i Psi_j Psi_k> :: Calculated using moment_ND.m
%
%   Ver0.0 : January 05, 2017 : Ajit Desai : Slow Speed (Full Cijk)
%                                          : Renamed to Kuf*_slow.m
%   Ver0.1 : January 10, 2017 : Ajit Desai : Medium Speed (NonZero Cijk)
%                                          : Renamed to Kuf*_medium.m
%   Ver1.0 : January 17, 2017 : Ajit Desai : Fast (nonZero Ckij)
%
%   Ver1.1 : May 02, 2017     : Ajit Desai : Modified to include Stochastic
%                                        Process instead of Random Variable
%
%%: Steps to follow
% 1: open 'mySquare.geo' in editor: adjust 'lc' for mesh density(0.122=133)
% 2: open 'mySquare.geo' in GMSH: generate mesh and save as 'mySquare.msh'
% 3: From MATLAB run 'Preprocessor.m': to generate necessory input files
% 4: Run this code: outputs will be save in 'Dim*_p*.dat' file
%%:------------------------------------------------------------------------
clearvars
tic

disp('Cleaning previous data...')
!rm -rf ../data/Amats/subdom*
!rm -rf ../data/Bvecs/*.dat

%%: Read mesh data :
ndom = dlmread('../data/num_partition.dat', '');

for pn = 1:ndom     %% pn: partition number
    %%: Points
    tempPt2  = ['../data/rpoints00',num2str(pn)];
    tempRead1 = dlmread([tempPt2 '.dat']);
    p = tempRead1';
    %%: Elements
    tempEg2 = ['../data/redges00',num2str(pn)];
    test2 = load([tempEg2 '.dat']);
    test3 = isempty(test2);
    if test3 == 1
        e = [];
    else
        tempRead2 = dlmread([tempEg2 '.dat']);
        e = tempRead2';
    end
    %%: Triangls
    tempT2  = ['../data/rtriangles00',num2str(pn)];
    tempRead3 = dlmread([tempT2 '.dat']);
    t = tempRead3';
    
    %%: Mesh-dimensions
    tempMd2  = ['../data/meshdim00',num2str(pn)]; % meshdim = [nPnt nEdg nTra nPar]
    meshdim12 = dlmread([tempMd2 '.dat']);
    
    nP = meshdim12(1);                    % number of points
    nE = meshdim12(2);                    % number of edges
    nT = meshdim12(3);                    % number of triangles
    nB = meshdim12(4);                    % number of global boundary nodes
    
    
    %%% Assemble Matrix & Vector for each sub-domain
    Amat = AssembleMatrix(p, e, t, 'PoissonModi', [], 0);
    bvec = AssembleVector(p, e, t, 'PoissonModi', [], 0);
    
    %%% Local stiffness matrices: oneLevel: notations are similar to Waad's Thesis:
    [Aii, Aig, Agi, Agg, Fi, Fg] = oneLevelSchur(Amat, bvec, nP, nB);
    
    
    if (pn < 10)
        createD = horzcat('../data/Amats/subdom000',num2str(pn));
    else
        createD = horzcat('../data/Amats/subdom00',num2str(pn));
    end
        mkdir (createD)
        APii  = [createD,'/Aii'];
        dlmwrite([APii '.dat'], Aii, '\t');
        APgi  = [createD,'/Agi'];
        dlmwrite([APgi '.dat'], Agi, '\t');
        APgg  = [createD,'/Agg'];
        dlmwrite([APgg '.dat'], Agg, '\t');
       
end

%%: Deterministic Vector Assebly (f=1)
% b = AssembleVector(p, e, t, 'PoissonModi', [], 0, 0);

%Bvecs  = ['../data/Bvecs/Bvecs00',num2str(pn)];
%dlmwrite([Bvecs '.dat'], b, '\t');

%%: Solving Intrusive System of Equations
%u_intru = As\fs;


toc

%%-------------------------------------------------------------------------



%% Intrusive PCE Code to solve Ku = f : Gaussian RV
% %%: GRV :: K = K0 PS0 + K1 PS1 :: K0*1 + K1*x = Kmn + Ksg*x
% %%: Kmg, Ksg :: Mean and SD of Underaying Gaussian RV
% %%: f = 1 :: Constant
% %%: u = SUM[ui*PSi] :: PSi's are Standard Hermite Polynomial
% %%: Cijk = <PSi PSj PSk> :: Calulated using trippleProd.m
% %%: nord - Maximum order of expansion to represent solution (u)
% %%-------------------------------------------------------------------------
%
% clear all
% nord = 8;
% nkle = 2;
% npce = nord+1;
% u = zeros(npce,nord);
%
%     str = strcat('cijk_L2_8');
%     load(str)
%
% for p = 2:nord
%     nP = p+1;
%     f = zeros(nP,1);
%     KL = zeros(nP,1);
%
%     f(1) = 1;
%     KL(1) = 1;
%     KL(2) = 0.2;
%
%     for i=1:nP
%         for j=1:nP
%             aa = 0;
%             for k=1:nkle
%                 Kijk = cijk(i,j,k)*KL(k);
%                 aa = aa+Kijk;
%
%             end
%             KK(i,j) = aa;
%
%         end
%     end
%
%     u_intru = inv(KK)*f;
%     u(1:nP,p) = u_intru;
% end
% u = u';
% u = u(2:end,:)
%
%%-------------------------------------------------------------------------

