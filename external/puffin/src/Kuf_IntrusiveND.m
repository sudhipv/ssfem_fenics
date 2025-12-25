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

for pn = 1:4     %% pn: partition number
    %%: Points
    tempPt2  = ['../../../data/meshData/points000',num2str(pn)]; %%FortraData
    %tempPt2  = ['../data/rpoints00',num2str(pn)]    %% MatlabData
    tempRead1 = dlmread([tempPt2 '.dat']);
    p = tempRead1';
    %%: Elements 
    tempEg2 = ['../../../data/meshData/edges000',num2str(pn)];
    %tempEg2 = ['../data/redges00',num2str(pn)];
    test2 = load([tempEg2 '.dat']);
    test3 = isempty(test2);
    if test3 == 1
        e = [];
    else
        tempRead2 = dlmread([tempEg2 '.dat']);
        e = tempRead2';
    end
    %%: Triangls 
    tempT2  = ['../../../data/meshData/triangles000',num2str(pn)];
    %tempT2  = ['../data/rtriangles00',num2str(pn)];
    tempRead3 = dlmread([tempT2 '.dat']);
    t = tempRead3';
    
    %%: Mesh-dimensions
    tempMd2  = ['../../../data/meshData/meshdim000',num2str(pn)]; % meshdim = [nPnt nEdg nTra nPar]
    %tempMd2  = ['../data/meshdim00',num2str(pn)]; % meshdim = [nPnt nEdg nTra nPar]
    meshdim12 = dlmread([tempMd2 '.dat']);
    
    nP = meshdim12(1);                    % number of points
    nE = meshdim12(2);                    % number of edges
    nT = meshdim12(3);                    % number of triangles
    nB = meshdim12(4);                    % number of global boundary nodes

    %%: Input PCE Terms: nPCEin
    ndim = 3;                    %% Number of Underalaying RVs 
    KLEord = 2;                  %% Order of Input PCE
    nPCEin = factorial(KLEord+ndim)./(factorial(KLEord).*factorial(ndim));

    %%: Output PCE Terms: nPCEout: KLEord+1 and ndim same as input PCE  
    nord = 3;                %% Maximum order of solution PCE
    nPCEout = factorial(nord+ndim)./(factorial(nord).*factorial(ndim));
    
    %%: Based on input/output PCE : Calculates Cijk : Tripple Prod
    [ijk, cijk]  = moment_ND_kij(ndim, nord, KLEord);
    
    %%: NonZero Cijk Index calculation
    % indexi = 1;
    
    %%: Stochastic Matrix Assebly
    % As = zeros(nd*nPCEout,nd*nPCEout);
    
    for k=1:nPCEin
        
        %%: Random Process Case    
        A = AssembleMatrix(p, e, t, 'PoissonModi', [], 0, k);  
        
        [Aii, Agi, Agg] = oneLevelDDMat(A, nP, nB);
        
        NAgi = [];
        for ii = 1:length(Agi(1,:))
            for jj=1:length(Agi(:,1))
                NAgi = [NAgi; Agi(jj,ii)];
            end
        end
        
        %%: Stochastic Matrix Assembly Initiated
        % Kijk = zeros(nd*nPCEout,nd*nPCEout);
        % for i=1:nPCEout
        %    for j=1:nPCEout
                
        %        if (k == ijk(indexi,1)) && (i == ijk(indexi,2)) && (j == ijk(indexi,3))
                    
                    %%: Replace it to directly populate As(i,j) = As(i,j)+ Cijk*A
        %            Kijk( ((i-1)*nd+1):(i*nd), ((j-1)*nd+1):(j*nd) ) = cijk(indexi) * A;
                    
        if (k < 10)
                    createD = horzcat('../data/Amats/subdom000',num2str(pn));
                    mkdir (createD)
                    AD  = [createD,'/AD',num2str(k)];
                    dlmwrite([AD '.dat'], A, '\t');
                    ADii  = [createD,'/ADii000',num2str(k)];
                    dlmwrite([ADii '.dat'], Aii, '\t');
                    ADgi  = [createD,'/ADgi000',num2str(k)];
                    dlmwrite([ADgi '.dat'], NAgi, '\t');
                    ADgg  = [createD,'/ADgg000',num2str(k)];
                    dlmwrite([ADgg '.dat'], Agg, '\t');
        else
                    createD = horzcat('../data/Amats/subdom000',num2str(pn));
                    mkdir (createD)
                    AD  = [createD,'/AD',num2str(k)];
                    dlmwrite([AD '.dat'], A, '\t');
                    ADii  = [createD,'/ADii00',num2str(k)];
                    dlmwrite([ADii '.dat'], Aii, '\t');
                    ADgi  = [createD,'/ADgi00',num2str(k)];
                    dlmwrite([ADgi '.dat'], NAgi, '\t');
                    ADgg  = [createD,'/ADgg00',num2str(k)];
                    dlmwrite([ADgg '.dat'], Agg, '\t');
        end
            
        %            indexi = indexi+1;
        %        end
        
        %    end
        % end
        
        %%: Stochastic Matrix Final: use spy(As) for plotting
        % As = As + Kijk;
    end
    
                      
    %%: RHS Assembly for Source Vector Initiated
    % fs = zeros(nd*nPCEout,1);
    
    %%: Deterministic Vector Assebly (f=1)
    % b = AssembleVector(p, e, t, 'PoissonModi', [], 0, 0);
    
    %Bvecs  = ['../data/Bvecs/Bvecs00',num2str(pn)];
    %dlmwrite([Bvecs '.dat'], b, '\t');
    
    %%: Stochastic Vector Assebly
    %fs(1:nd) = b;
    
    %%: Solving Intrusive System of Equations
    %u_intru = As\fs;
    
    %%: Sorting PCE Coefficient of Solution Process
    %     u = [];
    %     for i = 1:nPCEin
    %         u0 = u_intru( ((i-1)*nd)+1:i*nd );
    %         u = [u u0];
    %         %     figure(i); clf
    %         %     pdesurf(p,t,u0)
    %         %     shading faceted
    %         %     title('Computed solution')
    %     end
    %
    %     %%: Writting Solution Coefficient to a .dat File
    %     str = strcat('Dim',int2str(ndim),'_p',int2str(nord),'.dat');
    %     dlmwrite(str, u, '\t');
    %     %%: dlmread('str','');
    
    
    %%% For MALLOC CALCULATION WE CAN USE SOMETHING LIKE THIS
    %%% BUT THIS NEED PCE basis for PCEOUT not PCEin
%     for k=1:nPCEout
%         
%         %%: Random Process Case
%         A = AssembleMatrix(p, e, t, 'PoissonModi', [], 0, k);
%         
%         [Aii, Agi, Agg] = oneLevelDDMat(A, nP, nB);
%         
%         NAgi = [];
%         for ii = 1:length(Agi(1,:))
%             for jj=1:length(Agi(:,1))
%                 NAgi = [NAgi; Agi(jj,ii)];
%             end
%         end
%         
%         if (k < 10)
%             createD = horzcat('../data/Malloc/subdom000',num2str(pn));
%             mkdir (createD)
%             ADii  = [createD,'/ADii000',num2str(k)];
%             dlmwrite([ADii '.dat'], Aii, '\t');
%             ADgi  = [createD,'/ADgi000',num2str(k)];
%             dlmwrite([ADgi '.dat'], NAgi, '\t');
%             ADgg  = [createD,'/ADgg000',num2str(k)];
%             dlmwrite([ADgg '.dat'], Agg, '\t');
%         else
%             createD = horzcat('../data/Malloc/subdom000',num2str(pn));
%             mkdir (createD)
%             ADii  = [createD,'/ADii00',num2str(k)];
%             dlmwrite([ADii '.dat'], Aii, '\t');
%             ADgi  = [createD,'/ADgi00',num2str(k)];
%             dlmwrite([ADgi '.dat'], NAgi, '\t');
%             ADgg  = [createD,'/ADgg00',num2str(k)];
%             dlmwrite([ADgg '.dat'], Agg, '\t');
%         end
%     end
    
    
end
toc

exit

% 
% %%: To plot PCE coefficient
% figure(2)
% for i=1:nPCEin
%     subplot(1,nPCEin,i)
%     pdesurf(p,t,u(:,i))
% end

%%-------------------------------------------------------------------------
%%: Plot first few coefficients
%
% err1 = U1 - u1';
% enorm = max(abs(err1));
% disp(['Maximum norm error: ' num2str(enorm)])
% figure(8); clf
% pdesurf(p,t,err1)
% shading faceted
% title('Error')
%
% %err = zeros(nord,npce);
% for i = 1:3
% res = u(end:end,i)-u(1:end,i);
% error(:,i) = res/u(end:end,i);
% end
% 
% err = abs(error);
% figure(1); hold on;
% plot(err(:,1),'r')
% figure(2); hold on;
% plot(err(:,2),'r')
% figure(3); hold on;
% plot(err(:,3),'r')
% 
% figure(4); semilogy(err(:,1),'r')
% figure(5); semilogy(err(:,2),'r')
% figure(6); semilogy(err(:,3),'r')
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

