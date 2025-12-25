function [nP, nE, nT, nH, nB] = getMeshDim(pn)
% To extract mesh dimention for each partition                             Aprl 2,2014/Ajit 
% "pn" is subdomain number (partition number)

%%-    Date:       Author:   Comments/Modifications:     
%   Aprl/02/2014     AD       Original

    %tempMd2  = ['../data/meshdim00',num2str(pn)]; % meshdim = [nPnt nEdg nTra nPar]
    tempMd2  = ['../meshdim',pad(num2str(pn),4,'left','0')];
    meshdim12 = dlmread([tempMd2 '.dat']);

    nP = meshdim12(1);                    % number of points
    nE = meshdim12(2);                    % number of edges
    nT = meshdim12(3);                    % number of traingles
    nH = meshdim12(4);                    % number of tetrahedrons
    nB = meshdim12(5);                    % number of boundary nodes

end

%%%%%%% END %%%%%%%
