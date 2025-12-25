function [p, e, t] = getPET(pn)
%%% To read the Points elements and triangles for each subdomain           Aprl 2,2014/Ajit 
%%% "pn" is subdomain number (partition number)

%%-    Date:       Author:   Comments/Modifications:     
%   Aprl/02/2014     AD       Original

tempPt2  = ['../rpoints000',num2str(pn)];
tempRead1 = dlmread([tempPt2 '.dat']);
p = tempRead1';

%%%--------------------------------------------------------------------
%%% Exception handling
tempEg2 = ['../redges000',num2str(pn)];
test2 = load([tempEg2 '.dat']);
test3 = isempty(test2);

if test3 == 1
    % disp('floating sub-domain');
    e = [];
%%%--------------------------------------------------------------------
else
    tempRead2 = dlmread([tempEg2 '.dat']);
    e = tempRead2';
end

tempT2  = ['../rtriangles000',num2str(pn)];
tempRead3 = dlmread([tempT2 '.dat']);
t = tempRead3';

end

%%%%%%% END %%%%%%%