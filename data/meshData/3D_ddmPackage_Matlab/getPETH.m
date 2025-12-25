function [p, e, t, h] = getPETH(pn)
%%% To read the Points elements and triangles for each subdomain           Aprl 2,2014/Ajit 
%%% "pn" is subdomain number (partition number)

%%-    Date:       Author:   Comments/Modifications:     
%   Aprl/02/2014     AD       Original

nZ = 4;
%tempPt2  = ['../data/rpoints00',num2str(pn)];
tempPt2  = ['../rpoints',pad(num2str(pn),nZ,'left','0')];
tempRead1 = dlmread([tempPt2 '.dat']);
p = tempRead1';

%%%--------------------------------------------------------------------
%%% Exception handling
%tempEg2 = ['../data/redges00',num2str(pn)];
tempEg2 = ['../redges',pad(num2str(pn),nZ,'left','0')];
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

%tempT2  = ['../data/rtriangles00',num2str(pn)];
tempT2  = ['../rtriangles',pad(num2str(pn),nZ,'left','0')];
tempRead3 = dlmread([tempT2 '.dat']);
t = tempRead3';

%tempH2  = ['../data/rtetrahedrons00',num2str(pn)];
tempH2  = ['../rtetrahedron',pad(num2str(pn),nZ,'left','0')];
tempRead4 = dlmread([tempH2 '.dat']);
h = tempRead4';

end

%%%%%%% END %%%%%%%