function outputY = getub(inputX, m1, nB, nBG)   %%(pid, nb, nBG, inputX)
% outputY = zeros(1, nBG);

%tempNodes2  = ['rbnodes00',num2str(m1)];  % Where, i is the domain number
%tempNodes2  = ['bnodes00',num2str(m1)];  % Where, i is the domain number

tempNodes2  = ['../rbnodes',pad(num2str(m1),4,'left','0')];
%%% NOTE: we need to use re-arranged (nbones[r,c]) interface nodes  
%%% because, we also re-arrange nodes -> Nnodes[i,r,c] which is 
%%% later used for re-assembling final solution vector
bnodes = dlmread([tempNodes2 '.dat']);

for i = 1:nB
    outputY(i) = inputX(bnodes(i));
end
outputY = outputY';

end

%%%%%%% END %%%%%%%