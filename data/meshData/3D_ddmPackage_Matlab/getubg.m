function outputY = getubg(inputX, m1, nB, nBG)   %%(pid, nb, nBG, inputX)
outputY = zeros(1, nBG);

%%% NOTE: we need to use re-arranged (nbones[r,c]) interface nodes  
%%% because, we also re-arrange nodes -> Nnodes[i,r,c] which is 
%%% later used for re-assembling final solution vector
%tempNodes2  = ['nbnodes00',num2str(m1)]; % Where, m1 is the domain number
tempNodes2  = ['../rbnodes',pad(num2str(m1),4,'left','0')];
bnodes = dlmread([tempNodes2 '.dat']);

for i = 1:nB                             % Where, nB is number of local boundary nodes
    outputY(bnodes(i)) = inputX(i);
end
outputY = outputY';

end

%%%%%%% END %%%%%%%