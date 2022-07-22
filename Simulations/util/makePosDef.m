function [G,err] = makePosDef(G)
% pass a covariance matrix of any size and this function attempts to make
% it positive definite, if we can't then the err is non-zero. 
% there may be better ways to ensure this.

    cnt = 1;
    scalingFac = max(max(triu(abs(nonzeros(G)),1)));
    while any(real(eig(G)) < 0) 
        G = G + scalingFac*eye(length(G)); %Changing weight on diagonal
%      G = G.*(1.001*eye(length(G)));
    cnt = cnt +1;
        if cnt == 1e6
            break
        end
    end
    [~,err] = chol(G);