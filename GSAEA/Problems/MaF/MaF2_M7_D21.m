function [f,g] = MaF2_M7_D21(x)
if nargin == 0
    prob.nf = 7;
    prob.ng = 0;
    prob.nx = 3*prob.nf;
    prob.cast = [];
    for i = 1:prob.nx
        prob.bounds(i,:) = [0,1];
    end
    f = prob;
    g = [];
else
    [f,g] = MaF2_2d_true(x);
end
return

function [f,gg] = MaF2_2d_true(PopDec)
[N,D] = size(PopDec);
M     = 7;
g     = zeros(N,M);
for m = 1 : M
    if m < M
        g(:,m) = sum(((PopDec(:,M+(m-1)*floor((D-M+1)/M):M+m*floor((D-M+1)/M)-1)/2+1/4)-0.5).^2,2);
    else
        g(:,m) = sum(((PopDec(:,M+(M-1)*floor((D-M+1)/M):D)/2+1/4)-0.5).^2,2);
    end
end
f = (1+g).*fliplr(cumprod([ones(size(g,1),1),cos((PopDec(:,1:M-1)/2+1/4)*pi/2)],2)).*[ones(size(g,1),1),sin((PopDec(:,M-1:-1:1)/2+1/4)*pi/2)];
gg = [];
return