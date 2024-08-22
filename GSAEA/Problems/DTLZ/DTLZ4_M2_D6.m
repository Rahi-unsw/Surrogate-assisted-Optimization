function [f,g] = DTLZ4_M2_D6(x)
if nargin == 0
    prob.nf = 2;
    prob.ng = 0;
    prob.nx = 3*prob.nf;
    prob.cast = [];
    for i = 1:prob.nx
        prob.bounds(i,:) = [0,1];
    end
    f = prob;
    g = [];
else
    [f,g] = DTLZ4_2d_true(x);
end
return

function [f,gg] = DTLZ4_2d_true(PopDec)
M      = 2;
PopDec(:,1:M-1) = PopDec(:,1:M-1).^100;
g      = sum((PopDec(:,M:end)-0.5).^2,2);
f = repmat(1+g,1,M).*fliplr(cumprod([ones(size(g,1),1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(PopDec(:,M-1:-1:1)*pi/2)];
gg      = [];
return
