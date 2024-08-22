function [f,g] = DTLZ2_M5_D15(x)
if nargin == 0
    prob.nf = 5;
    prob.ng = 0;
    prob.nx = 3*prob.nf;
    prob.cast = [];
    for i = 1:prob.nx
        prob.bounds(i,:) = [0,1];
    end
    f = prob;
    g = [];
else
    [f,g] = DTLZ2_2d_true(x);
end
return

function [f,g] = DTLZ2_2d_true(PopDec)
M      = 5;
gg     = sum((PopDec(:,M:end)-0.5).^2,2);
f      = repmat(1+gg,1,M).*fliplr(cumprod([ones(size(gg,1),1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(size(gg,1),1),sin(PopDec(:,M-1:-1:1)*pi/2)];
g      = [];
return
