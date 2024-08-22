function [f,g] = MaF1_M7_D21(x)
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
    [f,g] = MaF1_2d_true(x);
end
return

function [f,g] = MaF1_2d_true(PopDec)
M      = 7;
g      = sum((PopDec(:,M:end)-0.5).^2,2);
f = repmat(1+g,1,M) - repmat(1+g,1,M).*fliplr(cumprod([ones(size(g,1),1),PopDec(:,1:M-1)],2)).*[ones(size(g,1),1),1-PopDec(:,M-1:-1:1)];
g = [];
return
