function [f,g] = MaF4_M3_D9(x)
if nargin == 0
    prob.nf = 3;
    prob.ng = 0;
    prob.nx = 3*prob.nf;
    prob.cast = [];
    for i = 1:prob.nx
        prob.bounds(i,:) = [0,1];
    end
    f = prob;
    g = [];
else
    [f,g] = MaF4_2d_true(x);
end
return

function [f,gg] = MaF4_2d_true(PopDec)
[N,D]  = size(PopDec);
M      = 3;
g      = (D-M+1+sum((PopDec(:,M:end)-0.5).^2-cos(20.*pi.*(PopDec(:,M:end)-0.5)),2));
PopObj = repmat(1+g,1,M) - repmat(1+g,1,M).*fliplr(cumprod([ones(N,1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(N,1),sin(PopDec(:,M-1:-1:1)*pi/2)];
f = PopObj.*repmat(2.^(1:M),N,1);
gg = [];
return