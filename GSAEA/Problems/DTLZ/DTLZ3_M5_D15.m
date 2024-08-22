function [f,g] = DTLZ3_M5_D15(x)
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
    [f,g] = DTLZ3_2d_true(x);
end
return

function [f,g] = DTLZ3_2d_true(PopDec)
[N,D]  = size(PopDec);
M      = 5;
gg     = (D-M+1+sum((PopDec(:,M:end)-0.5).^2-cos(20.*pi.*(PopDec(:,M:end)-0.5)),2));%100* (Removed in the modified form)
f      = repmat(1+gg,1,M).*fliplr(cumprod([ones(N,1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(N,1),sin(PopDec(:,M-1:-1:1)*pi/2)];
g      = [];
return
