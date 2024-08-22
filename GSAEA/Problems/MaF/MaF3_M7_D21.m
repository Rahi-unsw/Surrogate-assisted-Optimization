function [f,g] = MaF3_M7_D21(x)
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
    [f,g] = MaF3_2d_true(x);
end
return

function [f,gg] = MaF3_2d_true(PopDec)
M      = 7;
g      = (size(PopDec,2)-M+1+sum((PopDec(:,M:end)-0.5).^2-cos(20.*pi.*(PopDec(:,M:end)-0.5)),2));
f = repmat(1+g,1,M).*fliplr(cumprod([ones(size(g,1),1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(PopDec(:,M-1:-1:1)*pi/2)];
f = [f(:,1:M-1).^4,f(:,M).^2];
gg = [];
return