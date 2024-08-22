function [f,g] = MaF7_M2_D6(x)
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
    [f,g] = MaF7_2d_true(x);
end
return

function [PopObj,gg] = MaF7_2d_true(PopDec)
M               = 2;
PopObj          = zeros(size(PopDec,1),M);
g               = 1+9*mean(PopDec(:,M:end),2);
PopObj(:,1:M-1) = PopDec(:,1:M-1);
PopObj(:,M)     = (1+g).*(M-sum(PopObj(:,1:M-1)./(1+repmat(g,1,M-1)).*(1+sin(3*pi.*PopObj(:,1:M-1))),2));
gg = [];
return