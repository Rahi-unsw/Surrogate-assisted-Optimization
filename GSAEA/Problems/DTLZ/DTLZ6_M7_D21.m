function [f,g] = DTLZ6_M7_D21(x)
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
    [f,g] = DTLZ6_2d_true(x);
end
return

function [f,gg] = DTLZ6_2d_true(PopDec)
M      = 7;
g      = sum(PopDec(:,M:end).^0.1,2);
Temp   = repmat(g,1,M-2);
PopDec(:,2:M-1) = (1+2*Temp.*PopDec(:,2:M-1))./(2+2*Temp);
f = repmat(1+g,1,M).*fliplr(cumprod([ones(size(g,1),1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(PopDec(:,M-1:-1:1)*pi/2)];
gg      = [];
return
