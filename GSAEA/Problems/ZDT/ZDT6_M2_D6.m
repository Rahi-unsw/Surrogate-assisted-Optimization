function [f,g] = ZDT6_M2_D6(x)
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
    [f,g] = ZDT6_2d_true(x);
end
return

function [f,gg] = ZDT6_2d_true(PopDec)
f(:,1) = 1 - exp(-4*PopDec(:,1)).*sin(6*pi*PopDec(:,1)).^6;
g = 1 + 9*mean(PopDec(:,2:end),2).^0.25;
h = 1 - (f(:,1)./g).^2;
f(:,2) = g.*h;
gg = [];
return
