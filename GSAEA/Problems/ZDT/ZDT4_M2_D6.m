function [f,g] = ZDT4_M2_D6(x)
if nargin == 0
    prob.nf = 2;
    prob.ng = 0;
    prob.nx = 3*prob.nf;
    prob.cast = [];
    prob.bounds = [[0,zeros(1,prob.nx-1)-5]',[1,zeros(1,prob.nx-1)+5]'];
    f = prob;
    g = [];
else
    [f,g] = ZDT4_2d_true(x);
end
return

function [f,gg] = ZDT4_2d_true(PopDec)
f(:,1) = PopDec(:,1);
g = 1 + 10*(size(PopDec,2)-1) + sum(PopDec(:,2:end).^2-10*cos(4*pi*PopDec(:,2:end)),2);
h = 1 - (f(:,1)./g).^0.5;
f(:,2) = g.*h;
gg = [];
return
