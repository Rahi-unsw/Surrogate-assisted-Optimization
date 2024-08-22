function [f,g] = MaF8_M3_D2(x)
if nargin == 0
    prob.nf = 3;
    prob.ng = 0;
    prob.nx = 2;
    prob.cast = [];
    for i = 1:prob.nx
        prob.bounds(i,:) = [-100,100];
    end
    f = prob;
    g = [];
else
    [f,g] = MaF8_M3_D2_true(x);
end
return

function [f,gg] = MaF8_M3_D2_true(PopDec)
% Generate vertexes
M   = 3;
Points  = [];
[thera,rho] = cart2pol(0,1);
[Points(:,1),Points(:,2)] = pol2cart(thera-(1:M)*2*pi/M,rho);
f     = pdist2(PopDec,Points);
gg = [];
return