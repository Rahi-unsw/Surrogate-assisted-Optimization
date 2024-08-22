function [f,g] = TR3(x)
if nargin == 0
    prob.nf = 2;
    prob.ng = 0;
    prob.nx = 3;
    prob.cast = [];
    for i = 1:prob.nx
        prob.bounds(i,:) = [-5,5];
    end
    f = prob;
    g = [];
else
    [f,g] = TR3_true(x);
end
return

function [f,g] = TR3_true(x)
g = x(:,2).^2 + x(:,3).^2;
f(:,1)=-x(:,1) + g ;
f(:,2)=tanh(x(:,1)) + g;
g = [];
return
