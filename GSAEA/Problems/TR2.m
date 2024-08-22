function [f,g] = TR2(x)
if nargin == 0
    prob.nf = 2;
    prob.ng = 0;
    prob.nx = 2;
    prob.cast = [];
    for i = 1:prob.nx
        prob.bounds(i,:) = [0,1];
    end
    f = prob;
    g = [];
else
    [f,g] = TR2_true(x);
end
return

function [f,g] = TR2_true(x)
f(:,1) = x(:,1);
L = (1-x(:,1).^6).^(1/6);
U = 1;
f(:,2) =0.5*(L+x(:,2).*U);
g = [];
return
