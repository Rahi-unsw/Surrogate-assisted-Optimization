function[f,g] = MW8(x)
if nargin == 0
	prob.nx = 10;
	prob.nf = 3;
	prob.ng = 1;
    prob.bounds(1:prob.nx,1) = zeros(prob.nx,1);
    prob.bounds(1:prob.nx,2) = ones(prob.nx,1);
	f = prob;
else
	[f,g] = MW8_true(x);
end
return


function [PopObj,PopCon] = MW8_true(x)
D = size(x,2); M = 3;
z = 1 - exp(-10*(x(:,M:end) - ((M:D) - 1)/D).^2);
g = sum((1.5 + (0.1/D)*z.^2 - 1.5*cos(2*pi*z)),2);

% Objective functions
PopObj = repmat(1+g,1,M).*flip(cumprod([ones(size(x,1),1),cos(x(:,1:M-1)*pi/2)],2),2).*[ones(size(x,1),1),sin(x(:,M-1:-1:1)*pi/2)];

% Constraint functions : PopCon(x) = G(x) <= 0 where <= is feasible
l      = asin(PopObj(:,end)./sqrt(sum(PopObj.^2,2)));
PopCon = sum(PopObj.^2,2) - (1.25 - 0.5*sin(6*l).^2).^2;
return