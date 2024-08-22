function[f,g] = MW13(x)
if nargin == 0
	prob.nx = 10;
	prob.nf = 2;
	prob.ng = 2;
    prob.bounds(1:prob.nx,1) = zeros(prob.nx,1);
    prob.bounds(1:prob.nx,2) = ones(prob.nx,1);
	f = prob;
else
	[f,g] = MW13_true(x);
end
return


function [PopObj,PopCon] = MW13_true(x)
D = size(x,2); M = 2;
z = 1 - exp(-10*(x(:,M:end) - (repmat(M:D,size(x,1),1) - 1)/D).^2);
g = 1 + sum((1.5 + (0.1/D)*z.^2 - 1.5*cos(2*pi*z)),2);

% Objective functions
PopObj(:,1) = g.*x(:,1)*1.5;
PopObj(:,2) = g.*(5 - exp(PopObj(:,1)./g) - abs(0.5*sin(3*pi*PopObj(:,1)./g)));

% Constraint functions : PopCon(x) = G(x) <= 0 where <= is feasible
PopCon(:,1) = (5 - exp(PopObj(:,1)) - 0.5*sin(3*pi*PopObj(:,1)) - PopObj(:,2)).*(5 - (1 + 0.4*PopObj(:,1)) - 0.5*sin(3*pi*PopObj(:,1)) - PopObj(:,2));
PopCon(:,2) = -(5 - (1 + PopObj(:,1) + 0.5*PopObj(:,1).^2) - 0.5*sin(3*pi*PopObj(:,1)) - PopObj(:,2)).*(5 - (1 + 0.7*PopObj(:,1)) - 0.5*sin(3*pi*PopObj(:,1)) - PopObj(:,2));
return