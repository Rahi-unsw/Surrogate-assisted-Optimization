function[f,g] = MW11(x)
if nargin == 0
	prob.nx = 10;
	prob.nf = 2;
	prob.ng = 4;
    prob.bounds(1:prob.nx,1) = zeros(prob.nx,1);
    prob.bounds(1:prob.nx,2) = ones(prob.nx,1);
	f = prob;
else
	[f,g] = MW11_true(x);
end
return


function [PopObj,PopCon] = MW11_true(x)
D = size(x,2); M = 2;
g = 1 + sum(2*(x(:,M:end) + (x(:,M-1:end-1) - 0.5).^2 - 1).^2,2);

% Objective functions
PopObj(:,1) = g.*x(:,1)*sqrt(1.9999);
PopObj(:,2) = g.*sqrt(2 - (PopObj(:,1)./g).^2);

% Constraint functions : PopCon(x) = G(x) <= 0 where <= is feasible
PopCon(:,1) = -(3 - PopObj(:,1).^2 - PopObj(:,2)).*(3 - 2*PopObj(:,1).^2 - PopObj(:,2));
PopCon(:,2) = (3 - 0.625*PopObj(:,1).^2 - PopObj(:,2)).*(3 - 7*PopObj(:,1).^2 - PopObj(:,2));
PopCon(:,3) = -(1.62 - 0.18*PopObj(:,1).^2 - PopObj(:,2)).*(1.125 - 0.125*PopObj(:,1).^2 - PopObj(:,2));
PopCon(:,4) = (2.07 - 0.23*PopObj(:,1).^2 - PopObj(:,2)).*(0.63 - 0.07*PopObj(:,1).^2 - PopObj(:,2));
return