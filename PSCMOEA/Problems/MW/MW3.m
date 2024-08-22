function[f,g] = MW3(x)
if nargin == 0
	prob.nx = 10;
	prob.nf = 2;
	prob.ng = 2;
    prob.bounds(1:prob.nx,1) = zeros(prob.nx,1);
    prob.bounds(1:prob.nx,2) = ones(prob.nx,1);
	f = prob;
else
	[f,g] = MW3_true(x);
end
return


function [PopObj,PopCon] = MW3_true(x)
D = size(x,2); M = 2;
g = 1 + sum(2*(x(:,M:end) + (x(:,M-1:end-1) - 0.5).^2 - 1).^2,2);

% Objective functions
PopObj(:,1) = x(:,1);
PopObj(:,2) = g.*(1 - PopObj(:,1)./g);

% Constraint functions : PopCon(x) = G(x) <= 0 where <= is feasible
l           = sqrt(2)*PopObj(:,2) - sqrt(2)*PopObj(:,1);
PopCon(:,1) = sum(PopObj,2) - 1.05 - 0.45*sin(0.75*pi*l).^6;
PopCon(:,2) = 0.85 - sum(PopObj,2) + 0.3*sin(0.75*pi*l).^2;
return