function[f,g] = MW7(x)
if nargin == 0
	prob.nx = 10;
	prob.nf = 2;
	prob.ng = 2;
    prob.bounds(1:prob.nx,1) = zeros(prob.nx,1);
    prob.bounds(1:prob.nx,2) = ones(prob.nx,1);
	f = prob;
else
	[f,g] = MW7_true(x);
end
return


function [PopObj,PopCon] = MW7_true(x)
D = size(x,2); M = 2;
g = 1 + sum(2*(x(:,M:end) + (x(:,M-1:end-1) - 0.5).^2 - 1).^2,2);

% Objective functions
PopObj(:,1) = g.*x(:,1);
PopObj(:,2) = g.*sqrt(1 - (PopObj(:,1)./g).^2);

% Constraint functions : PopCon(x) = G(x) <= 0 where <= is feasible
l           = atan(PopObj(:,2)./PopObj(:,1));
PopCon(:,1) = PopObj(:,1).^2 + PopObj(:,2).^2 - (1.2+0.4*sin(4*l).^16).^2;
PopCon(:,2) = (1.15-0.2*sin(4*l).^8).^2 - PopObj(:,1).^2 - PopObj(:,2).^2;
return