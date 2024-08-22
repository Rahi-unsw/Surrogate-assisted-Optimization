function[f,g] = MW1(x)
if nargin == 0
	prob.nx = 10;
	prob.nf = 2;
	prob.ng = 1;
    prob.bounds(1:prob.nx,1) = zeros(prob.nx,1);
    prob.bounds(1:prob.nx,2) = ones(prob.nx,1);
	f = prob;
else
	[f,g] = MW1_true(x);
end
return


function [PopObj,PopCon] = MW1_true(x)
D = size(x,2); M = 2;
g = 1 + sum(1 - exp(-10*((x(:,M:end).^(D-M)) - 0.5 - (repmat(M:D,size(x,1),1) - 1)/(2*D)).^2),2);

% Objective functions
PopObj(:,1) = x(:,1);
PopObj(:,2) = g.*(1 - 0.85*PopObj(:,1)./g);

% Constraint functions : PopCon(x) = G(x) <= 0 where <= is feasible
l           = sqrt(2)*PopObj(:,2) - sqrt(2)*PopObj(:,1);
PopCon      = sum(PopObj,2) - 1 - 0.5*sin(2*pi*l).^8;
return