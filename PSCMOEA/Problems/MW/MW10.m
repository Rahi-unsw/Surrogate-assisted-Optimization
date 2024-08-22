function[f,g] = MW10(x)
if nargin == 0
	prob.nx = 10;
	prob.nf = 2;
	prob.ng = 3;
    prob.bounds(1:prob.nx,1) = zeros(prob.nx,1);
    prob.bounds(1:prob.nx,2) = ones(prob.nx,1);
	f = prob;
else
	[f,g] = MW10_true(x);
end
return


function [PopObj,PopCon] = MW10_true(x)
D = size(x,2); M = 2;
z = 1 - exp(-10*(x(:,M:end) - (repmat(M:D,size(x,1),1) - 1)/D).^2);
g = 1 + sum((1.5 + (0.1/D)*z.^2 - 1.5*cos(2*pi*z)),2);

% Objective functions
PopObj(:,1) = g.*(x(:,1).^D);
PopObj(:,2) = g.*(1 - (PopObj(:,1)./g).^2);

% Constraint functions : PopCon(x) = G(x) <= 0 where <= is feasible
PopCon(:,1) = -(2 - 4*PopObj(:,1).^2 - PopObj(:,2)).*(2 - 8*PopObj(:,1).^2 - PopObj(:,2));
PopCon(:,2) = (2 - 2*PopObj(:,1).^2 - PopObj(:,2)).*(2 - 16*PopObj(:,1).^2 - PopObj(:,2));
PopCon(:,3) = (1 - PopObj(:,1).^2 - PopObj(:,2)).*(1.2 - 1.2*PopObj(:,1).^2 - PopObj(:,2));
return