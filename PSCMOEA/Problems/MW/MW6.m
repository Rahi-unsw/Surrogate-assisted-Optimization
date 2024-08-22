function[f,g] = MW6(x)
if nargin == 0
	prob.nx = 10;
	prob.nf = 2;
	prob.ng = 1;
    prob.bounds(1:prob.nx,1) = zeros(prob.nx,1);
    prob.bounds(1:prob.nx,2) = ones(prob.nx,1);
	f = prob;
else
	[f,g] = MW6_true(x);
end
return


function [PopObj,PopCon] = MW6_true(x)
D = size(x,2); M = 2;
z = 1 - exp(-10*(x(:,M:end) - (repmat(M:D,size(x,1),1) - 1)/D).^2);
g = 1 + sum((1.5 + (0.1/D)*z.^2 - 1.5*cos(2*pi*z)),2);

% Objective functions
PopObj(:,1) = g.*x(:,1)*1.0999;
PopObj(:,2) = g.*sqrt(1.1*1.1 - (PopObj(:,1)./g).^2);

% Constraint functions : PopCon(x) = G(x) <= 0 where <= is feasible
l           = cos(6*atan(PopObj(:,2)./PopObj(:,1)).^4).^10;
PopCon      = (PopObj(:,1)./(1+0.15*l)).^2 + (PopObj(:,2)./(1+0.75*l)).^2 - 1;
return