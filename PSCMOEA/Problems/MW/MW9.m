function[f,g] = MW9(x)
if nargin == 0
	prob.nx = 10;
	prob.nf = 2;
	prob.ng = 1;
    prob.bounds(1:prob.nx,1) = zeros(prob.nx,1);
    prob.bounds(1:prob.nx,2) = ones(prob.nx,1);
	f = prob;
else
	[f,g] = MW9_true(x);
end
return


function [PopObj,PopCon] = MW9_true(x)
D = size(x,2); M = 2;
g = 1 + sum(1 - exp(-10*((x(:,M:end).^(D-M)) - 0.5 - (repmat(M:D,size(x,1),1) - 1)/(2*D)).^2),2);

% Objective functions
PopObj(:,1) = g.*x(:,1);
PopObj(:,2) = g.*(1 - (PopObj(:,1)./g).^0.6);

% Constraint functions : PopCon(x) = G(x) <= 0 where <= is feasible
T1          = (1-0.64*PopObj(:,1).^2-PopObj(:,2)).*(1-0.36*PopObj(:,1).^2-PopObj(:,2));
T2          = 1.35.^2 - (PopObj(:,1)+0.35).^2 - PopObj(:,2);
T3          = 1.15.^2 - (PopObj(:,1)+0.15).^2 - PopObj(:,2);
PopCon      = min(T1,T2.*T3);
return