function[f,g] = MW12(x)
if nargin == 0
	prob.nx = 10;
	prob.nf = 2;
	prob.ng = 2;
    prob.bounds(1:prob.nx,1) = zeros(prob.nx,1);
    prob.bounds(1:prob.nx,2) = ones(prob.nx,1);
	f = prob;
else
	[f,g] = MW12_true(x);
end
return


function [PopObj,PopCon] = MW12_true(x)
D = size(x,2); M = 2;
g = 1 + sum(1 - exp(-10*((x(:,M:end).^(D-M)) - 0.5 - (repmat(M:D,size(x,1),1) - 1)/(2*D)).^2),2);

% Objective functions
PopObj(:,1) = g.*x(:,1);
PopObj(:,2) = g.*(0.85 - 0.8*(PopObj(:,1)./g) - 0.08*abs(sin(3.2*pi*(PopObj(:,1)./g))));

% Constraint functions : PopCon(x) = G(x) <= 0 where <= is feasible
PopCon(:,1) = (1 - 0.8*PopObj(:,1) - PopObj(:,2) + 0.08*sin(2*pi*(PopObj(:,2) - PopObj(:,1)/1.5))).*(1.8 - 1.125*PopObj(:,1) - PopObj(:,2) + 0.08*sin(2*pi*(PopObj(:,2)/1.8 - PopObj(:,1)/1.6)));
PopCon(:,2) = -(1 - 0.625*PopObj(:,1) - PopObj(:,2) + 0.08*sin(2*pi*(PopObj(:,2) - PopObj(:,1)/1.6))).*(1.4 - 0.875*PopObj(:,1) - PopObj(:,2) + 0.08*sin(2*pi*(PopObj(:,2)/1.4 - PopObj(:,1)/1.6)));
return