function[f,g] = MW5(x)
if nargin == 0
	prob.nx = 10;
	prob.nf = 2;
	prob.ng = 3;
    prob.bounds(1:prob.nx,1) = zeros(prob.nx,1);
    prob.bounds(1:prob.nx,2) = ones(prob.nx,1);
	f = prob;
else
	[f,g] = MW5_true(x);
end
return


function [PopObj,PopCon] = MW5_true(x)
D = size(x,2); M = 2;
g = 1 + sum(1 - exp(-10*((x(:,M:end).^(D-M)) - 0.5 - (repmat(M:D,size(x,1),1) - 1)/(2*D)).^2),2);

% Objective functions
PopObj(:,1) = g.*x(:,1);
PopObj(:,2) = g.*sqrt(1 - (PopObj(:,1)./g).^2);

% Constraint functions : PopCon(x) = G(x) <= 0 where <= is feasible
l1          = atan(PopObj(:,2)./PopObj(:,1));
l2          = 0.5*pi - 2*abs(l1-0.25*pi);
PopCon(:,1) = PopObj(:,1).^2 + PopObj(:,2).^2 - (1.7-0.2*sin(2*l1)).^2;
PopCon(:,2) = (1+0.5*sin(6*l2.^3)).^2 - PopObj(:,1).^2 - PopObj(:,2).^2;
PopCon(:,3) = (1-0.45*sin(6*l2.^3)).^2 - PopObj(:,1).^2 - PopObj(:,2).^2;
return