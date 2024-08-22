function[f,g] = MW4(x)
if nargin == 0
	prob.nx = 10;
	prob.nf = 3;
	prob.ng = 1;
    prob.bounds(1:prob.nx,1) = zeros(prob.nx,1);
    prob.bounds(1:prob.nx,2) = ones(prob.nx,1);
	f = prob;
else
	[f,g] = MW4_true(x);
end
return


function [PopObj,PopCon] = MW4_true(x)
D = size(x,2); M = 3;
g      = sum(1 - exp(-10*((x(:,M:end).^(D-M)) - 0.5 - ((M:D) - 1)/(2*D)).^2),2);

% Objective functions
PopObj = repmat(1+g,1,M).*flip(cumprod([ones(size(x,1),1),x(:,1:M-1)],2),2).*[ones(size(x,1),1),1-x(:,M-1:-1:1)];

% Constraint functions : PopCon(x) = G(x) <= 0 where <= is feasible
l      = PopObj(:,end) - sum(PopObj(:,1:(end-1)),2);
PopCon = sum(PopObj,2) - (1 + 0.4*sin(2.5*pi*l).^8);
return