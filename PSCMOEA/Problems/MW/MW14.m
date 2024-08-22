function[f,g] = MW14(x)
if nargin == 0
	prob.nx = 10;
	prob.nf = 3;
	prob.ng = 1;
    prob.bounds(1:prob.nx,1) = zeros(prob.nx,1);
    prob.bounds(1:prob.nx,2) = ones(prob.nx,1);
	f = prob;
else
	[f,g] = MW14_true(x);
end
return


function [PopObj,PopCon] = MW14_true(x)
D = size(x,2); M = 3;
x      = 1.5*x;
g      = sum(2*(x(:,M:end) + (x(:,M-1:end-1) - 0.5).^2 - 1).^2,2);

% Objective functions
PopObj(:,1:M-1) = x(:,1:M-1);
PopObj(:,M)     = ((1+g)/(M-1)).*sum(6 - exp(PopObj(:,1:M-1)) - 1.5*sin(1.1*pi*PopObj(:,1:M-1).^2),2);

% Constraint functions : PopCon(x) = G(x) <= 0 where <= is feasible
a           = 1 + PopObj(:,1:M-1) + 0.5*PopObj(:,1:M-1).^2 + 1.5*sin(1.1*pi*PopObj(:,1:M-1).^2);
PopCon      = PopObj(:,M) - 1/(M-1)*sum(6.1 - a,2);
return