function[f,g] = CF1(x)
if nargin == 0
	prob.nx = 10;
	prob.nf = 2;
	prob.ng = 1;
    prob.bounds(1:prob.nx,1) = zeros(prob.nx,1);
    prob.bounds(1:prob.nx,2) = ones(prob.nx,1);
	f = prob;
else
	[f,g] = CF1_true(x);
end
return


function [PopObj,PopCon] = CF1_true(X)
prob = feval('CF1');
D = size(X,2); M = 2;
X  = max(min(X,repmat(prob.bounds(:,2)',size(X,1),1)),repmat(prob.bounds(:,1)',size(X,1),1));
J1 = 3 : 2 : D;
J2 = 2 : 2 : D;
PopObj(:,1) = X(:,1)   + 2*mean((X(:,J1)-repmat(X(:,1),1,length(J1)).^(0.5*(1+3*(repmat(J1,size(X,1),1)-2)/(D-2)))).^2,2);
PopObj(:,2) = 1-X(:,1) + 2*mean((X(:,J2)-repmat(X(:,1),1,length(J2)).^(0.5*(1+3*(repmat(J2,size(X,1),1)-2)/(D-2)))).^2,2);
PopCon      = 1 - PopObj(:,1) - PopObj(:,2) + abs(sin(10*pi*(PopObj(:,1)-PopObj(:,2)+1)));
return