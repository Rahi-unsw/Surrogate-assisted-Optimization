function[f,g] = CF3(x)
if nargin == 0
	prob.nx = 10;
	prob.nf = 2;
	prob.ng = 1;
    prob.bounds(1:prob.nx,1) = [0;zeros(prob.nx-1,1)-2];
    prob.bounds(1:prob.nx,2) = [1;ones(prob.nx-1,1)+2];
	f = prob;
else
	[f,g] = CF3_true(x);
end
return


function [PopObj,PopCon] = CF3_true(X)
prob = feval('CF3');
X  = max(min(X,repmat(prob.bounds(:,2)',size(X,1),1)),repmat(prob.bounds(:,1)',size(X,1),1));
D = size(X,2); M = 2;
J1 = 3 : 2 : D;
J2 = 2 : 2 : D;
Y  = X - sin(6*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
PopObj(:,1) = X(:,1)      + 2/length(J1)*(4*sum(Y(:,J1).^2,2)-2*prod(cos(20*Y(:,J1)*pi./sqrt(repmat(J1,size(X,1),1))),2)+2);
PopObj(:,2) = 1-X(:,1).^2 + 2/length(J2)*(4*sum(Y(:,J2).^2,2)-2*prod(cos(20*Y(:,J2)*pi./sqrt(repmat(J2,size(X,1),1))),2)+2);
PopCon      = 1 - PopObj(:,2) - PopObj(:,1).^2 + sin(2*pi*(PopObj(:,1).^2-PopObj(:,2)+1));
return