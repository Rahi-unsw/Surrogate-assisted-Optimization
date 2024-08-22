function[f,g] = CF9(x)
if nargin == 0
	prob.nx = 10;
	prob.nf = 3;
	prob.ng = 1;
    prob.bounds(1:prob.nx,1) = [0;0;zeros(prob.nx-2,1)-2];
    prob.bounds(1:prob.nx,2) = [1;1;ones(prob.nx-2,1)+2];
	f = prob;
else
	[f,g] = CF9_true(x);
end
return


function [PopObj,PopCon] = CF9_true(X)
prob = feval('CF9');
X  = max(min(X,repmat(prob.bounds(:,2)',size(X,1),1)),repmat(prob.bounds(:,1)',size(X,1),1));
D  = size(X,2);
J1 = 4 : 3 : D;
J2 = 5 : 3 : D;
J3 = 3 : 3 : D;
Y  = X - 2*repmat(X(:,2),1,D).*sin(2*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
PopObj(:,1) = cos(0.5*X(:,1)*pi).*cos(0.5*X(:,2)*pi) + 2*mean(Y(:,J1).^2,2);
PopObj(:,2) = cos(0.5*X(:,1)*pi).*sin(0.5*X(:,2)*pi) + 2*mean(Y(:,J2).^2,2);
PopObj(:,3) = sin(0.5*X(:,1)*pi)                     + 2*mean(Y(:,J3).^2,2);
PopCon      = 1 - (PopObj(:,1).^2+PopObj(:,2).^2)./(1-PopObj(:,3).^2) + 3*sin(2*pi*((PopObj(:,1).^2-PopObj(:,2).^2)./(1-PopObj(:,3).^2)+1));
return