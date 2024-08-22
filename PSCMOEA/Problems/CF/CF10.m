function[f,g] = CF10(x)
if nargin == 0
	prob.nx = 10;
	prob.nf = 3;
	prob.ng = 1;
    prob.bounds(1:prob.nx,1) = [0;0;zeros(prob.nx-2,1)-2];
    prob.bounds(1:prob.nx,2) = [1;1;ones(prob.nx-2,1)+2];
	f = prob;
else
	[f,g] = CF10_true(x);
end
return


function [PopObj,PopCon] = CF10_true(X)
prob = feval('CF10');
X  = max(min(X,repmat(prob.bounds(:,2)',size(X,1),1)),repmat(prob.bounds(:,1)',size(X,1),1));
D  = size(X,2);
J1 = 4 : 3 : D;
J2 = 5 : 3 : D;
J3 = 3 : 3 : D;
Y  = X - 2*repmat(X(:,2),1,D).*sin(2*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
PopObj(:,1) = cos(0.5*X(:,1)*pi).*cos(0.5*X(:,2)*pi) + 2*mean(4*Y(:,J1).^2-cos(8*pi*Y(:,J1))+1,2);
PopObj(:,2) = cos(0.5*X(:,1)*pi).*sin(0.5*X(:,2)*pi) + 2*mean(4*Y(:,J2).^2-cos(8*pi*Y(:,J2))+1,2);
PopObj(:,3) = sin(0.5*X(:,1)*pi)                     + 2*mean(4*Y(:,J3).^2-cos(8*pi*Y(:,J3))+1,2);
PopCon      = 1 - (PopObj(:,1).^2+PopObj(:,2).^2)./(1-PopObj(:,3).^2) + sin(2*pi*((PopObj(:,1).^2-PopObj(:,2).^2)./(1-PopObj(:,3).^2)+1));
return