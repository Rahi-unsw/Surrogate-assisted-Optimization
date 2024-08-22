function[f,g] = CF6(x)
if nargin == 0
	prob.nx = 10;
	prob.nf = 2;
	prob.ng = 2;
    prob.bounds(1:prob.nx,1) = [0;zeros(prob.nx-1,1)-2];
    prob.bounds(1:prob.nx,2) = [1;ones(prob.nx-1,1)+2];
	f = prob;
else
	[f,g] = CF6_true(x);
end
return


function [PopObj,PopCon] = CF6_true(X)
prob = feval('CF6');
D  = size(X,2);
J1 = 3 : 2 : D;
J2 = 2 : 2 : D;
Y       = zeros(size(X));
Y(:,J1) = X(:,J1) - 0.8*repmat(X(:,1),1,length(J1)).*cos(6*pi*repmat(X(:,1),1,length(J1))+repmat(J1,size(X,1),1)*pi/D);
Y(:,J2) = X(:,J2) - 0.8*repmat(X(:,1),1,length(J2)).*sin(6*pi*repmat(X(:,1),1,length(J2))+repmat(J2,size(X,1),1)*pi/D);
PopObj(:,1) = X(:,1)        + sum(Y(:,J1).^2,2);
PopObj(:,2) = (1-X(:,1)).^2 + sum(Y(:,J2).^2,2);
PopCon(:,1) = -X(:,2) + 0.8*X(:,1).*sin(6*pi*X(:,1)+2*pi/size(X,2)) + sign(0.5*(1-X(:,1))-(1-X(:,1)).^2).*sqrt(abs(0.5*(1-X(:,1))-(1-X(:,1)).^2));
PopCon(:,2) = -X(:,4) + 0.8*X(:,1).*sin(6*pi*X(:,1)+4*pi/size(X,2)) + sign(0.25*sqrt(1-X(:,1))-0.5*(1-X(:,1))).*sqrt(abs(0.25*sqrt(1-X(:,1))-0.5*(1-X(:,1))));
return