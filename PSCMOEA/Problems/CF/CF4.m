function[f,g] = CF4(x)
if nargin == 0
	prob.nx = 10;
	prob.nf = 2;
	prob.ng = 1;
    prob.bounds(1:prob.nx,1) = [0;zeros(prob.nx-1,1)-2];
    prob.bounds(1:prob.nx,2) = [1;ones(prob.nx-1,1)+2];
	f = prob;
else
	[f,g] = CF4_true(x);
end
return


function [PopObj,PopCon] = CF4_true(X)
prob = feval('CF4');
D = size(X,2); M = 2;
J1 = 3 : 2 : D;
J2 = 2 : 2 : D;
Y  = X - sin(6*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
h           = Y.^2;
temp        = Y(:,2) < 3/2*(1-sqrt(1/2));
h(temp,2)   = abs(Y(temp,2));
h(~temp,2)  = 0.125 + (Y(~temp,2)-1).^2;
PopObj(:,1) = X(:,1)   + sum(h(:,J1),2);
PopObj(:,2) = 1-X(:,1) + sum(h(:,J2),2);
t      = X(:,2) - sin(6*pi*X(:,1)+2*pi/size(X,2)) - 0.5*X(:,1) + 0.25;
PopCon = -t./(1+exp(4*abs(t)));
return