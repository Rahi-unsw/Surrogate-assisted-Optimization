function[f,g] = CF2(x)
if nargin == 0
	prob.nx = 10;
	prob.nf = 2;
	prob.ng = 1;
    prob.bounds(1:prob.nx,1) = [0;zeros(prob.nx-1,1)-1];
    prob.bounds(1:prob.nx,2) = ones(prob.nx,1);
	f = prob;
else
	[f,g] = CF2_true(x);
end
return


function [PopObj,PopCon] = CF2_true(X)
prob = feval('CF2');
D = size(X,2); M = 2;
X  = max(min(X,repmat(prob.bounds(:,2)',size(X,1),1)),repmat(prob.bounds(:,1)',size(X,1),1));
J1 = 3 : 2 : D;
J2 = 2 : 2 : D;
PopObj(:,1) = X(:,1)         + 2*mean((X(:,J1)-sin(6*pi*repmat(X(:,1),1,length(J1))+repmat(J1,size(X,1),1)*pi/D)).^2,2);
PopObj(:,2) = 1-sqrt(X(:,1)) + 2*mean((X(:,J2)-cos(6*pi*repmat(X(:,1),1,length(J2))+repmat(J2,size(X,1),1)*pi/D)).^2,2);
t           = PopObj(:,2) + sqrt(PopObj(:,1)) - sin(2*pi*(sqrt(PopObj(:,1))-PopObj(:,2)+1)) - 1;
PopCon      = -t./(1+exp(4*abs(t)));
return

