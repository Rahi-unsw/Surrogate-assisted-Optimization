function[f,g] = LIRCMOP13(x)
if nargin == 0
	prob.nx = 10;
	prob.nf = 3;
	prob.ng = 2;
    prob.bounds(1:prob.nx,1) = zeros(prob.nx,1);
    prob.bounds(1:prob.nx,2) = ones(prob.nx,1);
	f = prob;
else
	[f,g] = LIRCMOP13_true(x);
end
return


function [PopObj,PopCon] = LIRCMOP13_true(X)
prob = feval('LIRCMOP13');
%% Calculate objective values
X = max(min(X,repmat(prob.bounds(:,2)',size(X,1),1)),repmat(prob.bounds(:,1)',size(X,1),1));
[popsize,variable_length] = size(X);
sum1 = zeros(popsize,1);
for j = 3 : variable_length
    sum1 = sum1+10*(X(:,j)-0.5).^2;
end
PopObj(:,1) = (1.7057+sum1).*cos(0.5*pi*X(:,1)).*cos(0.5*pi*X(:,2));
PopObj(:,2) = (1.7057+sum1).*cos(0.5*pi*X(:,1)).*sin(0.5*pi*X(:,2));
PopObj(:,3) = (1.7057+sum1).*sin(0.5*pi*X(:,1));
gx          =  PopObj(:,1).^2+PopObj(:,2).^2+PopObj(:,3).^2;
PopCon(:,1) = (gx-9).*(4-gx);
PopCon(:,2) = (gx-3.61).*(3.24-gx);
return
