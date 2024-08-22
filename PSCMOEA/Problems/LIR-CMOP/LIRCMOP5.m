function[f,g] = LIRCMOP5(x)
if nargin == 0
	prob.nx = 10;
	prob.nf = 2;
	prob.ng = 2;
    prob.bounds(1:prob.nx,1) = zeros(prob.nx,1);
    prob.bounds(1:prob.nx,2) = ones(prob.nx,1);
	f = prob;
else
	[f,g] = LIRCMOP5_true(x);
end
return


function [PopObj,PopCon] = LIRCMOP5_true(X)
prob = feval('LIRCMOP5');
%% Calculate objective values
X = max(min(X,repmat(prob.bounds(:,2)',size(X,1),1)),repmat(prob.bounds(:,1)',size(X,1),1));
[popsize,variable_length] = size(X);
sum1 = zeros(popsize,1);
sum2 = zeros(popsize,1);
for j = 2 : variable_length
    if mod(j,2) == 1
        sum1 = sum1+(X(:,j)-sin((0.5*j/variable_length*pi)*X(:,1))).^2;
    else
        sum2 = sum2+(X(:,j)-cos((0.5*j/variable_length*pi)*X(:,1))).^2;
    end
end
gx          = 0.7057;
PopObj(:,1) = X(:,1)+10*sum1+gx;
PopObj(:,2) = 1-X(:,1).^0.5+10.*sum2+gx;
PopCon = Constraint(PopObj);
return


function PopCon = Constraint(PopObj)
p     = [1.6,2.5];
q     = [1.6,2.5];
a     = [2,2];
b     = [4,8];
r     = 0.1;
theta = -0.25 * pi;
for k = 1 : 2
    PopCon(:,k) = r - ((PopObj(:,1)-p(k))*cos(theta)-(PopObj(:,2)-q(k))*sin(theta)).^2/(a(k)^2) - ...
        ((PopObj(:,1)-p(k))*sin(theta)+(PopObj(:,2)-q(k))*cos(theta)).^2/(b(k)^2);
end
return

