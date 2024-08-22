function[f,g] = LIRCMOP10(x)
if nargin == 0
	prob.nx = 10;
	prob.nf = 2;
	prob.ng = 2;
    prob.bounds(1:prob.nx,1) = zeros(prob.nx,1);
    prob.bounds(1:prob.nx,2) = ones(prob.nx,1);
	f = prob;
else
	[f,g] = LIRCMOP10_true(x);
end
return


function [PopObj,PopCon] = LIRCMOP10_true(X)
prob = feval('LIRCMOP10');
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
PopObj(:,1) = 1.7057*X(:,1).*(10*sum1+1);
PopObj(:,2) = 1.7057*(1-X(:,1).^0.5).*(10*sum2+1);
PopCon = Constraint(PopObj);
return


function PopCon = Constraint(PopObj)
p     = 1.1;
q     = 1.2;
a     = 2;
b     = 4;
r     = 0.1;
theta = -0.25 * pi;
alpha = 0.25 * pi;
PopCon(:,1)= r - (((PopObj(:,1)-p)*cos(theta)-(PopObj(:,2)-q)*sin(theta)).^2)/(a^2)-...
    (((PopObj(:,1)-p)*sin(theta)+(PopObj(:,2)-q)*cos(theta)).^2)/(b^2);
PopCon(:,2) = 1 - PopObj(:,1)*sin(alpha) - PopObj(:,2)*cos(alpha) + sin(4*pi*(PopObj(:,1)*cos(alpha)-PopObj(:,2)*sin(alpha)));
return