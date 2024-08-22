function[f,g] = LIRCMOP4(x)
if nargin == 0
	prob.nx = 10;
	prob.nf = 2;
	prob.ng = 6;
    prob.bounds(1:prob.nx,1) = zeros(prob.nx,1);
    prob.bounds(1:prob.nx,2) = ones(prob.nx,1);
	f = prob;
else
	[f,g] = LIRCMOP4_true(x);
end
return


function [PopObj,PopCon] = LIRCMOP4_true(X)
prob = feval('LIRCMOP4');
%% Calculate objective values
X = max(min(X,repmat(prob.bounds(:,2)',size(X,1),1)),repmat(prob.bounds(:,1)',size(X,1),1));
x_odd       = X(:,3:2:end);
x_even      = X(:,2:2:end);
len_odd     = size(x_odd,2);
len_even    = size(x_even,2);
g_1         = sum((x_odd - repmat(X(:,1),1,len_odd)).^2,2);
g_2         = sum((x_even - repmat(X(:,1),1,len_even)).^2,2);
PopObj(:,1) = X(:,1) + g_1;
PopObj(:,2) = 1 - sqrt(X(:,1)) + g_2;
eps = 1e-4;
PopCon(:,1) = (0.5 - g_1).*(0.51 - g_1) - eps;
PopCon(:,2) = -(0.5 - g_1).*(0.51 - g_1) - eps;
PopCon(:,3) = (0.5 - g_2).*(0.51 - g_2) - eps;
PopCon(:,4) = -(0.5 - g_2).*(0.51 - g_2) - eps;
PopCon(:,5) = 0.5 - sin(20 * pi * X(:,1)) - eps;
PopCon(:,6) = -0.5 + sin(20 * pi * X(:,1)) - eps;
return


