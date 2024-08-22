function[f,g] = LIRCMOP3(x)
if nargin == 0
	prob.nx = 10;
	prob.nf = 2;
	prob.ng = 3;
    prob.bounds(1:prob.nx,1) = zeros(prob.nx,1);
    prob.bounds(1:prob.nx,2) = ones(prob.nx,1);
	f = prob;
else
	[f,g] = LIRCMOP3_true(x);
end
return


function [PopObj,PopCon] = LIRCMOP3_true(X)
prob = feval('LIRCMOP3');
%% Calculate objective values
X = max(min(X,repmat(prob.bounds(:,2)',size(X,1),1)),repmat(prob.bounds(:,1)',size(X,1),1));
x_odd       = X(:,3:2:end);
x_even      = X(:,2:2:end);
len_odd     = size(x_odd,2);
len_even    = size(x_even,2);
g_1         = sum((x_odd - repmat(X(:,1),1,len_odd)).^2,2);
g_2         = sum((x_even - repmat(X(:,1),1,len_even)).^2,2);
PopObj(:,1) = X(:,1) + g_1;
PopObj(:,2) = 1 - X(:,1) .^ 2 + g_2;
PopCon(:,1) = (0.5 - g_1).*(0.51 - g_1);
PopCon(:,2) = (0.5 - g_2).*(0.51 - g_2);
PopCon(:,3) = 0.5 - sin(20 * pi * X(:,1));
return


