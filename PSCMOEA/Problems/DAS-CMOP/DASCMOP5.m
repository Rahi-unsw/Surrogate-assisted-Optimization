function[f,g] = DASCMOP5(x)
if nargin == 0
	prob.nx = 10;
	prob.nf = 2;
	prob.ng = 11;
    prob.bounds(1:prob.nx,1) = zeros(prob.nx,1);
    prob.bounds(1:prob.nx,2) = ones(prob.nx,1);
	f = prob;
else
	[f,g] = DASCMOP5_true(x);
end
return


function [PopObj,PopCon] = DASCMOP5_true(X)
%% Calculate objective values
x_all = X(:,2:1:end); D = size(X,2);
g_1   = sum((x_all - 0.5) .* (x_all - 0.5) - cos(20.0 * pi .* ( x_all - 0.5)),2);
sum1  = D - 1 + g_1;
PopObj(:,1) = X(:,1) + sum1;
PopObj(:,2) = 1.0 - sqrt(X(:,1)) + sum1;

%% Calculate constraint violations
PopCon = CalCon(X);
return


%% Calculate constraint violations
function PopCon = CalCon(X)
x_all  = X(:,2:1:end);D = size(X,2);
g_1    = sum((x_all - 0.5) .* (x_all - 0.5) - cos(20.0 * pi .* ( x_all - 0.5)),2);
sum1   = D - 1 + g_1;
PopCon = Constraint2(X(:,1),sum1);
return


function PopCon = Constraint2(X1,sum1)
% set the parameters of constraints
DifficultyFactors = [0.5,0.5,0.5];
% Type-I parameters
a = 20;
b = 2 * DifficultyFactors(1) - 1;
% Type-II parameters
d = 0.5;
if DifficultyFactors(2) == 0.0
    d = 0.0;
end
e = d - log(DifficultyFactors(2));
if isfinite(e) == 0
    e = 1e+30;
end
% Type-III parameters
r = 0.5 * DifficultyFactors(3);
% Calculate objective values
PopObj(:,1) = X1 + sum1;
PopObj(:,2) = 1 - X1 .^ 2 + sum1;
% Type-I constraints
PopCon(:,1) = b -  sin(a * pi * X1);
% Type-II constraints
PopCon(:,2) = -(e - sum1) .* (sum1 - d);
if DifficultyFactors(2) == 1.0
    PopCon(:,2) = 1e-4 - abs(sum1 - e);
end
% Type-III constraints
p_k = [0.0, 1.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 3.0];
q_k = [1.5, 0.5, 2.5, 1.5, 0.5, 3.5, 2.5, 1.5, 0.5];
a_k = 0.3;
b_k = 1.2;
theta_k = -0.25 * pi;
for k=1:length(p_k)
    PopCon(:,2+k) = r - ((PopObj(:,1) - p_k(k)) * cos(theta_k) - (PopObj(:,2) - q_k(k)) * sin(theta_k)).^2 ./ a_k -...
        ((PopObj(:,1) - p_k(k)) * sin(theta_k) + (PopObj(:,2) - q_k(k)) * cos(theta_k)).^2 ./ b_k;
end
return