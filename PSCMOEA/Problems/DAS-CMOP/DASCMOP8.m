function[f,g] = DASCMOP8(x)
if nargin == 0
	prob.nx = 10;
	prob.nf = 3;
	prob.ng = 7;
    prob.bounds(1:prob.nx,1) = zeros(prob.nx,1);
    prob.bounds(1:prob.nx,2) = ones(prob.nx,1);
	f = prob;
else
	[f,g] = DASCMOP8_true(x);
end
return


function [PopObj,PopCon] = DASCMOP8_true(X)
%% Calculate objective values
x_all = X(:,3:1:end);D = size(X,2);
g_1   = sum((x_all - 0.5) .* (x_all - 0.5) - cos(20.0 * pi .* ( x_all - 0.5)),2);
sum1  = D - 2 + g_1;
PopObj(:,1) = cos(0.5 * pi * X(:,1)) .* cos(0.5 * pi * X(:,2)) + sum1;
PopObj(:,2) = cos(0.5 * pi * X(:,1)) .* sin(0.5 * pi * X(:,2)) + sum1;
PopObj(:,3) = sin(0.5 * pi * X(:,1)) + sum1;

%% Calculate constraint violations
PopCon = CalCon(X);
return


%% Calculate constraint violations
function PopCon = CalCon(X)
x_all  = X(:,3:1:end);D = size(X,2);
g_1    = sum((x_all - 0.5) .* (x_all - 0.5) - cos(20.0 * pi .* ( x_all - 0.5)),2);
sum1   = D - 2 + g_1;
PopCon = Constraint3(X(:,1),X(:,2),sum1);
return


function PopCon = Constraint3(X1,X2,sum1)
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
PopObj(:,1) = X1 .* X2 + sum1;
PopObj(:,2) = X2 .* (1 - X1) +  sum1;
PopObj(:,3) = 1 - X2 + sum1;
% Type-I constraints
PopCon(:,1) = b - sin(a * pi * X1);
PopCon(:,2) = b - cos(a * pi * X2);
% Type-II constraints
PopCon(:,3) = -(e - sum1) .* (sum1 - d);
if DifficultyFactors(2) == 1.0
    PopCon(:,3) = 1e-4 - abs(sum1 - e);
end
% Type-III constraints
x_k = [1.0, 0.0, 0.0, 1.0 / sqrt(3.0)];
y_k = [0.0, 1.0, 0.0, 1.0 / sqrt(3.0)];
z_k = [0.0, 0.0, 1.0, 1.0 / sqrt(3.0)];
for k=1:length(x_k)
    PopCon(:,3+k) = r * r - ((PopObj(:,1) - x_k(k))).^2  -...
        ((PopObj(:,2) - y_k(k))).^2 - ((PopObj(:,3) - z_k(k))).^2;
end
return