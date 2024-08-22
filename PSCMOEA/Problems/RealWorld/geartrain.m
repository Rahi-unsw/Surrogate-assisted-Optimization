function[f,g] = geartrain(x)
if nargin == 0
	prob.nx = 4;
	prob.nf = 2;
	prob.ng = 1;
	prob.bounds(1,:) = [12,60];
	prob.bounds(2,:) = [12,60];
    prob.bounds(3,:) = [12,60];
    prob.bounds(4,:) = [12,60];
	f = prob;
else
	[f,g] = geartrain_true(x);
end
return


function [obj,g] = geartrain_true(x)
x1 = round(x(:,1));
x2 = round(x(:,2));
x3 = round(x(:,3));
x4 = round(x(:,4));

% First original objective function
obj(:,1) = abs(6.931 - ((x3 ./ x1) .* (x4 ./ x2)));

% Second original objective function (the maximum value among the four variables)
A = [x1 x2 x3 x4];
obj(:,2) = max(A,[],2);

% Constraint functions
g(:,1) = -(0.5 - (obj(:,1) ./ 6.931));
return