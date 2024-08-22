function[f,g] = discbrake(x)
if nargin == 0
	prob.nx = 4;
	prob.nf = 2;
	prob.ng = 4;
	prob.bounds(1,:) = [55,80];
	prob.bounds(2,:) = [75,110];
	prob.bounds(3,:) = [1000,3000];
    prob.bounds(4,:) = [11,20];
	f = prob;
else
	[f,g] = discbrake_true(x);
end
return


function [obj,g] = discbrake_true(x)
x1 = x(:,1);
x2 = x(:,2);
x3 = x(:,3);
x4 = x(:,4);

% First original objective function
obj(:,1) = (4.9e-5) .* (x2.^2 - x1.^2) .* (x4 - 1.0);

% Second original objective function
obj(:,2) = ((9.82e6) .* (x2.^2 - x1.^2)) ./ ((x3 .* x4) .* (x2.^3 - x1.^3));

% Reformulated objective functions
g(:,1) = -((x2 - x1) - 20.0);
g(:,2) = -(0.4 - (x3 ./ (3.14 .* (x2.^2 - x1.^2))));
g(:,3) = -(1.0 - (2.22e-3 .* x3 .* (x2.^3 - x1.^3)) ./ power((x2.^2 - x1.^2), 2));
g(:,4) = -(((2.66e-2 .* x3 .* x4 .* (x2.^3 - x1.^3)) ./ (x2.^2 - x1.^2)) - 900.0);
return