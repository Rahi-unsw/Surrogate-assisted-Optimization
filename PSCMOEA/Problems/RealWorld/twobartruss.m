function[f,g] = twobartruss(x)
if nargin == 0
	prob.nx = 3;
	prob.nf = 2;
	prob.ng = 3;
	prob.bounds(1,:) = [0.00001,100.0];
	prob.bounds(2,:) = [0.00001,100.0];
	prob.bounds(3,:) = [1.0,3.0];
	f = prob;
else
	[f,g] = twobartruss_true(x);
end
return


function [obj,g] = twobartruss_true(x)
x1 = x(:,1);
x2 = x(:,2);
x3 = x(:,3);

% First original objective function
obj(:,1) = (x1 .* sqrt(16.0 + (x3 .* x3))) + (x2 .* sqrt(1.0 + (x3 .* x3)));

% Second original objective function
obj(:,2) = (20.0 .* sqrt(16.0 + (x3 .* x3))) ./ ((x1 .* x3));

% Constraint functions
g(:,1) = -(0.1 - obj(:,1));
g(:,2) = -(100000.0 - obj(:,2));
g(:,3) = -(100000 - ((80.0 .* sqrt(1.0 + x3 .* x3)) ./ (x3 .* x2)));
return