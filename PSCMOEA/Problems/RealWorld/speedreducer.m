function[f,g] = speedreducer(x)
if nargin == 0
	prob.nx = 7;
	prob.nf = 2;
	prob.ng = 11;
	prob.bounds(1,:) = [2.6,3.6];
	prob.bounds(2,:) = [0.7,0.8];
	prob.bounds(3,:) = [17,28];
    prob.bounds(4,:) = [7.3,8.3];
    prob.bounds(5,:) = [7.3,8.3];
    prob.bounds(6,:) = [2.9,3.9];
    prob.bounds(7,:) = [5.0,5.5];
	f = prob;
else
	[f,g] = speedreducer_true(x);
end
return


function [obj,g] = speedreducer_true(x)
x1 = x(:,1);
x2 = x(:,2);
x3 = round(x(:,3));
x4 = x(:,4);
x5 = x(:,5);
x6 = x(:,6);
x7 = x(:,7);

% First original objective function (weight)
obj(:,1) = 0.7854 .* x1 .* x2.^2 .* (((10.0 .* x3.^2) ./ 3.0) + (14.933 .* x3) - 43.0934) - (1.508 .* x1 .* (x6.^2 + x7.^2)) + (7.477 .* (x6.^3 + x7.^3)) + (0.7854 .* ((x4 .* x6.^2) + (x5 .* x7.^2)));

% Second original objective function (stress)
tmpVar = power((745.0 .* x4) ./ (x2 .* x3), 2.0)  + 1.69e7;
obj(:,2) =  sqrt(tmpVar) ./ (0.1 .* x6.^3);

% Constraint functions
g(:,1) = -(-(1.0 ./ (x1 .* x2.^2 .* x3)) + (1.0 / 27.0));
g(:,2) = -(-(1.0 ./ (x1 .* x2.^2 .* x3.^2)) + (1.0 / 397.5));
g(:,3) = -(-((x4.^3) ./ (x2 .* x3 .* x6.^4)) + (1.0 / 1.93));
g(:,4) = -(-((x5.^3) ./ (x2 .* x3 .* x7.^4)) + (1.0 / 1.93));
g(:,5) = -(-(x2 .* x3) + 40.0);
g(:,6) = -(-(x1 ./ x2) + 12.0);
g(:,7) = -(-5.0 + (x1 ./ x2));
g(:,8) = -(-1.9 + x4 - (1.5 .* x6));
g(:,9) = -(-1.9 + x5 - (1.1 .* x7));
g(:,10) =  -(-obj(:,2) + 1300.0);
tmpVar = power(((745.0 * x5) ./ (x2 .* x3)), 2.0) + 1.575e8;
g(:,11) = -(-sqrt(tmpVar) ./ (0.1 .* x7.^3) + 1100.0);
return