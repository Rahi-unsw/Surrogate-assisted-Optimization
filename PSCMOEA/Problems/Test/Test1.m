function[f,g] = Test1(x)
if nargin == 0
	prob.nx = 10;
	prob.nf = 2;
	prob.ng = 1;
    prob.bounds(1:prob.nx,1) = zeros(prob.nx,1);
    prob.bounds(1:prob.nx,2) = ones(prob.nx,1);
	f = prob;
else
	[f,g] = Test1_true(x);
end
return


function [PopObj,PopCon] = Test1_true(PopDec)
%% Calculate objective values
g = 1 + 9*mean(PopDec(:,2:end),2);
PopObj(:,1) = PopDec(:,1).*g;
PopObj(:,2) = (1-PopDec(:,1)).*g;

%% Calculate constraint violations
PopCon = CalCon(PopDec);
return


%% Calculate constraint violations
function PopCon = CalCon(PopDec)
g = 1 + 9*mean(PopDec(:,2:end),2);
Dis    = abs(9-g);
%%%%% Type-I constraints
y1     = Dis.^2-0.25;
y2     = g.*(1.2+0.5*sin(Dis*pi));
PopCon = min([y1,y2],[],2);
return