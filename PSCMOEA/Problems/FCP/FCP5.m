function[f,g] = FCP5(x)
if nargin == 0
	prob.nx = 10;
	prob.nf = 2;
	prob.ng = 1;
    prob.bounds(1:prob.nx,1) = zeros(prob.nx,1);
    prob.bounds(1:prob.nx,2) = ones(prob.nx,1);
	f = prob;
else
	[f,g] = FCP5_true(x);
end
return


function [PopObj,PopCon] = FCP5_true(PopDec)
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
%%%%% Type-III constraints
c1     = log(sqrt((10*PopDec(:,1)-9).^2+(g-3).^2)+0.5);
c2     = log(sqrt((10*PopDec(:,1)-6).^2+(g-6).^2)+0.05);
c3     = (10*PopDec(:,1)-sqrt(2)).^2+(g-10).^2-2;
c4     = 1.2+sin(pi*sqrt(c3+2));
PopCon = min([c1,c2,c3,c4],[],2);
return