function [f,g] = MaF9_M3_D2(x)
if nargin == 0
    prob.nf = 3;
    prob.ng = 0;
    prob.nx = 2;
    prob.cast = [];
    for i = 1:prob.nx
        prob.bounds(i,:) = [-100,100];
    end
    f = prob;
    g = [];
else
    [f,g] = MaF9_M3_D2_true(x);
end
return

function [f,gg] = MaF9_M3_D2_true(PopDec)
% Generate vertexes
M   = 3;
% Generate vertexes
Points = [];
[thera,rho] = cart2pol(0,1);
[Points(:,1),Points(:,2)] = pol2cart(thera-(1:M)*2*pi/M,rho);
% Generate infeasible polygons
head = repmat((1:M)',ceil(M/2-2),1);
tail = repmat(1:ceil(M/2-2),M,1);
tail = head + tail(:);
Polygons = cell(1,length(head));
for i = 1 : length(Polygons)
    Polygons{i} = Points(mod((head(i):tail(i))-1,M)+1,:);
    Polygons{i} = [Polygons{i};repmat(2*Intersection(Points(mod([head(i)-1,head(i),tail(i),tail(i)+1]-1,M)+1,:)),size(Polygons{i},1),1)-Polygons{i}];
end
obj.Global.lower = [-100,-100];obj.Global.upper = [100,100];
PopDec = CalDec(obj,PopDec,Polygons,Points);
f = zeros(size(PopDec,1),size(Points,1));
for m = 1 : size(Points,1)
    f(:,m) = Point2Line(PopDec,Points(mod(m-1:m,size(Points,1))+1,:));
end
gg = [];
return


function PopDec = CalDec(obj,PopDec,Polygons,Points)
Infeasible = getInfeasible(PopDec,Polygons,Points);
while any(Infeasible)
    PopDec(Infeasible,:) = unifrnd(repmat(obj.Global.lower,sum(Infeasible),1),repmat(obj.Global.upper,sum(Infeasible),1));
    Infeasible           = getInfeasible(PopDec,Polygons,Points);
end
return


function r = Intersection(p)
if p(1,1) == p(2,1)
    r(1) = p(1,1);
    r(2) = p(3,2)+(r(1)-p(3,1))*(p(3,2)-p(4,2))/(p(3,1)-p(4,1));
elseif p(3,1) == p(4,1)
    r(1) = p(3,1);
    r(2) = p(1,2)+(r(1)-p(1,1))*(p(1,2)-p(2,2))/(p(1,1)-p(2,1));
else
    k1   = (p(1,2)-p(2,2))/(p(1,1)-p(2,1));
    k2   = (p(3,2)-p(4,2))/(p(3,1)-p(4,1));
    r(1) = (k1*p(1,1)-k2*p(3,1)+p(3,2)-p(1,2))/(k1-k2);
    r(2) = p(1,2)+(r(1)-p(1,1))*k1;
end
return

function Infeasible = getInfeasible(PopDec,Polygons,Points)
Infeasible = false(size(PopDec,1),1);
for i = 1 : length(Polygons)
    Infeasible = Infeasible | inpolygon(PopDec(:,1),PopDec(:,2),Polygons{i}(:,1),Polygons{i}(:,2));
end
Infeasible = Infeasible & ~inpolygon(PopDec(:,1),PopDec(:,2),Points(:,1),Points(:,2));
return

function Distance = Point2Line(PopDec,Line)
Distance = abs((Line(1,1)-PopDec(:,1)).*(Line(2,2)-PopDec(:,2))-(Line(2,1)-PopDec(:,1)).*(Line(1,2)-PopDec(:,2)))./sqrt((Line(1,1)-Line(2,1)).^2+(Line(1,2)-Line(2,2)).^2);
return