%% Plot the feasible regions of the MW series
function R = FR(prob_name,M)
switch prob_name
    case 'MW1'
        [x,y] = meshgrid(linspace(0,1,400),linspace(0,1.5,400));
        z     = nan(size(x));
        fes   = x + y - 1 - 0.5*sin(2*pi*(sqrt(2)*y-sqrt(2)*x)).^8 <= 0;
        z(fes & 0.85*x+y>=1) = 0;
        R = {x,y,z};

    case 'MW2'
        [x,y] = meshgrid(linspace(0,1,400),linspace(0,1.5,400));
        z     = nan(size(x));
        fes   = x + y - 1 - 0.5*sin(3*pi*(sqrt(2)*y-sqrt(2)*x)).^8 <= 0;
        z(fes & x+y>=1) = 0;
        R = {x,y,z};

    case 'MW3'
        [x,y] = meshgrid(linspace(0,1,400),linspace(0,1.5,400));
        z     = nan(size(x));
        fes1  = x + y - 1.05 - 0.45*sin(0.75*pi*(sqrt(2)*y-sqrt(2)*x)).^6 <= 0;
        fes2  = 0.85 - x - y + 0.3*sin(0.75*pi*(sqrt(2)*y-sqrt(2)*x)).^2 <= 0;
        z(fes1 & fes2 & x+y>=1) = 0;
        R = {x,y,z};

    case 'MW4'
        if M == 2
            [x,y] = meshgrid(linspace(0,1.2,400));
            z     = nan(size(x));
            fes   = x + y - 1 - 0.4*sin(2.5*pi*(y-x)).^8 <= 0;
            z(fes & x+y>=1) = 0;
            R = {x,y,z};
        elseif M == 3
            a = linspace(0,1,10)';
            x = a*a';
            y = a*(1-a');
            z = (1-a)*ones(size(a'));
            R = {x,y,z};
        else
            R = [];
        end

    case 'MW5'
        [x,y] = meshgrid(linspace(0,1.7,400));
        z     = nan(size(x));
        l1    = atan(y./x);
        l2    = 0.5*pi - 2*abs(l1-0.25*pi);
        fes1  = x.^2 + y.^2 - (1.7-0.2*sin(2*l1)).^2 <= 0;
        fes2  = (1+0.5*sin(6*l2.^3)).^2 - x.^2 - y.^2 <= 0;
        fes3  = (1-0.45*sin(6*l2.^3)).^2 - x.^2 - y.^2 <= 0;
        z(fes1 & fes2 & fes3 & x.^2+y.^2>=1) = 0;
        R = {x,y,z};

    case 'MW6'
        [x,y] = meshgrid(linspace(0,1.7,400));
        z     = nan(size(x));
        l     = cos(6*atan(y./x).^4).^10;
        fes   = (x./(1+0.15*l)).^2 + (y./(1+0.75*l)).^2 - 1 <= 0;
        z(fes & x.^2+y.^2>=1.21) = 0;
        R = {x,y,z};

    case 'MW7'
        [x,y] = meshgrid(linspace(0,1.5,400));
        z     = nan(size(x));
        l     = atan(y./x);
        fes1  = x.^2 + y.^2 - (1.2+0.4*sin(4*l).^16).^2 <= 0;
        fes2  = (1.15-0.2*sin(4*l).^8).^2 - x.^2 - y.^2 <= 0;
        z(fes1 & fes2 & x.^2+y.^2>=1) = 0;
        R = {x,y,z};

    case 'MW8'
        if M == 2
            [x,y] = meshgrid(linspace(0,1.2,400));
            z     = nan(size(x));
            fes   = x.^2 + y.^2 - (1.25-0.5*sin(6*asin(y./sqrt(x.^2+y.^2))).^2).^2 <= 0;
            z(fes & x.^2+y.^2>=1) = 0;
            R = {x,y,z};
        elseif M == 3
            a = linspace(0,pi/2,40)';
            x = sin(a)*cos(a');
            y = sin(a)*sin(a');
            z = cos(a)*ones(size(a'));
            fes     = 1 - (1.25-0.5*sin(6*asin(z)).^2).^2 <= 0;
            z(~fes) = nan;
            R = {x,y,z};
        else
            R = [];
        end

    case 'MW9'
        [x,y] = meshgrid(linspace(0,1.7,400));
        z     = nan(size(x));
        T1    = (1-0.64*x.^2-y).*(1-0.36*x.^2-y);
        T2    = 1.35.^2 - (x+0.35).^2 - y;
        T3    = 1.15.^2 - (x+0.15).^2 - y;
        fes   = min(T1,T2.*T3) <= 0;
        z(fes & x.^0.6+y>=1) = 0;
        R = {x,y,z};

    case 'MW10'
        [x,y] = meshgrid(linspace(0,1,400),linspace(0,1.5,400));
        z     = nan(size(x));
        fes1  = -(2-4*x.^2-y).*(2-8*x.^2-y) <= 0;
        fes2  = (2-2*x.^2-y).*(2-16*x.^2-y) <= 0;
        fes3  = (1-x.^2-y).*(1.2-1.2*x.^2-y) <= 0;
        z(fes1 & fes2 & fes3 & x.^2+y>=1) = 0;
        R = {x,y,z};

    case 'MW11'
        [x,y] = meshgrid(linspace(0,2.1,400));
        z     = nan(size(x));
        fes1  = -(3-x.^2-y).*(3-2*x.^2-y) <= 0;
        fes2  = (3-0.625*x.^2-y).*(3-7*x.^2-y) <= 0;
        fes3  = -(1.62-0.18*x.^2-y).*(1.125-0.125*x.^2-y) <= 0;
        fes4  = (2.07-0.23*x.^2-y).*(0.63-0.07*x.^2-y) <= 0;
        z(fes1 & fes2 & fes3 & fes4 & x.^2+y.^2>=2) = 0;
        R = {x,y,z};

    case 'MW12'
        [x,y] = meshgrid(linspace(0,2,400));
        z     = nan(size(x));
        fes1  = (1-0.8*x-y+0.08*sin(2*pi*(y-x/1.5))).*(1.8-1.125*x-y+0.08*sin(2*pi*(y/1.8-x/1.6))) <= 0;
        fes2  = -(1-0.625*x-y+0.08*sin(2*pi*(y-x/1.6))).*(1.4-0.875*x-y+0.08*sin(2*pi*(y/1.4-x/1.6))) <= 0;
        z(fes1 & fes2 & 0.8*x+0.08*abs(sin(3.2*pi*x))+y>=0.85) = 0;
        R = {x,y,z};

    case 'MW13'
        [x,y] = meshgrid(linspace(0,2,400),linspace(0,4.5,400));
        z     = nan(size(x));
        fes1  = (5-exp(x)-0.5*sin(3*pi*x)-y).*(5-(1+0.4*x)-0.5*sin(3*pi*x)-y) <= 0;
        fes2  = -(5-(1+x+0.5*x.^2)-0.5*sin(3*pi*x)-y).*(5-(1+0.7*x)-0.5*sin(3*pi*x)-y) <= 0;
        z(fes1 & fes2 & exp(x)+abs(0.5*sin(3*pi*x))+y>=5) = 0;
        R = {x,y,z};

    case 'MW14'
        if M == 2
            x      = linspace(0,1.5,100)';
            y      = 6 - exp(x) - 1.5*sin(1.1*pi*x.^2);
            nd     = NDSort([x,y],1)==1;
            x(~nd) = nan;
            R      = [x,y];
        elseif M == 3
            [x,y]  = meshgrid(linspace(0,1.5,20));
            z      = 1/2*(12-exp(x)-1.5*sin(1.1*pi*x.^2)-exp(y)-1.5*sin(1.1*pi*y.^2));
            nd     = reshape(NDSort([x(:),y(:),z(:)],1)==1,size(z));
            z(~nd) = nan;
            R      = {x,y,z};
        else
            R = [];
        end

    case 'DASCMOP1'
        [x,y] = meshgrid(linspace(0.5,2.5,400));
        X1    = (sqrt(1-4*(-x+y-1))-1)/2;
        sum1  = x - X1;
        fes   = all(Constraint1(X1(:),sum1(:))<=0,2);
        z     = nan(size(x));
        z(reshape(fes,size(z))) = 0;
        R = {x,y,z};

    case 'DASCMOP2'
        [x,y] = meshgrid(linspace(0.5,2.5,400));
        X1    = ((sqrt(1-4*(-x+y-1))-1)/2).^2;
        sum1  = x - X1;
        fes   = all(Constraint1(X1(:),sum1(:))<=0,2);
        z     = nan(size(x));
        z(reshape(fes,size(z))) = 0;
        R = {x,y,z};

    case 'DASCMOP3'
        [x,y] = meshgrid(linspace(0.5,2.5,400));
        X1    = ones(size(x));
        E     = inf(size(x));
        for a = 0 : 0.0001 : 1
            E1 = abs(a-1+sqrt(a)-0.5*abs(sin(5*pi*a))-x+y);
            X1(E1<E) = a;
            E = min(E,E1);
        end
        sum1 = x - X1;
        fes  = all(Constraint1(X1(:),sum1(:))<=0,2);
        z    = nan(size(x));
        z(reshape(fes,size(z))) = 0;
        R = {x,y,z};

    case 'DASCMOP4'
        [x,y] = meshgrid(linspace(0.5,2.5,400));
        X1    = (sqrt(1-4*(-x+y-1))-1)/2;
        sum1  = x - X1;
        fes   = all(Constraint2(X1(:),sum1(:))<=0,2);
        z     = nan(size(x));
        z(reshape(fes,size(z))) = 0;
        R = {x,y,z};

    case 'DASCMOP5'
        [x,y] = meshgrid(linspace(0.5,2.5,400));
        X1    = ((sqrt(1-4*(-x+y-1))-1)/2).^2;
        sum1  = x - X1;
        fes   = all(Constraint2(X1(:),sum1(:))<=0,2);
        z     = nan(size(x));
        z(reshape(fes,size(z))) = 0;
        R = {x,y,z};

    case 'DASCMOP6'
        [x,y] = meshgrid(linspace(0.5,2.5,400));
        X1    = ones(size(x));
        E     = inf(size(x));
        for a = 0 : 0.0001 : 1
            E1 = abs(a-1+sqrt(a)-0.5*abs(sin(5*pi*a))-x+y);
            X1(E1<E) = a;
            E = min(E,E1);
        end
        sum1 = x - X1;
        fes  = all(Constraint2(X1(:),sum1(:))<=0,2);
        z    = nan(size(x));
        z(reshape(fes,size(z))) = 0;
        R = {x,y,z};

    case 'DASCMOP7'
        [X1,X2] = meshgrid(linspace(0,1,100));
        x = X1.*X2;
        y = (1-X1).*X2;
        z = (1-X2);
        fes = all(Constraint3(X1(:),X2(:),0.5)<=0,2);
        z(reshape(~fes,size(z))) = nan;
        R = {x+0.5,y+0.5,z+0.5};

    case 'DASCMOP8'
        [X1,X2] = meshgrid(linspace(0,1,100));
        x = cos(0.5*pi*X1).*cos(0.5*pi*X2);
        y = cos(0.5*pi*X1).*sin(0.5*pi*X2);
        z = sin(0.5*pi*X1);
        fes = all(Constraint3(X1(:),X2(:),0.5)<=0,2);
        z(reshape(~fes,size(z))) = nan;
        R = {x+0.5,y+0.5,z+0.5};

    case 'DASCMOP9'
        [X1,X2] = meshgrid(linspace(0,1,100));
        x = cos(0.5*pi*X1).*cos(0.5*pi*X2);
        y = cos(0.5*pi*X1).*sin(0.5*pi*X2);
        z = sin(0.5*pi*X1);
        fes = all(Constraint3(X1(:),X2(:),0.5)<=0,2);
        z(reshape(~fes,size(z))) = nan;
        R = {x+0.5,y+0.5,z+0.5};

    case 'LIRCMOP1'
        R = [];

    case 'LIRCMOP2'
        R = [];

    case 'LIRCMOP3'
        R = [];

    case 'LIRCMOP4'
        R = [];

    case 'LIRCMOP5'
        [x,y] = meshgrid(linspace(0.7057,5,400));
        fes   = all(Constraint4([x(:),y(:)])<=0,2);
        z     = nan(size(x));
        z(reshape(fes,size(z)) & sqrt(x-0.7057)+y>=1.7057) = 0;
        R = {x,y,z};

    case 'LIRCMOP6'
        [x,y] = meshgrid(linspace(0.7057,5,400));
        fes   = all(Constraint5([x(:),y(:)])<=0,2);
        z     = nan(size(x));
        z(reshape(fes,size(z)) & (x-0.7057).^2+y>=1.7057) = 0;
        R = {x,y,z};

    case 'LIRCMOP7'
        [x,y] = meshgrid(linspace(0.7057,5,400));
        fes   = all(Constraint6([x(:),y(:)])<=0,2);
        z     = nan(size(x));
        z(reshape(fes,size(z)) & sqrt(x-0.7057)+y>=1.7057) = 0;
        R = {x,y,z};

    case 'LIRCMOP8'
        [x,y] = meshgrid(linspace(0.7057,5,400));
        fes   = all(Constraint7([x(:),y(:)])<=0,2);
        z     = nan(size(x));
        z(reshape(fes,size(z)) & sqrt(x-0.7057)+y>=1.7057) = 0;
        R = {x,y,z};

    case 'LIRCMOP9'
        [x,y] = meshgrid(linspace(0,5,400));
        fes   = all(Constraint8([x(:),y(:)])<=0,2);
        z     = nan(size(x));
        z(reshape(fes,size(z)) & (x/1.7057).^2+y/1.7057>=1) = 0;
        R = {x,y,z};

    case 'LIRCMOP10'
        [x,y] = meshgrid(linspace(0,5,400));
        fes   = all(Constraint9([x(:),y(:)])<=0,2);
        z     = nan(size(x));
        z(reshape(fes,size(z)) & sqrt(x/1.7057)+y/1.7057>=1) = 0;
        R = {x,y,z};

    case 'LIRCMOP11'
        [x,y] = meshgrid(linspace(0,5,400));
        fes   = all(Constraint10([x(:),y(:)])<=0,2);
        z     = nan(size(x));
        z(reshape(fes,size(z)) & sqrt(x/1.7057)+y/1.7057>=1) = 0;
        R = {x,y,z};

    case 'LIRCMOP12'
        [x,y] = meshgrid(linspace(0,5,400));
        fes   = all(Constraint11([x(:),y(:)])<=0,2);
        z     = nan(size(x));
        z(reshape(fes,size(z)) & (x/1.7057).^2+y/1.7057>=1) = 0;
        R = {x,y,z};

    case 'LIRCMOP13'
        a = linspace(0,pi/2,10)';
        R = {sin(a)*cos(a')*1.7057,sin(a)*sin(a')*1.7057,cos(a)*ones(size(a'))*1.7057};

    case 'LIRCMOP14'
        a = linspace(0,pi/2,10)';
        R = {sin(a)*cos(a')*sqrt(3.0625),sin(a)*sin(a')*sqrt(3.0625),cos(a)*ones(size(a'))*sqrt(3.0625)};

    case 'FCP1'
        R = [];

    case 'Test1'
        [x1,y1] = meshgrid(linspace(0,1,400),linspace(0,1,400));
        g = 1 + 9*y1;
        x = x1.*g;
        y = (1-x1).*g;
        %% Calculate constraint violations
        Dis    = abs(9-g);
        %%%%% Type-I constraints
        y1     = Dis.^2-0.25;
        y2     = g.*(1.2+0.5*sin(Dis*pi));
        PopCon = zeros(size(x1));
        for i = 1:size(x1,1)
            for j = 1:size(x1,1)
                PopCon(i,j) = min(y1(i,j),y2(i,j));
            end
        end
        z     = nan(size(x1));
        fes   = PopCon <= 0;
        z(fes) = 0;
        R = {x,y,z};

    case 'Test2'
        [x1,y1] = meshgrid(linspace(0,1,400),linspace(0,1,400));
        g = 1 + 9*y1;
        t = mod(floor(100*g),2);
        g = g + t.*(g-9).^2;
        x = cos(0.5*pi*x1).*g;
        y = sin(0.5*pi*x1).*g;

        %% Calculate constraint violations
        Dis    = abs(9-g);
        %%%%% Type-II constraints
        y1     = Dis.^2-0.25;
        y2     = 1./(g+1e-6).*(1.2+0.5.*sin(Dis*pi));
        PopCon = zeros(size(x1));
        for i = 1:size(x1,1)
            for j = 1:size(x1,1)
                PopCon(i,j) = min(y1(i,j),y2(i,j));
            end
        end

        z     = nan(size(x1));
        fes   = PopCon <= 0;
        z(fes) = 0;
        R = {x,y,z};

    case 'FCP2'
        R = [];

    case 'FCP3'
        R = [];

    case 'FCP4'
        R = [];

    case 'FCP5'
        R = [];

    otherwise
        R = [];
end
end



function [FrontNo,MaxFNo] = NDSort(varargin)
%NDSort - Do non-dominated sorting by efficient non-dominated sort.
%
%   FrontNo = NDSort(F,s) does non-dominated sorting on F, where F is the
%   matrix of objective values of a set of individuals, and s is the number
%   of individuals to be sorted at least. FrontNo(i) denotes the front
%   number of the i-th individual. The individuals have not been sorted are
%   assigned a front number of inf.
%
%   FrontNo = NDSort(F,C,s) does non-dominated sorting based on constrained
%   domination, where C is the matrix of constraint values of the
%   individuals. In this case, feasible solutions always dominate
%   infeasible solutions, and one infeasible solution dominates another
%   infeasible solution if the former has a smaller overall constraint
%   violation than the latter.
%
%   In particular, s = 1 indicates finding only the first non-dominated
%   front, s = size(F,1)/2 indicates sorting only half the population
%   (which is often used in the algorithm), and s = inf indicates sorting
%   the whole population.
%
%   [FrontNo,K] = NDSort(...) also returns the maximum front number besides
%   inf.
%
%   Example:
%       [FrontNo,MaxFNo] = NDSort(PopObj,1)
%       [FrontNo,MaxFNo] = NDSort(PopObj,PopCon,inf)

%------------------------------- Reference --------------------------------
% [1] X. Zhang, Y. Tian, R. Cheng, and Y. Jin, An efficient approach to
% nondominated sorting for evolutionary multiobjective optimization, IEEE
% Transactions on Evolutionary Computation, 2015, 19(2): 201-213.
% [2] X. Zhang, Y. Tian, R. Cheng, and Y. Jin, A decision variable
% clustering based evolutionary algorithm for large-scale many-objective
% optimization, IEEE Transactions on Evolutionary Computation, 2018, 22(1):
% 97-112.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

PopObj = varargin{1};
[N,M]  = size(PopObj);
if nargin == 2
    nSort  = varargin{2};
else
    PopCon = varargin{2};
    nSort  = varargin{3};
    Infeasible           = any(PopCon>0,2);
    PopObj(Infeasible,:) = repmat(max(PopObj,[],1),sum(Infeasible),1) + repmat(sum(max(0,PopCon(Infeasible,:)),2),1,M);
end
if M < 3 || N < 500
    % Use efficient non-dominated sort with sequential search (ENS-SS)
    [FrontNo,MaxFNo] = ENS_SS(PopObj,nSort);
else
    % Use tree-based efficient non-dominated sort (T-ENS)
    [FrontNo,MaxFNo] = T_ENS(PopObj,nSort);
end
end

function [FrontNo,MaxFNo] = ENS_SS(PopObj,nSort)
[PopObj,~,Loc] = unique(PopObj,'rows');
Table   = hist(Loc,1:max(Loc));
[N,M]   = size(PopObj);
FrontNo = inf(1,N);
MaxFNo  = 0;
while sum(Table(FrontNo<inf)) < min(nSort,length(Loc))
    MaxFNo = MaxFNo + 1;
    for i = 1 : N
        if FrontNo(i) == inf
            Dominated = false;
            for j = i-1 : -1 : 1
                if FrontNo(j) == MaxFNo
                    m = 2;
                    while m <= M && PopObj(i,m) >= PopObj(j,m)
                        m = m + 1;
                    end
                    Dominated = m > M;
                    if Dominated || M == 2
                        break;
                    end
                end
            end
            if ~Dominated
                FrontNo(i) = MaxFNo;
            end
        end
    end
end
FrontNo = FrontNo(:,Loc);
end

function [FrontNo,MaxFNo] = T_ENS(PopObj,nSort)
[PopObj,~,Loc] = unique(PopObj,'rows');
Table     = hist(Loc,1:max(Loc));
[N,M]     = size(PopObj);
FrontNo   = inf(1,N);
MaxFNo    = 0;
Forest    = zeros(1,N);
Children  = zeros(N,M-1);
LeftChild = zeros(1,N) + M;
Father    = zeros(1,N);
Brother   = zeros(1,N) + M;
[~,ORank] = sort(PopObj(:,2:M),2,'descend');
ORank     = ORank + 1;
while sum(Table(FrontNo<inf)) < min(nSort,length(Loc))
    MaxFNo = MaxFNo + 1;
    root   = find(FrontNo==inf,1);
    Forest(MaxFNo) = root;
    FrontNo(root)  = MaxFNo;
    for p = 1 : N
        if FrontNo(p) == inf
            Pruning = zeros(1,N);
            q = Forest(MaxFNo);
            while true
                m = 1;
                while m < M && PopObj(p,ORank(q,m)) >= PopObj(q,ORank(q,m))
                    m = m + 1;
                end
                if m == M
                    break;
                else
                    Pruning(q) = m;
                    if LeftChild(q) <= Pruning(q)
                        q = Children(q,LeftChild(q));
                    else
                        while Father(q) && Brother(q) > Pruning(Father(q))
                            q = Father(q);
                        end
                        if Father(q)
                            q = Children(Father(q),Brother(q));
                        else
                            break;
                        end
                    end
                end
            end
            if m < M
                FrontNo(p) = MaxFNo;
                q = Forest(MaxFNo);
                while Children(q,Pruning(q))
                    q = Children(q,Pruning(q));
                end
                Children(q,Pruning(q)) = p;
                Father(p) = q;
                if LeftChild(q) > Pruning(q)
                    Brother(p)   = LeftChild(q);
                    LeftChild(q) = Pruning(q);
                else
                    bro = Children(q,LeftChild(q));
                    while Brother(bro) < Pruning(q)
                        bro = Children(q,Brother(bro));
                    end
                    Brother(p)   = Brother(bro);
                    Brother(bro) = Pruning(q);
                end
            end
        end
    end
end
FrontNo = FrontNo(:,Loc);
end


function PopCon = Constraint1(X1,sum1)
% set the parameters of constraints
DifficultyFactors = [0,0.5,0.5];
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
PopCon(:,1) = b - sin(a * pi * X1);
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
end


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
end



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
end


function PopCon = Constraint4(PopObj)
p     = [1.6,2.5];
q     = [1.6,2.5];
a     = [2,2];
b     = [4,8];
r     = 0.1;
theta = -0.25 * pi;
for k = 1 : 2
    PopCon(:,k) = r - ((PopObj(:,1)-p(k))*cos(theta)-(PopObj(:,2)-q(k))*sin(theta)).^2/(a(k)^2) - ...
        ((PopObj(:,1)-p(k))*sin(theta)+(PopObj(:,2)-q(k))*cos(theta)).^2/(b(k)^2);
end
end


function PopCon = Constraint5(PopObj)
p     = [1.8,2.8];
q     = [1.8,2.8];
a     = [2,2];
b     = [8,8];
r     = 0.1;
theta = -0.25 * pi;
for k = 1 : 2
    PopCon(:,k) = r - ((PopObj(:,1)-p(k))*cos(theta)-(PopObj(:,2)-q(k))*sin(theta)).^2/(a(k)^2) -...
        ((PopObj(:,1)-p(k))*sin(theta)+(PopObj(:,2)-q(k))*cos(theta)).^2/(b(k)^2);
end
end


function PopCon = Constraint6(PopObj)
p     = [1.2,2.25,3.5];
q     = [1.2,2.25,3.5];
a     = [2,2.5,2.5];
b     = [6,12,10];
r     = 0.1;
theta = -0.25*pi;
for k = 1 : 3
    PopCon(:,k) = r - ((PopObj(:,1)-p(k))*cos(theta)-(PopObj(:,2)-q(k))*sin(theta)).^2/(a(k)^2) -...
        ((PopObj(:,1)-p(k))*sin(theta)+(PopObj(:,2)-q(k))*cos(theta)).^2/(b(k)^2);
end
end

function PopCon = Constraint7(PopObj)
p     = [1.2,2.25,3.5];
q     = [1.2,2.25,3.5];
a     = [2,2.5,2.5];
b     = [6,12,10];
r     = 0.1;
theta = -0.25*pi;
for k = 1 : 3
    PopCon(:,k) = r - ((PopObj(:,1)-p(k))*cos(theta)-(PopObj(:,2)-q(k))*sin(theta)).^2/(a(k)^2) -...
        ((PopObj(:,1)-p(k))*sin(theta)+(PopObj(:,2)-q(k))*cos(theta)).^2/(b(k)^2);
end
end

function PopCon = Constraint8(PopObj)
p     = 1.4;
q     = 1.4;
a     = 1.5;
b     = 6;
r     = 0.1;
theta = -0.25 * pi;
alpha = 0.25 * pi;
PopCon(:,1) = r - (((PopObj(:,1)-p)*cos(theta)-(PopObj(:,2)-q)*sin(theta)).^2)/(a^2)-...
    (((PopObj(:,1)-p)*sin(theta)+(PopObj(:,2)-q)*cos(theta)).^2)/(b^2);
PopCon(:,2) = 2 - PopObj(:,1)*sin(alpha) - PopObj(:,2)*cos(alpha) + sin(4*pi*(PopObj(:,1)*cos(alpha)-PopObj(:,2)*sin(alpha)));
end

function PopCon = Constraint9(PopObj)
p     = 1.1;
q     = 1.2;
a     = 2;
b     = 4;
r     = 0.1;
theta = -0.25 * pi;
alpha = 0.25 * pi;
PopCon(:,1)= r - (((PopObj(:,1)-p)*cos(theta)-(PopObj(:,2)-q)*sin(theta)).^2)/(a^2)-...
    (((PopObj(:,1)-p)*sin(theta)+(PopObj(:,2)-q)*cos(theta)).^2)/(b^2);
PopCon(:,2) = 1 - PopObj(:,1)*sin(alpha) - PopObj(:,2)*cos(alpha) + sin(4*pi*(PopObj(:,1)*cos(alpha)-PopObj(:,2)*sin(alpha)));
end


function PopCon = Constraint10(PopObj)
p     = 1.2;
q     = 1.2;
a     = 1.5;
b     = 5;
r     = 0.1;
theta = -0.25 * pi;
alpha = 0.25 * pi;
PopCon(:,1) = r - (((PopObj(:,1)-p)*cos(theta)-(PopObj(:,2)-q)*sin(theta)).^2)/(a^2)-...
    (((PopObj(:,1)-p)*sin(theta)+(PopObj(:,2)-q)*cos(theta)).^2)/(b^2);
PopCon(:,2) = 2.1 - PopObj(:,1)*sin(alpha) - PopObj(:,2)*cos(alpha) + sin(4*pi*(PopObj(:,1)*cos(alpha)-PopObj(:,2)*sin(alpha)));
end


function PopCon = Constraint11(PopObj)
p     = 1.6;
q     = 1.6;
a     = 1.5;
b     = 6;
r     = 0.1;
theta = -0.25 * pi;
alpha = 0.25 * pi;
PopCon(:,1) = r - (((PopObj(:,1)-p)*cos(theta)-(PopObj(:,2)-q)*sin(theta)).^2)/(a^2)-...
    (((PopObj(:,1)-p)*sin(theta)+(PopObj(:,2)-q)*cos(theta)).^2)/(b^2);
PopCon(:,2) = 2.5 - PopObj(:,1)*sin(alpha) - PopObj(:,2)*cos(alpha) + sin(4*pi*(PopObj(:,1)*cos(alpha)-PopObj(:,2)*sin(alpha)));
end