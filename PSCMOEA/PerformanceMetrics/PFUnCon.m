%% True pareto equations of the given problems
function P = PFUnCon(prob_name,N,M)
switch prob_name
    case 'MW1'
        R(:,1) = (0:1/(N-1):1)';
        R(:,2) = 1 - 0.85*R(:,1);
        P = R;

    case 'MW2'
        R(:,1) = (0:1/(N-1):1)';
        R(:,2) = 1 - R(:,1);
        P = R;

    case 'MW3'
        R(:,1)  = (0:1/(N-1):1)';
        R(:,2)  = 1 - R(:,1);
        P = R;

    case 'MW4'
        R = UniformPoint(N,M);
        P = R;

    case 'MW5'
        R(:,1) = (0:1/(N-1):1)';
        R(:,2) = 1 - R(:,1);
        R = R./repmat(sqrt(sum(R.^2,2)),1,2);
        P = R;

    case 'MW6'
        R(:,1) = (0:1/(N-1):1)';
        R(:,2) = 1 - R(:,1);
        R = R./repmat(sqrt(sum(R.^2,2)/1.21),1,2);
        P = R;

    case 'MW7'
        R(:,1) = (0:1/(N-1):1)';
        R(:,2) = 1 - R(:,1);
        R = R./repmat(sqrt(sum(R.^2,2)),1,2);
        R = R(NDSort(R,1)==1,:);
        P = R;

    case 'MW8'
        R = UniformPoint(N,M);
        R = R./repmat(sqrt(sum(R.^2,2)),1,M);
        P = R;

    case 'MW9'
        R(:,1)  = (0:1/(N-1):1)';
        R(:,2)  = 1 - R(:,1).^0.6;
        R = R(NDSort(R,1)==1,:);
        P = R;

    case 'MW10'
        R(:,1)  = (0:1/(N-1):1)';
        R(:,2)  = 1 - R(:,1).^2;
        R = R(NDSort(R,1)==1,:);
        P = R;

    case 'MW11'
        R(:,1)  = (0:1/(N-1):1)';
        R(:,2)  = 1 - R(:,1);
        R       = R./repmat(sqrt(sum(R.^2,2)/2),1,2);
        R = R(NDSort(R,1)==1,:);
        P = R;

    case 'MW12'
        R(:,1)  = (0:1/(N-1):1)';
        R(:,2)  = 0.85 - 0.8*R(:,1) - 0.08*abs(sin(3.2*pi*R(:,1)));
        P = R;

    case 'MW13'
        R(:,1)  = (0:1.5/(N-1):1.5)';
        R(:,2)  = 5 - exp(R(:,1)) - 0.5*abs(sin(3*pi*R(:,1)));
        R = R(NDSort(R,1)==1,:);
        P = R;

    case 'MW14'
        interval     = [0,0.731000,1.331000,1.500000];
        median       = (interval(2)-interval(1))/(interval(4)-interval(3)+interval(2)-interval(1));
        X            = UniformPoint(N,M-1,'grid');
        %         X(X<=median) = X(X<=median)*(interval(2)-interval(1))/median+interval(1);
        %         X(X>median)  = (X(X>median)-median)*(interval(4)-interval(3))/(1-median)+interval(3);
        R            = [X,1/(M-1)*sum(6 - exp(X) - 1.5*sin(1.1*pi*X.^2),2)];
        P = R;

    case 'CF1'
        R(:,1) = (0:1/20:1)';
        R(:,2) = 1 - R(:,1);
        P = R;

    case 'CF2'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - sqrt(R(:,1));
        P = R;

    case 'CF3'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - R(:,1).^2;
        P = R;

    case 'CF4'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - R(:,1);
        P = R;

    case 'CF5'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - R(:,1);
        P = R;

    case 'CF6'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = (1-R(:,1)).^2;
        P = R;

    case 'CF7'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = (1-R(:,1)).^2;
        P = R;

    case 'CF8'
        N      = ceil(N/5)*5;
        R      = zeros(N,3);
        R(:,3) = repmat(sin((0:1/(N/5-1):1).*pi/2)',5,1);
        for i = 0 : 4
            R(i*N/5+1:(i+1)*N/5,1) = sqrt(i/4*(1-R(i*N/5+1:(i+1)*N/5,3).^2));
        end
        R(:,2) = sqrt(max(1-R(:,1).^2-R(:,3).^2,0));
        P = R;

    case 'CF9'
        R = UniformPoint(N,3);
        R = R./repmat(sqrt(sum(R.^2,2)),1,3);
        P = R;

    case 'CF10'
        R = UniformPoint(N,3);
        R = R./repmat(sqrt(sum(R.^2,2)),1,3);
        P = R;

    case 'DASCMOP1'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - R(:,1).^2;
        P = R;

    case 'DASCMOP2'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - sqrt(R(:,1));
        P = R;

    case 'DASCMOP3'
        R = [0.5000,1.5000;0.5010,1.4762;0.5020,1.4710;0.5030,1.4688;0.5040,1.4681;0.6502,1.4652;0.7002,1.0541;
            0.9044,0.8986;1.1066,0.7729;1.3008,0.6114;1.5000,0.5000;0.9069,0.8951;1.1126,0.7727;0.9129,0.8950;
            1.1151,0.7690;0.9153,0.8914;1.1175,0.7653;1.1200,0.7616;0.9213,0.8913;1.1260,0.7613;1.1285,0.7576];
        P = R;

    case 'DASCMOP4'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - R(:,1).^2;
        P = R;

    case 'DASCMOP5'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - sqrt(R(:,1));
        P = R;

    case 'DASCMOP6'
        R = [0.5000,1.5000;0.5010,1.4762;0.5020,1.4710;0.5030,1.4688;0.5040,1.4681;0.6502,1.4652;0.7002,1.0541;
            0.9044,0.8986;1.1066,0.7729;1.3008,0.6114;1.5000,0.5000;0.9069,0.8951;1.1126,0.7727;0.9129,0.8950;
            1.1151,0.7690;0.9153,0.8914;1.1175,0.7653;1.1200,0.7616;0.9213,0.8913;1.1260,0.7613;1.1285,0.7576];
        P = R;

    case 'DASCMOP7'
        R = UniformPoint(N,3);
        P = R;

    case 'DASCMOP8'
        R = UniformPoint(N,3);
        P = R;

    case 'DASCMOP9'
        R = UniformPoint(N,3);
        P = R;

    case 'LIRCMOP1'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - R(:,1).^2;
        P = R;

    case 'LIRCMOP2'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - sqrt(R(:,1));
        P = R;

    case 'LIRCMOP3'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - R(:,1).^2;
        P = R;

    case 'LIRCMOP4'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - sqrt(R(:,1));
        P = R;

    case 'LIRCMOP5'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - sqrt(R(:,1));
        P = R;

    case 'LIRCMOP6'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - R(:,1).^2;
        P = R;

    case 'LIRCMOP7'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - sqrt(R(:,1));
        P = R;

    case 'LIRCMOP8'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - sqrt(R(:,1));
        P = R;

    case 'LIRCMOP9'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - R(:,1).^2;
        P = R;

    case 'LIRCMOP10'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - sqrt(R(:,1));
        P = R;

    case 'LIRCMOP11'
        R = [1.3965,0.1591;1.0430,0.5127;0.6894,0.8662;
            0.3359,1.2198;0.0106,1.6016;0,2.1910;1.8730,0];
        P = R;

    case 'LIRCMOP12'
        R = [1.6794,0.4419;1.3258,0.7955;0.9723,1.1490;2.0320,0.0990;
            0.6187,1.5026;0.2652,1.8562;0,2.2580;2.5690,0];
        P = R;

    case 'LIRCMOP13'
        R = UniformPoint(N,3);
        P = R;

    case 'LIRCMOP14'
        R = UniformPoint(N,3);
        P = R;

    case 'FCP1'
        P = UniformPoint(N,M);
    
    case 'FCP1Mod'
        t = 0.5*pi*(0:1/N:1)';
        P = [cos(t),sin(t)];

    case 'Test1'
        P = UniformPoint(N,M);

    case 'Test2'
        t = 0.5*pi*(0:1/N:1)';
        P = [cos(t),sin(t)];

    case 'FCP2'
        t = 0:1/N:1;
        f1 = cos(0.5*pi*t);
        f2 = sin(0.5*pi*t)+0.2*sin(4*pi*t);
        P = [f1',f2'];
        P = P(find(NDSort(P,1)==1),:);

    case 'FCP3'
        t = 0.5*pi*(0:1/N:1)';
        P = [cos(t),sin(t)];

    case 'FCP4'
        t = 0:1/N:1;
        f1 = 1-t;
        f2 = t+0.2*sin(4*pi*t);
        P = [f1',f2'];
        P = P(find(NDSort(P,1)==1),:);

    case 'FCP5'
        P = [];

    otherwise
        P = [];
end
end


function [W,N] = UniformPoint(N,M,method)
%UniformPoint - Generate a set of uniformly distributed points.
%
%   [W,L] = UniformPoint(N,M) returns approximately N uniformly distributed
%   points with M objectives on the unit hyperplane via the normal-boundary
%   intersection method with two layers. Note that the number of sampled
%   points L may be slightly smaller than the predefined size N due to the
%   need for uniformity.
%
%   [W,L] = UniformPoint(N,M,'ILD') returns approximately N uniformly
%   distributed points with M objectives on the unit hyperplane via the
%   incremental lattice design. Note that the number of sampled points L
%   may be slightly larger than the predefined size N due to the need for
%   uniformity.
%
%   W = UniformPoint(N,M,'MUD') returns exactly N uniformly distributed
%   points with M objectives on the unit hyperplane via the mixture uniform
%   design method.
%
%   [W,L] = UniformPoint(N,M,'grid') returns approximately N uniformly
%   distributed points with M objectives in the unit hypercube via the grid
%   sampling. Note that the number of sampled points L may be slighly
%   larger than the predefined size N due to the need for uniformity.
%
%   W = UniformPoint(N,M,'Latin') returns exactly N randomly distributed
%   points with M objectives in the unit hypercube via the Latin hypercube
%   sampling method.
%
%   Example:
%       [W,N] = UniformPoint(275,10)
%       [W,N] = UniformPoint(286,10,'ILD')
%       [W,N] = UniformPoint(102,10,'MUD')
%       [W,N] = UniformPoint(1000,3,'grid')
%       [W,N] = UniformPoint(103,10,'Latin')

%------------------------------- Reference --------------------------------
% [1] Y. Tian, X. Xiang, X. Zhang, R. Cheng, and Y. Jin, Sampling reference
% points on the Pareto fronts of benchmark multi-objective optimization
% problems, Proceedings of the IEEE Congress on Evolutionary Computation,
% 2018.
% [2] T. Takagi, K. Takadama, and H. Sato, Incremental lattice design
% of weight vector set, Proceedings of the Genetic and Evolutionary
% Computation Conference Companion, 2020, 1486-1494.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

if nargin < 3
    method = 'NBI';
end
[W,N] = feval(method,N,M);
end

function [W,N] = NBI(N,M)
H1 = 1;
while nchoosek(H1+M,M-1) <= N
    H1 = H1 + 1;
end
W = nchoosek(1:H1+M-1,M-1) - repmat(0:M-2,nchoosek(H1+M-1,M-1),1) - 1;
W = ([W,zeros(size(W,1),1)+H1]-[zeros(size(W,1),1),W])/H1;
if H1 < M
    H2 = 0;
    while nchoosek(H1+M-1,M-1)+nchoosek(H2+M,M-1) <= N
        H2 = H2 + 1;
    end
    if H2 > 0
        W2 = nchoosek(1:H2+M-1,M-1) - repmat(0:M-2,nchoosek(H2+M-1,M-1),1) - 1;
        W2 = ([W2,zeros(size(W2,1),1)+H2]-[zeros(size(W2,1),1),W2])/H2;
        W  = [W;W2/2+1/(2*M)];
    end
end
W = max(W,1e-6);
N = size(W,1);
end

function [W,N] = ILD(N,M)
I = M * eye(M);
W = zeros(1,M);
edgeW = W;
while size(W) < N
    edgeW = repmat(edgeW,M,1) + repelem(I,size(edgeW,1),1);
    edgeW = unique(edgeW,'rows');
    edgeW(min(edgeW,[],2)~=0,:) = [];
    W = [W+1;edgeW];
end
W = W./sum(W,2);
W = max(W,1e-6);
N = size(W,1);
end

function [W,N] = MUD(N,M)
X = GoodLatticePoint(N,M-1).^(1./repmat(M-1:-1:1,N,1));
X = max(X,1e-6);
W = zeros(N,M);
W(:,1:end-1) = (1-X).*cumprod(X,2)./X;
W(:,end)     = prod(X,2);
end

function [W,N] = grid(N,M)
gap = linspace(0,1,ceil(N^(1/M)));
eval(sprintf('[%s]=ndgrid(gap);',sprintf('c%d,',1:M)))
eval(sprintf('W=[%s];',sprintf('c%d(:),',1:M)))
N = size(W,1);
end

function [W,N] = Latin(N,M)
[~,W] = sort(rand(N,M),1);
W = (rand(N,M)+W-1)/N;
end

function Data = GoodLatticePoint(N,M)
hm           = find(gcd(1:N,N)==1);
udt          = mod((1:N)'*hm,N);
udt(udt==0)  = N;
nCombination = nchoosek(length(hm),M);
if nCombination < 1e4
    Combination = nchoosek(1:length(hm),M);
    CD2 = zeros(nCombination,1);
    for i = 1 : nCombination
        UT     = udt(:,Combination(i,:));
        CD2(i) = CalCD2(UT);
    end
    [~,minIndex] = min(CD2);
    Data = udt(:,Combination(minIndex,:));
else
    CD2 = zeros(N,1);
    for i = 1 : N
        UT     = mod((1:N)'*i.^(0:M-1),N);
        CD2(i) = CalCD2(UT);
    end
    [~,minIndex] = min(CD2);
    Data = mod((1:N)'*minIndex.^(0:M-1),N);
    Data(Data==0) = N;
end
Data = (Data-1)/(N-1);
end

function CD2 = CalCD2(UT)
[N,S] = size(UT);
X     = (2*UT-1)/(2*N);
CS1 = sum(prod(2+abs(X-1/2)-(X-1/2).^2,2));
CS2 = zeros(N,1);
for i = 1 : N
    CS2(i) = sum(prod((1+1/2*abs(repmat(X(i,:),N,1)-1/2)+1/2*abs(X-1/2)-1/2*abs(repmat(X(i,:),N,1)-X)),2));
end
CS2 = sum(CS2);
CD2 = (13/12)^S-2^(1-S)/N*CS1+1/(N^2)*CS2;
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




