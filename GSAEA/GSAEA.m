%------------------------------- Copyright ------------------------------------
% Copyright (c) MDO Group, UNSW, Australia. You are free to use the GSAEA and SASSEA for
% research purposes. All publications which use this code should acknowledge the
% use of "SASSEA", "GSAEA" and references:
% "K.H.Rahi, H.K.Singh, T.Ray, A steady-state algorithm for solving expensive multi-objective
% optimization problems with non-parallelizable evaluations, IEEE Transactions on Evolutionary Computation, 2022" and
% "K.H.Rahi, H.K.Singh, T.Ray, A generalized surrogate-assisted evolutionary
% algorithm for expensive multi-objective Optimization, IEEE Congress on Evolutionary Computation, 2023 (Accepted for publication)".
%------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Multirun: This function is used to run specified algorithms and problems across multiple trials.
%  Set the parameters first at the Params.m file properly.
%  Check that matlab path is added properly for all functions required
%  Now run the 'Multirun.m' script.
%% Coded by
%  Kamrul Hasan Rahi
%  k.rahi@student.adfa.edu.au; kamrulhasanme038@gmail.com
%  Last updated: March 24, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GSAEA - Generalized Surrogate Assisted Evolutionary Algorithm
function GSAEA(path)
load('Parameters.mat');

%% Setting the random seed
rng(param.seed,'twister');

%% Loading the problem
prob = feval(param.prob_name);

%% Set the population size
param.popsize = min(11*prob.nx - 1,param.MFE/2);

%% Set total cost
tot_cost = param.MFE*prob.nf;

%% Initialize population
pop = Initialize_pop(prob,param.popsize);
pop.x(:,prob.cast) = round(pop.x(:,prob.cast));

%% Evaluate population
[pop.muf,~] = feval(param.prob_name,pop.x);

%% Setting up the Archive
Archive.x = pop.x;Archive.muf = pop.muf;Archive.sf = zeros(size(pop.muf));
Archive.cost = repmat((1:size(Archive.muf,1))',1,size(Archive.muf,2));
Archive.eval_counter = ones(size(Archive.muf));cost = sum(Archive.cost,2);curr_cost = cost(end);
Archive.shadow_nd = [];

%% Main body
while curr_cost <= tot_cost
    %% Infill sampling
    [candidate,Archive] = InfillSelection(Archive,param,prob);
    
    %% Modification
    SampleNo = min(param.iterate,size(candidate.x,1));
    SelCanF = [];new.x = [];count = 0;
    while count < SampleNo && ~isempty(candidate.x)
        [final_id] = Probablistic_selection(candidate,Archive,SelCanF);
        candidate.x = candidate.x(final_id,:);candidate.muf = candidate.muf(final_id,:);candidate.sf = candidate.sf(final_id,:);
        SelCanF = [SelCanF;candidate.muf(1,:)];
        new.x = [new.x;candidate.x(1,:)];count = count + 1;
        candidate.x(1,:) = [];candidate.muf(1,:) = [];candidate.sf(1,:) = [];
    end
    
    if ~isempty(candidate.x)
        %% Evaluate the selected candidate
        new.x(:,prob.cast) = round(new.x(:,prob.cast));
        [new.muf,~] = feval(param.prob_name,new.x);
        
        %% Update the Archive
        Archive.x(end+1:end+SampleNo,:) = new.x;Archive.muf(end+1:end+SampleNo,:) = new.muf;
        Archive.sf(end+1:end+SampleNo,:) = zeros(SampleNo,size(pop.muf,2));
    end
    
    %% Updating shadow ND list
    [shadow] = CreateShadowList(Archive);
    Archive.shadow_nd = shadow;
    
    %% Modification
    Archive.cost(end+1:end+SampleNo,:) = [Archive.cost(end,:)+repmat((1:SampleNo)',1,size(Archive.muf,2))];
    Archive.eval_counter = [Archive.eval_counter;repmat(ones(1,size(Archive.muf,2)),SampleNo,1)];
    cost = sum(Archive.cost,2);curr_cost = cost(end);
    
    %% Displaying the run state
    disp(strcat(path,filesep,'NFE - ',num2str(100*curr_cost/tot_cost)));
end

%% Storing relevant information
save(strcat(param.prob_name,'_Archive.mat'),'Archive');
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% THE END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Initializing a population
function pop = Initialize_pop(prob,no_solutions)
pop.x = repmat(prob.bounds(:,1)',no_solutions,1)+repmat((prob.bounds(:,2)-prob.bounds(:,1))',no_solutions,1).*lhsdesign(no_solutions,prob.nx);
return


%% Infill selection algorithm
function [candidate,Archive] = InfillSelection(Archive,param,prob)
%% Train the global model
if size(Archive.x,1) == param.popsize
    theta = 5.*ones(prob.nx,prob.nf);
else
    theta = Archive.theta;
end

%% Global model
id_sparse = Sparse_x(Archive.x,[],param,prob);
x_trg = Archive.x(id_sparse,:);
y_trg = Archive.muf(id_sparse,:);
[gprMdl,theta] = GPR_train(x_trg,y_trg,theta);
Archive.theta = theta;

%% Normalization bound identification
[ideal,nadir] = NormBounds(Archive.muf);
Archive.ideal = ideal;
Archive.nadir = nadir + ((nadir-ideal).*0.1);

%% Running the EA on the model
pop = Run_EA(prob,param,Archive,gprMdl);

%% Subset selection based on distance or probabilistic metric
pop.x(:,prob.cast) = round(pop.x(:,prob.cast));
[candidate.x,id] = unique(pop.x,'rows','stable');
candidate.muf = pop.muf(id,:);candidate.sf = pop.sf(id,:);
[~,id1] = intersect(candidate.x,Archive.x,'rows','stable');
candidate.x(id1,:) = [];candidate.muf(id1,:) = [];candidate.sf(id1,:) = [];
return



%% Sparse neighbourhood solutions
function [id] = Sparse_x(archive_x,x,param,prob)
n_Ax = (archive_x-repmat(prob.bounds(:,1)',size(archive_x,1),1))./(repmat(prob.bounds(:,2)'-prob.bounds(:,1)',size(archive_x,1),1));
if ~isempty(x)
    % Checking sparsity in archive
    list = [];
    for i = 1:size(archive_x,1)
        others = setdiff(1:size(archive_x,1),i);
        d = min(pdist2(n_Ax(i,:),n_Ax(others,:)));
        if d > param.x_threshold
            list = [list,i];
        end
    end
    new_Ax = n_Ax(list,:);
    % Checking sparsity among the neighbours of the given point
    n_x = (x-prob.bounds(:,1)')./(prob.bounds(:,2)'-prob.bounds(:,1)');
    d = pdist2(n_x,new_Ax);
    [~,id1] = sort(d);
    id = list(id1);
    if ~isempty(id) && length(id) >= min(size(archive_x,1),param.x_neighbour*size(archive_x,2))
        id = id(1:min(size(archive_x,1),param.x_neighbour*size(archive_x,2)));
    else
        id = [];
    end
else
    % Checking sparsity in archive
    list = [];
    for i = 1:size(archive_x,1)
        others = setdiff(1:size(archive_x,1),i);
        d = min(pdist2(n_Ax(i,:),n_Ax(others,:)));
        if d > param.x_threshold
            list = [list,i];
        end
    end
    id = list;
end
return



%% Training model
function [gpr,theta] = GPR_train(x,y,varargin)
theta = varargin{:};
for i = 1:size(y,2)
    id = find(~isnan(y(:,i)));
    xx = x(id,:);
    yy = y(id,i);
    [xx,id1] = unique(xx,'rows','stable');
    yy = yy(id1,:);
    dmodel = dacefit(xx,yy,'regpoly1','corrgauss',theta(:,i)',1e-5.*ones(1,size(x,2)),100.*ones(1,size(x,2)));
    gpr{i} = dmodel;
    theta(:,i) = dmodel.theta';
end
return


%% Identify normalization bounds
function [ideal,nadir] = NormBounds(f)
% Normalization
[id_fronts,~,~,~] = E_NDSort_c(f);
% Finding ND set and ideal, nadir point
ids = 0;count = 1;tmp_fnd = [];
while ~isempty(ids)
    tmp_fnd = [tmp_fnd;f(id_fronts{count},:)];
    ideal = min(tmp_fnd,[],1);nadir = max(tmp_fnd,[],1);
    dif = (nadir-ideal);
    ids = find(dif < 1e-4);
    count = count + 1;
end
return



%% Running sub-EA
function pop = Run_EA(prob,param,archive,gprMdl)
%% Initialize the population for surrogate assisted search
s_pop = Initialize_pop(prob,param.EA_popsize);
[ids,~,~,~] = E_NDSort_c(archive.muf);
if length(ids{1}) < param.EA_popsize
    nd_x = archive.x(ids{1},:);
    s_pop.x = [nd_x;s_pop.x];
end
[mu,sigma] = GPR_predict(gprMdl,s_pop.x,archive);
s_pop.muf = mu(:,1:prob.nf);s_pop.sf = sigma(:,1:prob.nf);c_pop = [];

%% Predicted archive
s_archive.x = s_pop.x;s_archive.muf = s_pop.muf;s_archive.sf = s_pop.sf;

%% Reduction
[s_pop] = Reduce(param,s_pop,archive);

for i = 1:param.EA_generations
    %% Offspring generation
    c_pop.x = Generate_child_SBX_PM(s_pop,prob,param,archive);
    
    %% Fitness evaluation
    [mu,sigma] = GPR_predict(gprMdl,c_pop.x,archive);
    c_pop.muf = mu(:,1:prob.nf);
    c_pop.sf = sigma(:,1:prob.nf);
    
    %% Predicted archive
    s_archive.x = [s_archive.x;c_pop.x];s_archive.muf = [s_archive.muf;c_pop.muf];
    s_archive.sf = [s_archive.sf;c_pop.sf];
    s_pop.x = [s_pop.x;c_pop.x];s_pop.muf = [s_pop.muf;c_pop.muf];
    s_pop.sf = [s_pop.sf;c_pop.sf];
    
    %% Environmental Selection
    [s_pop] = Reduce(param,s_pop,archive);
end

%% Final population
% [pop] = Reduce(param,s_archive,archive); % Final pop from whole archive
pop = s_pop; % Final pop from parent + latest child
return



%% Prediction
function [y_prd,s_prd] = GPR_predict(gpr,x_tst,archive)
y_prd = [];s_prd = [];
if isfield(archive,'mug')==1
    response = [archive.muf archive.mug];
else
    response = [archive.muf];
end
for i = 1:length(gpr)
    [tmp_mu,tmp_sigma] = DacePredict(x_tst, gpr{i});
    
    % This is a check if the Augmented model prediction is sought
    if(size(x_tst,2)==size(archive.x,2))
        %Replacing the predicted values by true values if they have been evaluated
        for j = 1:size(x_tst,1)
            tmp = repmat(x_tst(j,:),size(archive.x,1),1)-archive.x;
            d = sum(tmp.*tmp,2);
            id = find(d==0);r=nan(1,size(archive.muf,2));
            % These Archive solutions have some information
            if ~isempty(id)
                for k=1:length(id)
                    for p=1:size(archive.muf,2)
                        if(~isnan(response(id(k),p)))
                            r(p)=response(id(k),p);
                        end
                    end
                end
                if(~isnan(r(i)))
                    tmp_mu(j) = r(i);
                    tmp_sigma(j) = 0;
                end
            end
        end
    end
    y_prd = [y_prd tmp_mu];
    s_prd = [s_prd tmp_sigma];
end
return



%% This function identified population for the next generation
function [pop] = Reduce(param,s_archive,archive)
N = param.EA_popsize;
% Initial ranking
[rank] = RankingDMirror(s_archive.muf,s_archive.sf,param.spacing,archive);

% Environmental selection
if length(rank) >= N
    pop.x = s_archive.x(rank(1:N),:);
    pop.muf = s_archive.muf(rank(1:N),:);
    pop.sf = s_archive.sf(rank(1:N),:);
else
    if mod(length(rank),2) == 1
        id = setdiff((1:size(s_archive.muf,1))',rank,'stable');
        rank = [rank;id(randperm(length(id),1))];
    end
    pop.x = s_archive.x(rank,:);
    pop.muf = s_archive.muf(rank,:);
    pop.sf = s_archive.sf(rank,:);
end
return



%% Simulated binary crossover and polynomial mutation
function [X_Child] = Generate_child_SBX_PM(pop,prob,param,archive)
LB = prob.bounds(:,1)';UB = prob.bounds(:,2)';
X_Child = genetic_operator(LB,UB,param,prob,pop.x);
% Uniqueness test of Offspring
[X_Child] = validity(archive.x,X_Child,prob,LB,UB);
return

function [X_child] = genetic_operator(LB,UB,def,prob,X_parent_ordered)
[X_child] = crossover_SBX_matrix(LB,UB,def,prob,X_parent_ordered);
[X_child] = mutation_POLY_matrix(LB,UB,def,prob,X_child);
return

% Simulated Binary Crossover (SBX) operator
function [y1, y2] = op_SBX_matrix(l_limit,u_limit, x1, x2, eta)
y1 = x1;
y2 = x2;
ipos = find(abs(x1-x2) > 1e-6);
if ~isempty(ipos)
    x1_op = x1(ipos);
    x2_op = x2(ipos);
    l_op = l_limit(ipos);
    u_op = u_limit(ipos);
    pos_swap = x2_op < x1_op;
    tmp = x1_op(pos_swap);
    x1_op(pos_swap) = x2_op(pos_swap);
    x2_op(pos_swap) = tmp;
    r = rand(size(x1_op));
    beta = 1 + (2*(x1_op - l_op)./(x2_op - x1_op));
    alpha = 2 - beta.^-(eta+1);
    betaq = (1./(2-r.*alpha)).^(1/(eta+1));
    betaq(r <= 1./alpha) = (r(r <= 1./alpha).*alpha(r <= 1./alpha)).^(1/(eta+1));
    y1_op = 0.5 * (x1_op + x2_op - betaq.*(x2_op - x1_op));
    
    beta = 1 + 2*(u_op - x2_op)./(x2_op - x1_op);
    alpha = 2 - beta.^-(eta+1);
    betaq = (1./(2-r.*alpha)).^(1/(eta+1));
    betaq(r <= 1./alpha) = (r(r <= 1./alpha).*alpha(r <= 1./alpha)).^(1/(eta+1));
    y2_op = 0.5 * (x1_op + x2_op + betaq.*(x2_op - x1_op));
    
    y1_op(y1_op < l_op) = l_op(y1_op < l_op);
    y1_op(y1_op > u_op) = u_op(y1_op > u_op);
    
    y2_op(y2_op < l_op) = l_op(y2_op < l_op);
    y2_op(y2_op > u_op) = u_op(y2_op > u_op);
    
    pos_swap = (rand(size(x1_op)) <= 0.5);
    tmp = y1_op(pos_swap);
    y1_op(pos_swap) = y2_op(pos_swap);
    y2_op(pos_swap) = tmp;
    
    y1(ipos) = y1_op;
    y2(ipos) = y2_op;
end
return

function [c, fn_evals] = crossover_SBX_matrix(LB,UB,def,prob,p)
eta=def.distribution_mutation;
c = p; % parent size = no_solution*no_variable.
fn_evals = 0;
A = rand(size(p,1)/2,1);
is_crossover =[(A <= def.prob_crossover)';(A <= def.prob_crossover)'];
p_cross = p(is_crossover,:);
[N,~] = size(p_cross);
c_cross = p_cross;
p1_cross = p_cross(1:2:N,:);
p2_cross = p_cross(2:2:N,:);
B = rand(size(p_cross,1)/2,prob.nx);
l_limit = repmat(LB,size(p_cross,1)/2,1);
u_limit = repmat(UB,size(p_cross,1)/2,1);
cross_pos = (B <= 0.5);
l_cross = l_limit(cross_pos);
u_cross = u_limit(cross_pos);
p1 = p1_cross(cross_pos);
p2 = p2_cross(cross_pos);
c1 = p1_cross;
c2 = p2_cross;
[y1, y2] = op_SBX_matrix(l_cross,u_cross,p1,p2,eta);
c1(cross_pos) = y1;
c2(cross_pos) = y2;
c_cross(1:2:N,:) = c1;
c_cross(2:2:N,:) = c2;
c(is_crossover,:) = c_cross;
return

% Polynomial mutation operator: Matrix Form
function [x] = op_POLY_matrix(LB,UB,x,def)
def.distribution_mutation;
x_min = LB;
x_max = UB;
pos_mu = find(x_max > x_min);
if ~isempty(pos_mu)
    x_mu = x(pos_mu);
    x_min_mu = x_min(pos_mu);
    x_max_mu = x_max(pos_mu);
    delta1 = (x_mu - x_min_mu)./(x_max_mu - x_min_mu);
    delta2 = (x_max_mu - x_mu)./(x_max_mu - x_min_mu);
    mut_pow = 1/(def.distribution_mutation+1);
    rand_mu = rand(size(delta2));
    xy = 1 - delta2;
    val = 2*(1 - rand_mu) + 2*(rand_mu - 0.5).*xy.^(def.distribution_mutation+1);
    deltaq = 1 - val.^mut_pow;
    xy(rand_mu <= 0.5) = 1 - delta1(rand_mu <= 0.5);
    val(rand_mu <= 0.5) = 2*rand_mu(rand_mu <= 0.5) + (1-2*rand_mu(rand_mu <= 0.5)).* xy(rand_mu <= 0.5).^(def.distribution_mutation+1);
    deltaq(rand_mu <= 0.5) = val(rand_mu <= 0.5).^mut_pow - 1;
    
    x_mu = x_mu + deltaq.*(x_max_mu - x_min_mu);
    x_mu(x_mu < x_min_mu) = x_min_mu(x_mu < x_min_mu);
    x_mu(x_mu > x_max_mu) = x_max_mu(x_mu > x_max_mu);
    
    x(pos_mu) = x_mu;
end
return

function [p, fn_evals] = mutation_POLY_matrix(LB,UB,def,prob, p)
fn_evals = 0;
A = rand(size(p,1),prob.nx);
l_limit = repmat(LB,size(p,1),1);
u_limit = repmat(UB,size(p,1),1);
p_mut = p(A <= def.prob_mutation);
l_mut = l_limit(A <= def.prob_mutation);
u_mut = u_limit(A <= def.prob_mutation);
p_mut = op_POLY_matrix(l_mut,u_mut,p_mut,def);
p(A <= def.prob_mutation) = p_mut;
return

%% Checking similarity with archive population with current child pop.
function [v_X_child] = validity(X_pop,X_child,prob,LB,UB)
N = size(X_child,1);
if ~isempty(X_child)
    X_child_nor = (X_child-repmat(LB,size(X_child,1),1))./repmat((UB-LB),size(X_child,1),1); % precision limit upto 5 digit with normalized solution
    X_child_r = round(X_child_nor,4);
    X_pop_nor = (X_pop-repmat(LB,size(X_pop,1),1))./repmat((UB-LB),size(X_pop,1),1);
    X_pop_r = round(X_pop_nor,4);
    [a,aid] = unique(X_child_r,'rows','stable');
    [b,b_id] = setdiff(a,X_pop_r,'rows','stable');
    if size(b,1)~=size(X_child,1)
        X_child_n = X_child(aid(b_id),:);
        N = N - size(X_child_n,1);
        new_pop = repmat(LB,N,1)+repmat((UB-LB),N,1).*lhsdesign(N,prob.nx);
        v_X_child = [X_child_n;new_pop];
    else
        v_X_child = X_child;
    end
else
    v_X_child = [];
end
return



%% Pulling from mirror points of reference points in normalized space
function [order] = RankingDMirror(mu,sigma,spacing,archive)
w = Direction_vector(size(mu,2),spacing);

%% Modification
ideal = archive.ideal;
nadir = archive.nadir;

%% Getting the solutions in scaled space
mu_N = (mu-repmat(ideal,size(mu,1),1))./repmat(nadir-ideal,size(mu,1),1);
factor = (nadir-ideal).^2;
sigma_N = sigma./repmat(factor,size(mu,1),1);

%% Ordering
order = assignmentmirror(mu_N,sigma_N,w);
return


%% Solution assignement to reference direction
function rank = assignmentmirror(ff,sigma,w)
w_p = w - 1;N = size(w,1);
rank = [];reference = (1:size(w,1))';
id_all = (1:size(ff,1))';
while ~isempty(reference) && ~isempty(id_all) && N > 0
    Angle = [];
    for i = 1:length(reference)
        tmp_w = w(reference(i),:) - w_p(reference(i),:);
        tmp_f = (ff(id_all,:) - w_p(reference(i),:));
        Angle(:,i)  = acos(1-pdist2(tmp_f,tmp_w,'cosine'));
    end
    [~,associate] = min(Angle,[],2);
    tmp = [];
    for i = 1:length(reference)
        id = find(associate == i);
        id_tmp = setdiff(id_all(id),rank,'stable');
        if ~isempty(id_tmp)
            t_w = w(reference(i),:) - w_p(reference(i),:);
            t_f = (ff(id_tmp,:) - w_p(reference(i),:));
            [mu_p,sigma_p] = Multivariate_Projection(t_w,t_f,sigma(id_tmp,:));
            score = CalculateScore(mu_p,sigma_p);
            [~,id_score] = max(score);
            rank = [rank;id_tmp(id_score)];
            N = N - 1;
        else
            tmp = [tmp,i];
        end
    end
    id_all = setdiff(id_all,rank,'stable');
    reference = reference(tmp);
end
return


%% Calculate round-robin for a single column with matrix formulation
function mean_score = CalculateScore(fval,fsigma)
ff = repmat(fval',size(fval,1),1);
ssigma = repmat(fsigma',size(fval,1),1);
s = zeros(size(ff));
b_f = repmat(fval,1,size(fval,1));
b_sigma = repmat(fsigma,1,size(fval,1));
denom = (sqrt(2.*(b_sigma.^2+ssigma.^2)));
id1 = find(denom ~= 0);
if ~isempty(id1)
    s(id1) = 0.5+0.5.*erf((ff(id1)-b_f(id1))./(sqrt(2.*(b_sigma(id1).^2+ssigma(id1).^2))));
    id2 = setdiff(1:size(s,1)*size(s,2),id1,'stable');
    if ~isempty(id2)
        id3 = find(b_f(id2) <= ff(id2));
        s(id2(id3)) = 1; % Violation of i is less than j: Winner is i
    end
else
    id2 = find(b_f <= ff);
    s(id2) = 1;
end
mean_score = sum(s,2)./size(s,2);
return


%% Multivariate Projection into Univariate Gaussian
function [mu_projected,sigma_projected] = Multivariate_Projection(V,mu,sigma)
V = V./sqrt(sum(V.^2,2));
mu_trans = mu';
for k = 1:size(mu,1)
    sigma_trans{k} = eye(size(mu,2)).*sigma(k,:);
end

for j = 1:size(mu,1)
    for k = 1:size(V,1)
        mu_projected(j,k)= V(k,:)*mu_trans(:,j);
        sigma_projected(j,k) = V(k,:)*sigma_trans{j}*V(k,:)';
    end
end
return


%% Probabilistic subset selection
function [final_id] = Probablistic_selection(candidate,Archive,SelCanF)
if ~isempty(candidate)
    mu = candidate.muf;sigma = candidate.sf;
    sigma(sigma < 1e-4) = 1e-4;
    F = Archive.muf;IDall = 1:size(mu,1);
    [id_fronts,~,~,~] = E_NDSort_c(F);
    nd_ids = id_fronts{1};
    fnd = F(nd_ids,:);
    
    %% Setting the normalization bounds
    ideal = Archive.ideal;nadir = Archive.nadir;
    
    %% Normalization
    mu_N = (mu-repmat(ideal,size(mu,1),1))./repmat(nadir-ideal,size(mu,1),1);
    factor = (nadir-ideal).^2;
    sigma_N = sigma./repmat(factor,size(mu,1),1);
    
    %% Non-dominated set of predicted solution
    IL = [];
    for i = 1:size(mu,1)
        tmu = [mu(i,:);fnd];
        [fronts,~,~,~] = E_NDSort_c(tmu);
        if ismember(1,fronts{1})
            IL = [IL;i];
        end
    end
    if isempty(IL)
        [fronts,~,~,~] = E_NDSort_c(mu);
        IL = fronts{1};
    end
    
    %% Final candidates
    tmp_F = [fnd;Archive.shadow_nd;SelCanF];
    tmp_FN = (tmp_F-repmat(ideal,size(tmp_F,1),1))./repmat(nadir-ideal,size(tmp_F,1),1);
    if ~isempty(IL)
        for i = 1:length(IL)
            for j = 1:size(tmp_FN,1)
                ref = tmp_FN(j,:);
                gm = gmdistribution(mu_N(IL(i),:),sigma_N(IL(i),:));
                dist(i,j) = mahal(gm,ref);
            end
            distance(i) = min(dist(i,:));
        end
        [~,id_dis] = sort(distance,'descend');
        final_id = IDall(IL(id_dis));
    end
else
    error('There is no candidate.');
end
return



%% Creating shadow ND list
function [shadow] = CreateShadowList(Archive)
% Check whether the newly evaluated soln is dominated or not
F = Archive.muf;shadow = Archive.shadow_nd;
[id_fronts,~,~,~] = E_NDSort_c(F);
% Normalization bound identification
[ideal,nadir] = NormBounds(F);
nadir = nadir + ((nadir-ideal).*0.1);
F_N = (F-ideal)./(nadir-ideal);
ff = F_N(1:end-1,:); ff_new = F_N(end,:);
if ~ismember(size(F,1),id_fronts{1})
    shadow = [shadow;F(end,:)];
    % Check whether the closest archived solution is dominated or not
    dis = sum((repmat(ff_new,size(ff,1),1) - ff).^2,2);
    [~,id_min] = min(dis);
    if ~ismember(id_min,id_fronts{1})
        shadow = [shadow;F(id_min,:)];
    end
else
    % Check whether the closest archived solution is dominated or not
    dis = sum((repmat(ff_new,size(ff,1),1) - ff).^2,2);
    [~,id_min] = min(dis);
    if ~ismember(id_min,id_fronts{1})
        shadow = [shadow;F(id_min,:)];
    end
end
shadow = unique(shadow,'rows','stable');
return



function [wall]=Direction_vector(M,pall)
count = 1;
wall = [];
for k = 1:numel(pall)
    p = pall(k);
    NumPoints=zeros(1,M-1);
    for i=1:(M-1)
        NumPoints(i)=p;
    end
    % Partition the first objective
    Beta1 = [0:(NumPoints(1))]'/(NumPoints(1));
    % Save the previous values
    Beta = Beta1;
    for i = 2:M-1
        % Compute the combination i.e. p*(p-1)*(p-2)*,-----,*0
        ki = round((1-sum(Beta1,2))*(NumPoints(i)));
        Beta = [];
        for j =1:size(Beta1,1)
            % Compute each subvector of (0,1,...p)/p,(0,1,...p-1)/p,...
            BetaVec = [0:ki(j)]'/(NumPoints(i));
            numOfreplications = length(BetaVec); % identify the length
            % Replicate each of the previous values in the equal size to all the subvectors
            Beta = [Beta; [repmat(Beta1(j,:), numOfreplications,1) BetaVec] ];
        end
        Beta1 = Beta;
    end
    % Compute the last objective values
    BetaVec = 1 - sum(Beta1,2);
    w= [Beta BetaVec];%include the last objective values
    w = w*count + (1-count)/M;
    wall = [wall;w];
    count = count/2;
end
return