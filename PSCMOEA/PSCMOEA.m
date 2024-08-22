%% PSCMOEA: Constrained Multiobjective Evolutionary Algorithm based on Probabilistic selection for expensive CMOPs.
% Normalization: Bounds are created based on true feasible ND and true infeasible ND which are ND to all feasible ND (when archive has mixed solutions), same as SASSEA (when Whole archive is feasible), ideal and nadir are min and max of all F in the Archive when whole archive is Infeasible;
% SubEA (Decomposition based): Unconstrained search and Constrained search; Initialization - random + Pressure given by IDEA ranked true archived solution, ranking - decomposed along mirrored RVs, subgroup ranking is done by lexicographic ranking based on PF followed by PD;
% Infill selection: Selection of RPs, true feasible ND and true infeasible ND which are ND to all feasible ND (when archive has mixed solutions), same as SASSEA (when Whole archive is feasible), When the whole archive is infeasible, select the top candidate solution with lexicographic ranking based on PF followed by PD + keep the track of the reference directions untill bound/all the reference vectors are finished;
% Shadow mechanism: Shadow mechanism is not applied until one feasible solution found;

function PSCMOEA(path)
load('Parameters.mat');

%% Setting the random seed
rng(param.seed,'twister');

%% Loading the problem
prob = feval(param.prob_name);

%% Set the population size
param.popsize = 11*prob.nx - 1;
param.ConstrFlag = 2; % Search considering the constraints 

%% Set total evaluation budget and current evaluation
TFE = param.MFE;FE = param.popsize;

%% Initialize population
pop = Initialize_pop(prob,param.popsize);

%% Evaluate population
[pop.muf,pop.mug] = feval(param.prob_name,pop.x);

%% Setting up the Archive
Archive.x = pop.x;Archive.muf = pop.muf;Archive.sf = zeros(size(pop.muf));
if prob.ng > 0
    Archive.mug = pop.mug;Archive.sg = zeros(size(pop.mug));
end
Archive.SubEAInitFlag = 1;
Archive.SearchFlag = param.ConstrFlag;
Archive.boundtrack = [];

%% RV tracking initialization: When the archive is fully infeasible
cv = sum(max(Archive.mug,0),2);id = find(cv == 0);
if isempty(id)
    Archive.flag = 1;
    Archive.ideal = []; Archive.nadir = [];
    Archive.RVs = [];
    w = Direction_vector(prob.nf,param.spacing);
    Archive.NumRVs = size(w,1);
else
    Archive.flag = 0;
    Archive.RVs = [];
end

%% Find the correlation between F and CV of the infeasible archive
idinfeas = setdiff(1:size(Archive.x,1),id,'stable');
if ~isempty(idinfeas)
    CorPop.muf = Archive.muf(idinfeas,:);
    CorPop.mug = Archive.mug(idinfeas,:);
    [FRank,GRank] = FGRank(CorPop);
    RankCor = corr(FRank,GRank,'type','Kendall');
    Archive.RankCor = RankCor;
else
    Archive.RankCor = 0;
end

%% Intiate the shadow ND: when the archive has feasible sol, keep the infeasible sols in the shadow ND which are ND to the feasible ND; keep it empty when infeasible archive
if ~isempty(id)
    [id_fronts,~,~,~] = E_NDSort_c(Archive.muf);
    count = 0;
    for i = 1:length(id_fronts)
        frontids = id_fronts{i};
        if sum(ismember(frontids,id)) ~= 0
            count = i;
        end
    end
    shadow_ids = [];
    for i = 1:count
        frontids = id_fronts{i};
        logvec = find(ismember(frontids,id)==0);
        shadow_ids = [shadow_ids;frontids(logvec)];
    end
    Archive.shadow_nd = Archive.muf(shadow_ids,:);
else
    Archive.shadow_nd = [];
end

%% Visualization data generation
if param.visualization == 1
    SplitStr = split(param.prob_name,"_");
    prob = feval(param.prob_name);
    [P] = PF(SplitStr{1},1000,prob.nf);
    [RR] = PFUnCon(SplitStr{1},1000,prob.nf);
    R = FR(param.prob_name,prob.nf);
    if ~isempty(R)
        xx = R{1};yy = R{2};zz = R{3};
    else
        xx = [];yy = []; zz = [];
    end
    idinfeas = setdiff(1:size(Archive.x,1),id,'stable');
    f_feas = Archive.muf(id,:); f_infeas = Archive.muf(idinfeas,:);
end

%% Main body: Optimization loop starts
while FE < TFE
    %% Visualization
    if param.visualization == 1
        s1 = surf(xx,yy,zz);hold on;
        s1.EdgeColor = 'none';s1.FaceAlpha = 0.7;
        if prob.nf < 3
            view(0,90);
            plot(RR(:,1),RR(:,2),'k-');
            plot(P(:,1),P(:,2),'m.');
            plot(f_infeas(:,1),f_infeas(:,2),'k.');
            plot(f_feas(:,1),f_feas(:,2),'r*');
        else
            view(47,19);
            plot3(RR(:,1),RR(:,2),RR(:,3),'k.');
            plot3(P(:,1),P(:,2),P(:,3),'m.');
            plot3(f_infeas(:,1),f_infeas(:,2),f_infeas(:,3),'k.');
            plot3(f_feas(:,1),f_feas(:,2),f_feas(:,3),'r*');
        end
    end

    %% Infill sampling
    [candidate,Archive] = InfillSelection(Archive,param,prob);

    %% New Infill evaluation
    if ~isempty(candidate.x)
        %% Evaluate the selected candidate
        new.x = candidate.x(1,:);
        [new.muf,new.mug] = feval(param.prob_name,new.x);

        %% Update the Archive
        Archive.x(end+1,:) = new.x;Archive.muf(end+1,:) = new.muf;
        Archive.sf(end+1,:) = zeros(1,size(pop.muf,2));
        if prob.ng > 0
            Archive.mug(end+1,:) = new.mug;Archive.sg(end+1,:) = zeros(1,size(Archive.sg,2));
        end
        Archive.SubEAInitFlag = 1;
        FE = FE + 1; % Update the evaluation number

        %% Search switch flag
        flag = SwitchSearchMechanism(Archive,param);
        param.ConstrFlag = flag;
        Archive.SearchFlag = [Archive.SearchFlag;param.ConstrFlag];
     
        %% Updating shadow ND list
        cv = sum(max(Archive.mug,0),2);id = find(cv == 0);
        if ~isempty(id)
            [shadow] = CreateShadowList(Archive);
            Archive.shadow_nd = shadow;
            Archive.flag = 0;Archive.RVs = [];
        end

        %% Update the rank correlation
        idinfeas = setdiff(1:size(Archive.x,1),id,'stable');
        if ~isempty(idinfeas)
            CorPop.muf = Archive.muf(idinfeas,:);
            CorPop.mug = Archive.mug(idinfeas,:);
            [FRank,GRank] = FGRank(CorPop);
            RankCor = corr(FRank,GRank,'type','Kendall');
            Archive.RankCor = RankCor;
        else
            Archive.RankCor = 0;
        end
    else
        Archive.SubEAInitFlag = 0;
    end

    %% Displaying the run state
    disp(strcat(path,filesep,'NFE - ',num2str(100*FE/TFE)));

    %% Visualization
    if param.visualization == 1
        idinfeas = setdiff(1:size(Archive.x,1),id,'stable');
        f_feas = Archive.muf(id,:);f_infeas = Archive.muf(idinfeas,:);
        if prob.nf < 3
            plot(new.muf(:,1),new.muf(:,2),'go','LineWidth',1.5);
        else
            plot3(new.muf(:,1),new.muf(:,2),new.muf(:,3),'go','LineWidth',1.5);
        end
        hold off;
        drawnow;
    end
end

%% Visualization
if param.visualization == 1
    cv = sum(max(Archive.mug,0),2);id = find(cv == 0);
    idinfeas = setdiff(1:size(Archive.x,1),id,'stable');
    f_feas = Archive.muf(id,:);f_infeas = Archive.muf(idinfeas,:);
    s1 = surf(xx,yy,zz);hold on;
    s1.EdgeColor = 'none';s1.FaceAlpha = 0.7;
    if prob.nf < 3
        view(0,90);
        plot(RR(:,1),RR(:,2),'k-');
        plot(P(:,1),P(:,2),'m.');
        plot(f_feas(:,1),f_feas(:,2),'o','MarkerEdgeColor','b','MarkerFaceColor','b');
        plot(f_infeas(:,1),f_infeas(:,2),'k.');
    else
        view(47,19);
        plot3(RR(:,1),RR(:,2),RR(:,3),'k.');
        plot3(P(:,1),P(:,2),P(:,3),'m.');
        plot3(f_feas(:,1),f_feas(:,2),f_feas(:,3),'o','MarkerEdgeColor','b','MarkerFaceColor','b');
        plot3(f_infeas(:,1),f_infeas(:,2),f_infeas(:,3),'k.');
    end
    hold off;
    % drawnow;
    filename = strcat(param.prob_name,'_Archive_plot');
    saveas(gcf,filename);
    print(filename,'-depsc2','-r300');
    saveas(gcf,strcat(filename,'.jpg'));
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
    theta = 5.*ones(prob.nx,prob.nf+prob.ng);
else
    theta = Archive.theta;
end

%% Global model
id_sparse = Sparse_x(Archive.x,[],param,prob);
x_trg = Archive.x(id_sparse,:);
y_trg = Archive.muf(id_sparse,:);
if prob.ng > 0
    y_trg = [y_trg,Archive.mug(id_sparse,:)];
end
[gprMdl,theta] = GPR_train(x_trg,y_trg,theta,prob);
Archive.theta = theta;

%% Normalization bound identification (only feasible, mix of feasible and infeasible, only infeasible)
[ideal,nadir] = BoundIdentification(Archive);
if Archive.flag == 1 && ~isempty(Archive.ideal) && ~isempty(Archive.nadir) % Check whether RV tracking is continued or not
    tmp = pdist2(ideal,Archive.ideal) + pdist2(nadir,Archive.nadir);
    if tmp == 0 % Check whether norm bound is unchanged or not
        Archive.flag = 1;
    else
        Archive.flag = 0;Archive.RVs = [];
    end
end
Archive.ideal = ideal;Archive.nadir = nadir;

%% Visualization
if param.visualization == 1
    plot_data = [Archive.ideal,Archive.nadir];
    for i=1:size(plot_data,1)
        a = plot_data(i,1);b = plot_data(i,2);c = plot_data(i,3);d = plot_data(i,4);
        pgon = polyshape([a a c c],[b d d b]);
        pg = plot(pgon);
        pg.FaceColor = [0.50,0.50,0.50];
        pg.FaceAlpha = 0;
        pg.EdgeColor = [0,0,0];
        pg.LineWidth = 1.2;
    end
end

%% Visualization
if param.visualization == 1
%     for i = 1:size(pop.mug,1)
%         [ProdPF(i,:),IndPF(i,:)] = Prob_feas(pop.mug(i,:),pop.sg(i,:));
%     end
%     if size(Archive.muf,2) < 3
%         colormap('jet');
%         scatter(pop.muf(:,1),pop.muf(:,2),20,ProdPF,'filled');
%     else
%         scatter3(pop.muf(:,1),pop.muf(:,2),pop.muf(:,3),20,ProdPF,'filled');
%     end
    colors = {'k'};
    w = Direction_vector(prob.nf,param.spacing);
    w_p = w - 1;
    w_p = Archive.ideal + (Archive.nadir - Archive.ideal).*w_p;
    w2 = w + 1;
    w2 = Archive.ideal + (Archive.nadir - Archive.ideal).*w2;
    for j = 1:size(w,1)
        line([w_p(j,1) w2(j,1)],[w_p(j,2) w2(j,2)],'Color',colors{1},'linestyle',':','linewidth',1);
    end
end

%% Running the EAs based on the flag: Unconstrained and Constrained search
pop = Run_EA(prob,param,Archive,gprMdl);

%% Subset selection based on distance or probabilistic metric
candidate.x = [pop.x];
candidate.muf = [pop.muf];candidate.sf = [pop.sf];
candidate.mug = [pop.mug];candidate.sg = [pop.sg];
candidate.RVs = pop.RVs;
[candidate.x,id] = unique(candidate.x,'rows','stable');
candidate.muf = candidate.muf(id,:);candidate.sf = candidate.sf(id,:);
if isfield(candidate,'mug') == 1
    candidate.mug = candidate.mug(id,:);candidate.sg = candidate.sg(id,:);candidate.RVs = candidate.RVs(id);
end
[~,id1] = intersect(candidate.x,Archive.x,'rows','stable');
candidate.x(id1,:) = [];candidate.muf(id1,:) = [];candidate.sf(id1,:) = [];
if isfield(candidate,'mug') == 1
    candidate.mug(id1,:) = [];candidate.sg(id1,:) = [];candidate.RVs(id1) = [];
end

if param.visualization == 1
    if size(Archive.muf,2) < 3
        colormap('jet');
        plot(candidate.muf(:,1),candidate.muf(:,2),'bo');
    else
        plot3(candidate.muf(:,1),candidate.muf(:,2),candidate.muf(:,3),'bo');
    end
end

%% Infill selection
if ~isempty(candidate.x)
    [final_id,Archive] = Probablistic_selection(candidate,Archive,prob,param);
    candidate.x = candidate.x(final_id,:);candidate.muf = candidate.muf(final_id,:);candidate.sf = candidate.sf(final_id,:);
    if isfield(Archive,'mug') == 1
        candidate.mug = candidate.mug(final_id,:);candidate.sg = candidate.sg(final_id,:);
    end
end
return


%% Normalization bound identification: For full infeasible, feasible and mix of feasible/infeasible archive
function [ideal,nadir] = BoundIdentification(Archive)
cv = sum(max(Archive.mug,0),2);idfeas = find(cv == 0);
idinfeas = setdiff(1:size(Archive.x,1),idfeas,'stable');
if isempty(idfeas)
    % If full infeasible; full F space is considered for search - min and max of F 
    ideal = min(Archive.muf,[],1);
    nadir = max(Archive.muf,[],1);
elseif isempty(idinfeas)
    % If full feasible; bound is constructed based on the feasible ND solutions - extension flexibility is present
    f_feas = Archive.muf(idfeas,:);
    [ideal,nadir] = NormBounds(f_feas);
    nadir = nadir + ((nadir-ideal).*0.1);
else
    % If mix of feasible/infeasible in the archive; take feasible ND and infeasible solutions which are ND to at least one feasible ND and
    % take their min and max - bound extension flexibility is present
    FeasF = Archive.muf(idfeas,:);InFeasF = Archive.muf(idinfeas,:);
    [id_fronts1,~,~,~] = E_NDSort_c(FeasF);
    Fnd = FeasF(id_fronts1{1},:);
    IL = [];
    for i = 1:size(InFeasF,1)
        tmu = [InFeasF(i,:);Fnd];
        [fronts,~,~,~] = E_NDSort_c(tmu);
        if ismember(1,fronts{1}) && length(fronts{1}) > 1 % size(Fnd,1)
            IL = [IL;i];
        end
    end
    MixedFF = [Fnd;InFeasF(IL,:)];
    if size(MixedFF,1) == 1
        [id_fronts2,~,~,~] = E_NDSort_c(InFeasF);
        InFnd = InFeasF(id_fronts2{1},:);
        MixedFF = [Fnd;InFnd];
    end
    ideal = min(MixedFF,[],1);
    nadir = max(MixedFF,[],1);
    nadir = nadir + ((nadir-ideal).*0.1);
    % When the bound in any direction becomes almost a point
    dif = (nadir - ideal);
    ids = find(dif < 1e-4);
    if ~isempty(ids)
        minF = min(Archive.muf,[],1);
        maxF = max(Archive.muf,[],1);
        ext = nadir + ((maxF-minF).*0.1);
        nadir(ids) = ext(ids);
        Archive.boundtrack = [Archive.boundtrack,size(Archive.muf,1)];
    end
end
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
theta = varargin{1};prob = varargin{2};LB = prob.bounds(:,1)'; UB = prob.bounds(:,2)';
for i = 1:size(y,2)
    id = find(~isnan(y(:,i)));
    xx = x(id,:);
    yy = y(id,i);
    xn = (xx-repmat(LB,size(xx,1),1))./repmat((UB-LB),size(xx,1),1); % precision limit upto 5 digit with normalized solution
    xr = round(xn,6);
    [~,id1] = unique(xr,'rows','stable');
    xx = xx(id1,:);
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
ids = 0;tmp_fnd = [];
for i = 1:length(id_fronts)
    tmp_fnd = [tmp_fnd;f(id_fronts{i},:)];
    ideal = min(tmp_fnd,[],1);nadir = max(tmp_fnd,[],1);
    dif = (nadir-ideal);
    ids = find(dif < 1e-4);
    if isempty(ids)
        break;
    end
end
if ~isempty(ids)
    nadir = ideal + 1e-4;
end
return



%% Running sub-EA
function pop = Run_EA(prob,param,archive,gprMdl)
%% Initialize the population for surrogate assisted search
s_pop = Initialize_pop(prob,param.EA_popsize); % Random + IDEA ranked true archived solutions for initialization
if archive.SubEAInitFlag == 1
    [rank] = IDEASort(archive, param);
    if length(rank) > param.EA_popsize
        xx = archive.x(rank(1:param.EA_popsize),:);
    else
        xx = archive.x(rank,:);
    end
    s_pop.x = [xx;s_pop.x];
end
[mu,sigma] = GPR_predict(gprMdl,s_pop.x,archive);
s_pop.muf = mu(:,1:prob.nf);s_pop.sf = sigma(:,1:prob.nf);
if prob.ng > 0
    s_pop.mug = mu(:,(prob.nf+1):(prob.nf+prob.ng));s_pop.sg = sigma(:,(prob.nf+1):(prob.nf+prob.ng));
end
c_pop = [];

%% Predicted archive
s_archive.x = s_pop.x;s_archive.muf = s_pop.muf;s_archive.sf = s_pop.sf;
if prob.ng > 0
    s_archive.mug = s_pop.mug;s_archive.sg = s_pop.sg;
end

%% Reduction
count = 0;
[s_pop] = Reduce(param,s_pop,archive,count);

for i = 1:param.EA_generations
    %% Offspring generation
    c_pop.x = Generate_child_SBX_PM(s_pop,prob,param,archive);

    %% Fitness evaluation
    [mu,sigma] = GPR_predict(gprMdl,c_pop.x,archive);
    c_pop.muf = mu(:,1:prob.nf);c_pop.sf = sigma(:,1:prob.nf);
    if prob.ng > 0
        c_pop.mug = mu(:,(prob.nf+1):(prob.nf+prob.ng));c_pop.sg = sigma(:,(prob.nf+1):(prob.nf+prob.ng));
    end

    %% Predicted archive
    s_archive.x = [s_archive.x;c_pop.x];
    s_archive.muf = [s_archive.muf;c_pop.muf];s_archive.sf = [s_archive.sf;c_pop.sf];
    s_pop.x = [s_pop.x;c_pop.x];
    s_pop.muf = [s_pop.muf;c_pop.muf];s_pop.sf = [s_pop.sf;c_pop.sf];
    if prob.ng > 0
        s_archive.mug = [s_archive.mug;c_pop.mug];s_archive.sg = [s_archive.sg;c_pop.sg];
        s_pop.mug = [s_pop.mug;c_pop.mug];s_pop.sg = [s_pop.sg;c_pop.sg];
    end

    %% Environmental Selection
    count = i;
    [s_pop] = Reduce(param,s_pop,archive,count);
end

%% Final population
% [pop] = Reduce(param,s_archive,archive); % Final pop from whole archive
pop = s_pop; % Final pop from parent + latest child
return



%% IDEA sort with CVsum
function [rank] = IDEASort(Archive, param)
X_pop = Archive.x; F_pop = Archive.muf; G_pop = Archive.mug;
f_new = new_objective(G_pop);
N = size(X_pop,1);
inf_ratio = (param.infeas_ratio*param.EA_popsize);
G1_pop = G_pop;
if(~isempty(G_pop))
    G1_pop(G1_pop <0) = 0;
    cvsum = nansum(G1_pop,2);
    feasible = find(cvsum == 0);
    infeasible = setdiff((1:N)',feasible);
else
    feasible = (1:N)';
    infeasible = [];
    cvsum = zeros(N,1);
end

% sort feasible solutions
ranks1 = [];
if ~isempty(feasible)
    feasiblepop = F_pop(feasible,:);
    if size(F_pop,2) < 2
        [~,id] = sort(feasiblepop); % For single objective
    else
        [ids,~,~,~] = E_NDSort_c(feasiblepop); % For Multi-objective
        id = cell2mat(ids);
    end
    ranks1 = feasible(id);
end

% sort infeasible solutions
ranks2 = [];
if ~isempty(infeasible)
    if length(infeasible) == N
        cv = cvsum(infeasible);
        [~,I] = sort(cv);
        ranks2 = infeasible(I);
    else
        combinedF_pop = [F_pop f_new];
        [I,~,~,~] = E_NDSort_c(combinedF_pop(infeasible));
        InfId = cell2mat(I);
        ranks2 = infeasible(InfId);
    end
end

assert(length(ranks1)+length(ranks2) == N);
if size(ranks1,2) > 1
    ranks1 = ranks1';
elseif size(ranks2,2) > 1
    ranks2 = ranks2';
end

pop_rank = zeros(N,1);
N_inf = min(numel(infeasible), round(inf_ratio));
pop_rank(1:N_inf) = ranks2(1:N_inf);
leftranks2 = setdiff(ranks2,ranks2(1:N_inf),'stable');
pop_rank(N_inf+1:end) = [ranks1;leftranks2];
% promote infeasible solutions
rank = pop_rank(1:N);
return



%% Considering summation of constraint violation as new objective
function [f_new] = new_objective(Gpop)
N = size(Gpop,1);
% constraint violation measure (rank based)
tmp = [Gpop > 0];
c = Gpop .* tmp;
con_f = zeros(N, size(Gpop,2));

for j = 1:size(Gpop,2)
    [tmp, I] = sort(c(:,j));
    cur_rank = 0;
    cur_val = 0;
    for i = 1:N
        if tmp(i) > cur_val
            cur_rank = cur_rank+1;
            cur_val = tmp(i);
        end
        con_f(I(i), j) = cur_rank;
    end
end

f_new = sum(con_f, 2);
return




function [f_new] = new_rank(Gpop)
N = size(Gpop,1);
con_f = ones(N, size(Gpop,2));

for j = 1:size(Gpop,2)
    [tmp, I] = sort(Gpop(:,j));
    cur_rank = 1;
    cur_val = min(Gpop(:,j));
    for i = 1:N
        if tmp(i) > cur_val
            cur_rank = cur_rank+1;
            cur_val = tmp(i);
        end
        con_f(I(i), j) = cur_rank;
    end
end

f_new = sum(con_f, 2);
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
            id = find(d==0);r = nan(1,size(response,2));
            % These Archive solutions have some information
            if ~isempty(id)
                for k=1:length(id)
                    for p=1:size(response,2)
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
function [pop] = Reduce(param,s_archive,archive,gen)
N = param.EA_popsize;
% Initial ranking
[rank,Score,RVs] = RankingDMirror(s_archive,param,archive,gen);
% Environmental selection
if length(rank) >= N
    pop.x = s_archive.x(rank(1:N),:);
    pop.muf = s_archive.muf(rank(1:N),:);pop.sf = s_archive.sf(rank(1:N),:);
    if isfield(archive,'mug') == 1
        pop.mug = s_archive.mug(rank(1:N),:);pop.sg = s_archive.sg(rank(1:N),:);
    end
    pop.Score = Score(1:N);pop.RVs = RVs(1:N);
else
    if gen ~= param.EA_generations
        if mod(length(rank),2) == 1
            id = setdiff((1:size(s_archive.muf,1))',rank,'stable');
            rank = [rank;id(randperm(length(id),1))];
        end
    else
        pop.Score = Score;pop.RVs = RVs;
    end
    pop.x = s_archive.x(rank,:);
    pop.muf = s_archive.muf(rank,:);pop.sf = s_archive.sf(rank,:);
    if isfield(archive,'mug') == 1
        pop.mug = s_archive.mug(rank,:);pop.sg = s_archive.sg(rank,:);
    end
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
function [order,Score,RVs] = RankingDMirror(s_archive,param,archive,gen)
%% Ordering
[order,Score,RVs] = assignmentmirror(s_archive,param,archive,gen);
return


%% Computes the probability of feasibility
function [P, p] = Prob_feas(g_pred,g_sig)
% g_pred <=0 is feasible
p = zeros(1,length(g_pred));
for i = 1:length(g_pred)
    if g_sig(i) ~= 0
        z = (0-g_pred(i))./g_sig(i);
        p(i) = 1-0.5 * erfc(z./sqrt(2));
    else
        if g_pred(i) <= 0
            p(i) = 1;
        else
            p(i) = 0;
        end
    end
end
P = prod(p);
return


%% Solution assignement to reference direction
function [rank,Score,RVs] = assignmentmirror(s_archive,param,archive,gen)
spacing = param.spacing;
muf = s_archive.muf;sigmaf = s_archive.sf;

%% Generating the direction vectors
w = Direction_vector(size(muf,2),spacing);

%% Modification
ideal = archive.ideal;
nadir = archive.nadir;

%% Getting the solutions in scaled space
ff = (muf-repmat(ideal,size(muf,1),1))./repmat(nadir-ideal,size(muf,1),1);
factor = (nadir-ideal).^2;
sigma = sigmaf./repmat(factor,size(muf,1),1);

%% Direction vector
w_p = w - 1;N = size(w,1);
rank = [];reference = (1:size(w,1))';Score = [];RVs = [];
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
        id = find(associate == i); % Solution ids associated with the RVs
        id_tmp = setdiff(id_all(id),rank,'stable');
        if ~isempty(id_tmp)
            %% For Unconstrained: based on PD and for constrained: based on PF and PD
            pop.muf = ff(id_tmp,:);pop.sf = sigma(id_tmp,:);
            pop.mug = s_archive.mug(id_tmp,:);pop.sg = s_archive.sg(id_tmp,:);
            if param.ConstrFlag == 1 % for unconstrained search
                t_w = w(reference(i),:) - w_p(reference(i),:);
                t_f = (pop.muf - w_p(reference(i),:));
                [mu_p,sigma_p] = Multivariate_Projection(t_w,t_f,pop.sf);
                score = CalculateScore(mu_p,sigma_p);
                [FinalScore,id_score] = max(score);
                rank = [rank;id_tmp(id_score)];
                Score = [Score;FinalScore];
            else
                ww = w(reference(i),:);
                muf = pop.muf;sigmaf = pop.sf;mug = pop.mug;sigmag = pop.sg; 
                [FinalScore,order,~] = PFPDCombinedScore(muf,sigmaf,mug,sigmag,ww);
                rank = [rank;id_tmp(order(1))];
                Score = [Score;FinalScore(1)];% FinalScore(1)
            end
            RVs = [RVs;reference(i)];
            N = N - 1;
        else
            tmp = [tmp,i];
        end
    end
    id_all = setdiff(id_all,rank,'stable');
    reference = reference(tmp);
    if gen == param.EA_generations && length(rank) < round(0.5*size(w,1))
        break;
    end
end
return



%% Calculate round-robin for a single column with matrix formulation
function mean_score = CalculateScore(fval,fsigma)
if ~isempty(fval)
    if length(fval) > 1
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
                id3 = find(b_f(id2) < ff(id2));
                s(id2(id3)) = 1; 
                id4 = find(b_f(id2) == ff(id2));
                s(id2(id4)) = 0.5; 
            end
        else
            id2 = find(b_f < ff);
            s(id2) = 1;
            id3 = find(b_f == ff);
            s(id3) = 0.5;
        end
        logmat = ones(size(s)) - eye(size(s));
        s = s.*logmat;
        mean_score = sum(s,2)./(size(s,2)-1);
    else
        mean_score = 0.5;
    end
end
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
function [final_id,Archive] = Probablistic_selection(candidate,Archive,prob,param)
if ~isempty(candidate.x)
    if isfield(Archive,'mug') == 1
        cv = sum(max(Archive.mug,0),2);idfeas = find(cv == 0);
        if isempty(idfeas)
            muf = candidate.muf;sigmaf = candidate.sf;mug = candidate.mug;sigmag = candidate.sg;
            
            if Archive.flag == 0
                if param.ConstrFlag == 2         
                    [~,final_id,~] = PFPDCombinedScore(muf,sigmaf,mug,sigmag);
                else
                    [final_id] = FinalCandidateRanking(candidate);
                end
                Archive.RVs = [Archive.RVs;candidate.RVs(final_id(1))];
                if length(Archive.RVs) == Archive.NumRVs
                    Archive.flag = 0;
                    Archive.RVs = [];
                else
                    Archive.flag = 1;
                end
            else
                if sum(ismember(candidate.RVs,Archive.RVs)) == length(candidate.RVs)
                    Archive.flag = 0;
                    Archive.RVs = [];
                end
                if param.ConstrFlag == 2
                    [~,final_id,~] = PFPDCombinedScore(muf,sigmaf,mug,sigmag);
                else
                    [final_id] = FinalCandidateRanking(candidate);
                end           
                count = 1;curr_RV = [];dis_id = [];
                while isempty(curr_RV) && length(Archive.RVs) < Archive.NumRVs
                    tmp_rv = candidate.RVs(final_id(count));
                    if ismember(tmp_rv,Archive.RVs)
                        dis_id = [dis_id;final_id(count)];
                        if count < length(final_id)
                            count = count + 1;
                            continue;
                        else
                            final_id = final_id(count);
                            Archive.flag = 0;
                            Archive.RVs = [];
                            break;
                        end
                    else
                        final_id = setdiff(final_id,dis_id,'stable');
                        curr_RV = tmp_rv;
                        Archive.RVs = [Archive.RVs;curr_RV];
                    end
                end
                if length(Archive.RVs) == Archive.NumRVs
                    Archive.flag = 0;
                    Archive.RVs = [];
                else
                    Archive.flag = 1;
                end
            end

            %% Visualization
            if param.visualization == 1
                if size(Archive.muf,2) < 3
                    plot(candidate.muf(final_id(1),1),candidate.muf(final_id(1),2),'ko','LineWidth',2);
                else
                    plot3(candidate.muf(final_id(1),1),candidate.muf(final_id(1),2),candidate.muf(final_id(1),3),'ko','LineWidth',2);
                end
            end
        else
            [final_id] = ProbSelectMixFeasInFeas(candidate,Archive,param);

            %% Visualization
            if param.visualization == 1
                if size(Archive.muf,2) < 3
                    plot(candidate.muf(final_id(1),1),candidate.muf(final_id(1),2),'ko','LineWidth',2);
                else
                    plot3(candidate.muf(final_id(1),1),candidate.muf(final_id(1),2),candidate.muf(final_id(1),3),'ko','LineWidth',2);
                end
            end
        end
    else
        [final_id] = ProbSelectAllFeas(candidate,Archive,param);
    end
else
    error('There is no candidate.');
end
return


%% Initial ranking when sufficient ND solutions are not available
function  [rank] = FinalCandidateRanking(pop)
%% Split the feasible and infeasible solutions
if isfield(pop,'mug') == 1
    feas_id = (1:size(pop.mug,1))';
    infeas_id = setdiff((1:size(pop.mug,1))',feas_id,'stable');
else
    feas_id = (1:size(pop.muf,1))';
    infeas_id = [];
end

%% For the feasible solutions
ranks1 = [];
if ~isempty(feas_id)
    muf = pop.muf(feas_id,:);sf = pop.sf(feas_id,:);
    score = RoundRobin(muf,sf);
    [~,order] = sort(score,'descend');
    ranks1 = feas_id(order);
end

%% For the infeasible solutions
ranks2 = [];

%% Assigned percentage of infeasible solutions on top
rank = [ranks1;ranks2];
return





%% Calculate round-robin
function score = RoundRobin(fval,fsigma)
if size(fval,1) > 1
    s = ones(size(fval,1),size(fval,1));
    for i = 1:size(fval,2)
        score_matrix = PDScore(fval(:,i),fsigma(:,i));
        s = s.*score_matrix;
    end
    tmp_score = s';
    sub = s - tmp_score;
    id0 = find(sub == 0);
    s(id0) = 0.5;
    logmat = ones(size(s)) - eye(size(s));
    s = s.*logmat;
    score = sum(s,2)./(size(s,2)-1);
else
    score = 0.5;
end
return


%% Calculate round-robin for a single column with matrix formulation
function s = PDScore(fval,fsigma)
if ~isempty(fval)
    if length(fval) > 1
        ff = repmat(fval',length(fval),1);
        ssigma = repmat(fsigma',length(fval),1);
        s = zeros(size(ff));
        b_f = repmat(fval,1,length(fval));
        b_sigma = repmat(fsigma,1,length(fval));
        denom = (sqrt(2.*(b_sigma.^2+ssigma.^2)));
        id1 = find(denom ~= 0);
        if ~isempty(id1)
            s(id1) = 0.5+0.5.*erf((ff(id1)-b_f(id1))./(sqrt(2.*(b_sigma(id1).^2+ssigma(id1).^2))));
            id2 = setdiff(1:size(s,1)*size(s,2),id1,'stable');
            if ~isempty(id2)
                id3 = find(b_f(id2) < ff(id2));
                s(id2(id3)) = 1; 
                id4 = find(b_f(id2) == ff(id2));
                s(id2(id4)) = 0.5; 
            end
        else
            id2 = find(b_f < ff);
            s(id2) = 1;
            id3 = find(b_f == ff);
            s(id3) = 0.5;
        end
    else
        s = 0.5;
    end
else
    s = [];
end
return


%% Probabilistic selection whole archive is feasible
function [final_id] = ProbSelectAllFeas(candidate,Archive,param)
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
tmp_F = [fnd;Archive.shadow_nd];
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
%% Visualization
if param.visualization == 1
    Ref = tmp_F;
    if size(Archive.muf,2) < 3
        plot(Ref(:,1),Ref(:,2),'bo','LineWidth',2);
        plot(candidate.muf(final_id(1),1),candidate.muf(final_id(1),2),'ko','LineWidth',2);
    else
        plot3(Ref(:,1),Ref(:,2),Ref(:,3),'bo','LineWidth',2);
        plot3(candidate.muf(final_id(1),1),candidate.muf(final_id(1),2),candidate.muf(final_id(1),3),'ko','LineWidth',2);
    end
end
return




%% Probabilistic selection archive includes mix of feasible and infeasible solutions
function [final_id] = ProbSelectMixFeasInFeas(candidate,Archive,param)
mu = candidate.muf;sigma = candidate.sf;
sigma(sigma < 1e-4) = 1e-4;IDall = 1:size(mu,1);

%% Setting the reference points
cv = sum(max(Archive.mug,0),2);idfeas = find(cv == 0);
idinfeas = setdiff(1:size(Archive.x,1),idfeas,'stable');
FeasF = Archive.muf(idfeas,:);
[id_fronts1,~,~,~] = E_NDSort_c(FeasF);
Fnd = FeasF(id_fronts1{1},:);
if ~isempty(idinfeas)
    InFeasF = Archive.muf(idinfeas,:);
    IL = [];% DL = [];
    for i = 1:size(InFeasF,1)
        tmu = [InFeasF(i,:);Fnd];
        [fronts,~,~,~] = E_NDSort_c(tmu);
        if ismember(1,fronts{1}) && length(fronts{1}) > size(Fnd,1)
            IL = [IL;i];
        end
    end
    Ref = [Fnd;InFeasF(IL,:)];
else
    Ref = Fnd;
end

%% Setting the normalization bounds
ideal = Archive.ideal;nadir = Archive.nadir;

%% Normalization
mu_N = (mu-repmat(ideal,size(mu,1),1))./repmat(nadir-ideal,size(mu,1),1);
factor = (nadir-ideal).^2;
sigma_N = sigma./repmat(factor,size(mu,1),1);

%% Non-dominated set of predicted solution
IL = [];
for i = 1:size(mu,1)
    tmu = [mu(i,:);Ref];
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
tmp_F = [Ref;Archive.shadow_nd];
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

return




%% Creating shadow ND list
function [shadow] = CreateShadowList(Archive)
% Classify the feasible solutions
cv = sum(max(Archive.mug,0),2);id = find(cv == 0);
shadow = Archive.shadow_nd;
if ~ismember(length(cv),id)
    shadow = [shadow;Archive.muf(end,:)];
else
    % Check whether the newly evaluated soln is dominated or not
    F = Archive.muf(id,:);
    [id_fronts,~,~,~] = E_NDSort_c(F);
    % Normalization bound identification
    [ideal,nadir] = BoundIdentification(Archive);
    nadir = nadir + ((nadir-ideal).*0.1);
    F_N = (F-ideal)./(nadir-ideal);
    ff = F_N(1:end-1,:); ff_new = F_N(end,:);
    if ~ismember(size(F,1),id_fronts{1})
        shadow = [shadow;Archive.muf(id(end),:)];
        % Check whether the closest archived solution is dominated or not
        dis = sum((repmat(ff_new,size(ff,1),1) - ff).^2,2);
        [~,id_min] = min(dis);
        if ~ismember(id_min,id_fronts{1})
            shadow = [shadow;Archive.muf(id(id_min),:)];
        end
    else
        % Check whether the closest archived solution is dominated or not
        dis = sum((repmat(ff_new,size(ff,1),1) - ff).^2,2);
        [~,id_min] = min(dis);
        if ~ismember(id_min,id_fronts{1})
            shadow = [shadow;Archive.muf(id(id_min),:)];
        end
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


%% Compute the correlation between F and CV
function [FRank,GRank] = FGRank(pop)
[id_fronts,~,~,~] = E_NDSort_c(pop.muf);
FRank = zeros(size(pop.muf,1),1);
for i = 1:length(id_fronts)
    FRank(id_fronts{i}) = i;
end
cv = sum(max(pop.mug,0),2);
[GRank] = new_rank(cv);
return


function flag = SwitchSearchMechanism(Archive,param)
curr_flag = param.ConstrFlag;
cv = sum(max(Archive.mug,0),2);
idfeas = find(cv == 0);
mu_id = size(Archive.muf,1);
if ~isempty(idfeas)
    flag = 2;
else
    [~,idmin] = min(cv);
    if ismember(mu_id,idmin)
        flag = curr_flag;
    else
        if Archive.RankCor >= 0.27 && param.ConstrFlag == 2
            flag = 1; % For unconstrained search flag;
        else
            flag = 2;
        end
    end
end
return
