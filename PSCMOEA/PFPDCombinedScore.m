%% Rank a population based on their PF and PD score
% Three steps included in this ranking:
% -> Order the pop based on their score following the eqaution: P(sol A being feasible and better than sol B) = (PF_A.*(1-PF_B)) + (PF_A.*PF_B.*PD_AB) + ((1-PF_A).*(1-PF_B).*PF_AB);
% -> Now if any block of solutions with 0 score exists that means they are in the infeasible region, so order the block based on their cv and append
% them under the nonzero score block: this is to maintain Feasibility first strategy;
% -> If the whole block returns with a score of 0 and  order the block based on their PD score;
function [FinalScore,rank,score] = PFPDCombinedScore(muf,sigmaf,mug,sigmag,w)
if size(muf,1) > 1
    for i = 1:size(mug,1)
        [PF(i,:),~] = Prob_feas(mug(i,:),sigmag(i,:));
    end
    [~,~,muCV,sigmaCV] = Sum_Rectified_Gaussians(mug,sigmag);
    PFDomScore = PDScore(muCV,sigmaCV);
    if  nargin > 4
        w_m = w - 1;
        t_w = w - w_m;
        t_f = (muf - w_m);
        [mu_p,sigma_p] = Multivariate_Projection(t_w,t_f,sigmaf);
        PDmatrix = PDScore(mu_p,sigma_p);
    else
        PDmatrix = RoundRobin(muf,sigmaf);
    end
    PFMatrix = repmat(PF',size(muf,1),1);
    b_pf = repmat(PF,1,size(muf,1));
    score = (b_pf.*(1-PFMatrix)) + (b_pf.*PFMatrix.*PDmatrix) + ((1-b_pf).*(1-PFMatrix).*PFDomScore); % Compute the score
    logmat = ones(size(score)) - eye(size(score));
    score = score.*logmat;
    Score = sum(score,2)./(size(score,2)-1);
else
    Score = 0.5;score = 0.5;
end
[FinalScore,order] = sort(Score,'descend');
rank = order;
return


%% Computes the probability of feasibility
function [P, p] = Prob_feas(g_pred,g_sig)
% g_pred <= 0 is feasible
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


%%
function [mur,sigmar,muCV,sigmaCV] = Sum_Rectified_Gaussians(mug,sigmag)
% Rectified Gaussian distribution for CV of individual constraints
sigmag(sigmag == 0) = 1e-16; % To avoid NAN value
a = zeros(size(mug));
b = max((mug+6.*sigmag),0); % Allow six sigma limit
c = (a - mug)./sigmag;
d = (b - mug)./sigmag;
mu_t = (1./sqrt(2*pi).*(exp(-c.^2./2)-exp(-d.^2./2))) + (c./2.*(1+erf(c./sqrt(2)))) + ((d./2).*(1-erf(d./sqrt(2))));
var = sqrt((((mu_t.^2+1)./2).*(erf(d./sqrt(2))-erf(c./sqrt(2)))) - (1./sqrt(2*pi).*((d-2.*mu_t).*exp(-d.^2./2)-(c-2.*mu_t).*exp(-c.^2./2)))...
    + (((c-mu_t).^2./2).*(1+erf(c./sqrt(2)))) + (((d-mu_t).^2./2).*(1-erf(d./sqrt(2)))));
mur = mug+sigmag.*mu_t;
sigmar = max((sigmag.*var),0);
sigmar(sigmar <= 1e-16) = 0;
% Converting to an overall muCV and sigmaCV (Sum of Gaussians)
muCV = sum(mur,2);sigmaCV = sqrt(sum(sigmar.^2,2));
return


%% Calculate round-robin
function score = RoundRobin(fval,fsigma)
s = ones(size(fval,1),size(fval,1));
for i = 1:size(fval,2)
    score_matrix = PDScore(fval(:,i),fsigma(:,i));
    s = s.*score_matrix;
end
tmp_score = s';
sub = s - tmp_score;
id0 = find(sub == 0);
s(id0) = 0.5;
score = s;
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