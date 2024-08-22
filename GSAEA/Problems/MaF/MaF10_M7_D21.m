function [f,g] = MaF10_M7_D21(x)
if nargin == 0
    prob.nf = 7;
    prob.ng = 0;
    prob.nx = 3*prob.nf;
    prob.cast = [];
    prob.bounds(:,1) = zeros(prob.nx,1);
    prob.bounds(:,2) = (2 : 2 : 2*prob.nx)';
    f = prob;
    g = [];
else
    [f,g] = MaF10_2d_true(x);
end
return

function [f,gg] = MaF10_2d_true(PopDec)
[N,D] = size(PopDec);
M = 7;
K = M - 1;
L = D - K;
D = 1;
S = 2 : 2 : 2*M;
A = ones(1,M-1);

z01 = PopDec./repmat(2:2:size(PopDec,2)*2,N,1);

t1 = zeros(N,K+L);
t1(:,1:K)     = z01(:,1:K);
t1(:,K+1:end) = s_linear(z01(:,K+1:end),0.35);

t2 = zeros(N,K+L);
t2(:,1:K)     = t1(:,1:K);
t2(:,K+1:end) = b_flat(t1(:,K+1:end),0.8,0.75,0.85);

t3 = b_poly(t2,0.02);

t4 = zeros(N,M);
for i = 1 : M-1
    t4(:,i) = r_sum(t3(:,(i-1)*K/(M-1)+1:i*K/(M-1)),2*((i-1)*K/(M-1)+1):2:2*i*K/(M-1));
end
t4(:,M) = r_sum(t3(:,K+1:K+L),2*(K+1):2:2*(K+L));

x = zeros(N,M);
for i = 1 : M-1
    x(:,i) = max(t4(:,M),A(i)).*(t4(:,i)-0.5)+0.5;
end
x(:,M) = t4(:,M);

h      = convex(x);
h(:,M) = mixed(x);
f = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;
gg = [];
return


function Output = s_linear(y,A)
Output = abs(y-A)./abs(floor(A-y)+A);
return

function Output = b_flat(y,A,B,C)
Output = A+min(0,floor(y-B))*A.*(B-y)/B-min(0,floor(C-y))*(1-A).*(y-C)/(1-C);
Output = roundn(Output,-6);
return

function Output = b_poly(y,a)
Output = y.^a;
return

function Output = r_sum(y,w)
Output = sum(y.*repmat(w,size(y,1),1),2)./sum(w);
return

function Output = convex(x)
Output = fliplr(cumprod([ones(size(x,1),1),1-cos(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),1-sin(x(:,end-1:-1:1)*pi/2)];
return

function Output = mixed(x)
Output = 1-x(:,1)-cos(10*pi*x(:,1)+pi/2)/10/pi;
return