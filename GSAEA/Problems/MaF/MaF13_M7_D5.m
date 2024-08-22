function [f,g] = MaF13_M7_D5(x)
if nargin == 0
    prob.nf = 7;
    prob.ng = 0;
    prob.nx = 5;
    prob.cast = [];
    prob.bounds(:,1) = [zeros(1,2),zeros(1,prob.nx-2)-2]';
    prob.bounds(:,2) = [ones(1,2),zeros(1,prob.nx-2)+2]';
    f = prob;
    g = [];
else
    [f,g] = MaF13_M7_D5_true(x);
end
return

function [PopObj,gg] = MaF13_M7_D5_true(X)
[N,D] = size(X);
M = 7;
Y = X - 2*repmat(X(:,2),1,D).*sin(2*pi*repmat(X(:,1),1,D)+repmat(1:D,N,1)*pi/D);
PopObj(:,1) = sin(X(:,1)*pi/2)                   + 2*mean(Y(:,4:3:D).^2,2);
PopObj(:,2) = cos(X(:,1)*pi/2).*sin(X(:,2)*pi/2) + 2*mean(Y(:,5:3:D).^2,2);
PopObj(:,3) = cos(X(:,1)*pi/2).*cos(X(:,2)*pi/2) + 2*mean(Y(:,3:3:D).^2,2);
PopObj(:,4:M) = repmat(PopObj(:,1).^2+PopObj(:,2).^10+PopObj(:,3).^10+2*mean(Y(:,4:D).^2,2),1,M-3);
gg = [];
return