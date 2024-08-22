%% True pareto equations of the given problems
function P = PF(prob_name,N,M)
switch prob_name
    case 'MaF1'
        P = 1 - UniformPoint(N,M);
    case 'MaF2'
        P = UniformPoint(N,M);
        c = zeros(size(P,1),M-1);
        for i = 1 : size(P,1)
            for j = 2 : M
                temp = P(i,j)/P(i,1)*prod(c(i,M-j+2:M-1));
                c(i,M-j+1) = sqrt(1/(1+temp^2));
            end
        end
        if M > 5
            c = c.*(cos(pi/8)-cos(3*pi/8)) + cos(3*pi/8);
        else
            c(any(c<cos(3*pi/8)|c>cos(pi/8),2),:) = [];
        end
        P = fliplr(cumprod([ones(size(c,1),1),c(:,1:M-1)],2)).*[ones(size(c,1),1),sqrt(1-c(:,M-1:-1:1).^2)];
    case 'MaF3'
        P    = UniformPoint(N,M).^2;
        temp = sum(sqrt(P(:,1:end-1)),2) + P(:,end);
        P    = P./[repmat(temp.^2,1,size(P,2)-1),temp];
    case 'MaF4'
        P = UniformPoint(N,M);
        P = P./repmat(sqrt(sum(P.^2,2)),1,M);
        P = (1-P).*repmat(2.^(1:M),size(P,1),1);
    case 'MaF5'
        P = UniformPoint(N,M);
        P = P./repmat(sqrt(sum(P.^2,2)),1,M);
        P = P.*repmat(2.^(M:-1:1),size(P,1),1);
    case 'MaF6'
        I = 2;
        P = UniformPoint(N,I);
        P = P./repmat(sqrt(sum(P.^2,2)),1,size(P,2));
        P = [P(:,ones(1,M-size(P,2))),P];
        P = P./sqrt(2).^repmat(max([M-I,M-I:-1:2-I],0),size(P,1),1);
    case 'MaF7'
        interval     = [0,0.251412,0.631627,0.859401];
        median       = (interval(2)-interval(1))/(interval(4)-interval(3)+interval(2)-interval(1));
        X            = ReplicatePoint(N,M-1);
        X(X<=median) = X(X<=median)*(interval(2)-interval(1))/median+interval(1);
        X(X>median)  = (X(X>median)-median)*(interval(4)-interval(3))/(1-median)+interval(3);
        P            = [X,2*(M-sum(X/2.*(1+sin(3*pi.*X)),2))];
    case 'MaF8'
        % Generate vertexes
        Points  = [];
        [thera,rho] = cart2pol(0,1);
        [Points(:,1),Points(:,2)] = pol2cart(thera-(1:M)*2*pi/M,rho);
        [X,Y] = ndgrid(linspace(-1,1,ceil(sqrt(N))));
        ND    = inpolygon(X(:),Y(:),Points(:,1),Points(:,2));
        P     = pdist2([X(ND),Y(ND)],Points);
    case 'MaF9'
        if M == 3
            load('MaF9_M3_PF.mat');
        elseif M == 5
            load('MaF9_M5_PF.mat');
        elseif M == 7
            load('MaF9_M7_PF.mat');
        end
    case 'MaF10'
        P = UniformPoint(N,M);
        c = ones(size(P,1),M);
        for i = 1 : size(P,1)
            for j = 2 : M
                temp = P(i,j)/P(i,1)*prod(1-c(i,M-j+2:M-1));
                c(i,M-j+1) = (temp^2-temp+sqrt(2*temp))/(temp^2+1);
            end
        end
        x = acos(c)*2/pi;
        temp = (1-sin(pi/2*x(:,2))).*P(:,M)./P(:,M-1);
        a = 0 : 0.0001 : 1;
        E = abs(temp*(1-cos(pi/2*a))-1+repmat(a+cos(10*pi*a+pi/2)/10/pi,size(x,1),1));
        [~,rank] = sort(E,2);
        for i = 1 : size(x,1)
            x(i,1) = a(min(rank(i,1:10)));
        end
        P      = convex(x);
        P(:,M) = mixed(x);
        P      = repmat(2:2:2*M,size(P,1),1).*P;
    case 'MaF11'
        P = UniformPoint(N,M);
        c = ones(size(P,1),M);
        for i = 1 : size(P,1)
            for j = 2 : M
                temp = P(i,j)/P(i,1)*prod(1-c(i,M-j+2:M-1));
                c(i,M-j+1) = (temp^2-temp+sqrt(2*temp))/(temp^2+1);
            end
        end
        x = acos(c)*2/pi;
        temp = (1-sin(pi/2*x(:,2))).*P(:,M)./P(:,M-1);
        a = 0 : 0.0001 : 1;
        E = abs(temp*(1-cos(pi/2*a))-1+repmat(a.*cos(5*pi*a).^2,size(x,1),1));
        [~,rank] = sort(E,2);
        for i = 1 : size(x,1)
            x(i,1) = a(min(rank(i,1:10)));
        end
        P      = convex(x);
        P(:,M) = disc(x);
        P      = P(NDSort(P,1)==1,:);
        P      = repmat(2:2:2*M,size(P,1),1).*P;
    case 'MaF12'
        P = UniformPoint(N,M);
        P = P./repmat(sqrt(sum(P.^2,2)),1,M);
        P = repmat(2:2:2*M,size(P,1),1).*P;
    case 'MaF13'
        P = UniformPoint(N,3);
        P = P./repmat(sqrt(sum(P.^2,2)),1,3);
        P = [P,repmat(P(:,1).^2+P(:,2).^10+P(:,3).^10,1,M-3)];
    case 'MaF14'
        P = UniformPoint(N,M);
    case 'MaF15'
        P = UniformPoint(N,M);
        P = 1 - P./repmat(sqrt(sum(P.^2,2)),1,M);
    case 'TR1'
        x = linspace(0,1,N)';
        f = 0.5*(1-x.^0.3).^(1/0.3);
        ce = [x f];
        pall = N-1;
        [r_ori] = Direction_vector(M,pall);
        lr = -1*r_ori;rr = flipud(r_ori);
        rr = rr+5*ones(size(rr,1),2);
        r_ori = flipud(r_ori);p2 = [];
        for i = 1:size(lr,1)
            rr(i,1) = 1000;
            rr(i,2) = 1000*(r_ori(i,2)-lr(i,2))/(r_ori(i,1)-lr(i,1));
            l2 = [lr(i,1) rr(i,1);lr(i,2) rr(i,2)];
            P = InterX([ce(:,1)';ce(:,2)'],l2);
            if (~isempty(P))
                p2(i,1) = P(1);
                p2(i,2) = P(2);
            end
        end
        P = p2;

    case 'TR2'
        x = linspace(0,1,N)';
        f = 0.5*(1-x.^6).^(1/6);
        ce = [x f];
        pall = N-1;
        [r_ori] = Direction_vector(M,pall);
        lr = -1*r_ori;rr = flipud(r_ori);
        rr = rr+5*ones(size(rr,1),2);
        r_ori = flipud(r_ori);p2 = [];
        for i = 1:size(lr,1)
            rr(i,1) = 1000;
            rr(i,2) = 1000*(r_ori(i,2)-lr(i,2))/(r_ori(i,1)-lr(i,1));
            l2 = [lr(i,1) rr(i,1);lr(i,2) rr(i,2)];
            P = InterX([ce(:,1)';ce(:,2)'],l2);
            if (~isempty(P))
                p2(i,1) = P(1);
                p2(i,2) = P(2);
            end
        end
        P = p2;


    case 'TR3'
        x = linspace(-5,5,N)';
        f = tanh(-x);
        ce = [x f];LB = [-5 -1];UB = [5 1];
        n_ce = (ce - repmat(LB,size(ce,1),1))./(repmat(UB,size(ce,1),1) - repmat(LB,size(ce,1),1));
        pall = N-1;
        [r_ori] = Direction_vector(M,pall);
        lr = -1*r_ori;rr = flipud(r_ori);
        rr = rr+5*ones(size(rr,1),2);
        r_ori = flipud(r_ori);p2 = [];
        for i = 1:size(lr,1)
            rr(i,1) = 1000;
            rr(i,2) = 1000*(r_ori(i,2)-lr(i,2))/(r_ori(i,1)-lr(i,1));
            l2 = [lr(i,1) rr(i,1);lr(i,2) rr(i,2)];
            P = InterX([n_ce(:,1)';n_ce(:,2)'],l2);
            if (~isempty(P))
                p2(i,1) = P(1);
                p2(i,2) = P(2);
            end
        end
        p2 = repmat(LB,size(p2,1),1) + p2.*(repmat(UB,size(p2,1),1) - repmat(LB,size(p2,1),1));
        P = p2;

    case 'TR4'
        if N == 1000
            step = 51;
        elseif N == 100
            step = 11;
        else
            step = 101;
        end
        xx = [linspace(0,1,step)',linspace(0,1,step)'];
        [xmesh,ymesh] = meshgrid(xx(:,1),xx(:,2));
        x = [reshape(xmesh,size(xmesh,1)*size(xmesh,2),1) reshape(ymesh,size(ymesh,1)*size(ymesh,2),1)];
        P(:,1) = (cos(0.5*pi.*(x(:,1))).*cos(0.5*pi.*(x(:,2)))).^4;
        P(:,2) = (cos(0.5*pi.*(x(:,1))).*sin(0.5*pi.*(x(:,2)))).^4;
        P(:,M) = (sin(0.5*pi.*(x(:,1)))).^2;

    case 'MOPCH'
        P = UniformPoint(N,M);
        P = 1-(P.^0.35);

    case 'DTLZ1'
        P = UniformPoint(N,M)/2;

    case 'DTLZ2'
        P = UniformPoint(N,M);
        P = P./repmat(sqrt(sum(P.^2,2)),1,M);

    case 'DTLZ3'
        P = UniformPoint(N,M);
        P = P./repmat(sqrt(sum(P.^2,2)),1,M);

    case 'DTLZ4'
        P = UniformPoint(N,M);
        P = P./repmat(sqrt(sum(P.^2,2)),1,M);

    case 'DTLZ5'
        P = [0:1/(N-1):1;1:-1/(N-1):0]';
        P = P./repmat(sqrt(sum(P.^2,2)),1,size(P,2));
        P = [P(:,ones(1,M-2)),P];
        P = P./sqrt(2).^repmat([M-2,M-2:-1:0],size(P,1),1);

    case 'DTLZ6'
        P = [0:1/(N-1):1;1:-1/(N-1):0]';
        P = P./repmat(sqrt(sum(P.^2,2)),1,size(P,2));
        P = [P(:,ones(1,M-2)),P];
        P = P./sqrt(2).^repmat([M-2,M-2:-1:0],size(P,1),1);

    case 'ZDT1'
        P(:,1) = (0:1/(N-1):1)';
        P(:,2) = 1 - P(:,1).^0.5;

    case 'ZDT2'
        P(:,1) = (0:1/(N-1):1)';
        P(:,2) = 1 - P(:,1).^2;

    case 'ZDT3'
        P(:,1) = (0:1/(N-1):1)';
        P(:,2) = 1 - P(:,1).^0.5 - P(:,1).*sin(10*pi*P(:,1));
        P      = P(NDSort(P,1)==1,:);

    case 'ZDT4'
        P(:,1) = (0:1/(N-1):1)';
        P(:,2) = 1 - P(:,1).^0.5;

    case 'ZDT6'
        minf1  = 0.280775;
        P(:,1) = (minf1:(1-minf1)/(N-1):1)';
        P(:,2) = 1 - P(:,1).^2;

    case 'MW1'
        R(:,1) = (0:1/(N-1):1)';
        R(:,2) = 1 - 0.85*R(:,1);
        l = sqrt(2)*R(:,2) - sqrt(2)*R(:,1);
        c = 1 - R(:,1) - R(:,2) + 0.5*sin(2*pi*l).^8;
        R(c<0,:) = [];
        P = R;

    case 'MW2'
        R(:,1) = (0:1/(N-1):1)';
        R(:,2) = 1 - R(:,1);
        P = R;

    case 'MW3'
        R(:,1)  = (0:1/(N-1):1)';
        R(:,2)  = 1 - R(:,1);
        invalid = (0.85-R(:,1)-R(:,2)+0.3*sin(0.75*pi*sqrt(2)*(R(:,2)-R(:,1))).^2) > 0;
        while any(invalid)
            R(invalid,:) = R(invalid,:).*1.001;
            invalid = (0.85-R(:,1)-R(:,2)+0.3*sin(0.75*pi*sqrt(2)*(R(:,2)-R(:,1))).^2) > 0;
        end
        P = R;

    case 'MW4'
        R = UniformPoint(N,M);
        l = R(:,end) - sum(R(:,1:end-1),2);
        c = (1+0.4*sin(2.5*pi*l).^8) - sum(R,2);
        R(c<0,:) = [];
        P = R;

    case 'MW5'
        R = [0 1;0.3922 0.9199;0.4862 0.8739;0.5490 0.8358;0.5970 0.8023;0.6359 0.7719;0.6686 0.7436;0.6969 0.7174];
        R = [R;flip(R,2)];
        P = R;

    case 'MW6'
        R(:,1) = (0:1/(N-1):1)';
        R(:,2) = 1 - R(:,1);
        R = R./repmat(sqrt(sum(R.^2,2)/1.21),1,2);
        l = cos(6*atan(R(:,2)./R(:,1)).^4).^10;
        c = 1 - (R(:,1)./(1+0.15*l)).^2 - (R(:,2)./(1+0.75*l)).^2;
        R(c<0,:) = [];
        P = R;

    case 'MW7'
        R(:,1) = (0:1/(N-1):1)';
        R(:,2) = 1 - R(:,1);
        R = R./repmat(sqrt(sum(R.^2,2)),1,2);
        invalid = ((1.15-0.2*sin(4*atan(R(:,2)./R(:,1))).^8).^2-R(:,1).^2-R(:,2).^2) > 0;
        while any(invalid)
            R(invalid,:) = R(invalid,:).*1.001;
            invalid = ((1.15-0.2*sin(4*atan(R(:,2)./R(:,1))).^8).^2-R(:,1).^2-R(:,2).^2) > 0;
        end
        R = R(NDSort(R,1)==1,:);
        P = R;

    case 'MW8'
        R = UniformPoint(N,M);
        R = R./repmat(sqrt(sum(R.^2,2)),1,M);
        R(1-(1.25 - 0.5*sin(6*asin(R(:,end))).^2).^2>0,:) = [];
        P = R;

    case 'MW9'
        R(:,1)  = (0:1/(N-1):1)';
        R(:,2)  = 1 - R(:,1).^0.6;
        T1      = (1-0.64*R(:,1).^2-R(:,2)).*(1-0.36*R(:,1).^2-R(:,2));
        T2      = 1.35.^2 - (R(:,1)+0.35).^2 - R(:,2);
        T3      = 1.15.^2 - (R(:,1)+0.15).^2 - R(:,2);
        invalid = min(T1,T2.*T3) > 0;
        while any(invalid)
            R(invalid,:) = R(invalid,:).*1.001;
            T1      = (1-0.64*R(:,1).^2-R(:,2)).*(1-0.36*R(:,1).^2-R(:,2));
            T2      = 1.35.^2 - (R(:,1)+0.35).^2 - R(:,2);
            T3      = 1.15.^2 - (R(:,1)+0.15).^2 - R(:,2);
            invalid = min(T1,T2.*T3) > 0;
        end
        R = R(NDSort(R,1)==1,:);
        P = R;

    case 'MW10'
        R(:,1)  = (0:1/(N-1):1)';
        R(:,2)  = 1 - R(:,1).^2;
        c1      = (2 - 4*R(:,1).^2 - R(:,2)).*(2 - 8*R(:,1).^2 - R(:,2));
        c2      = (2 - 2*R(:,1).^2 - R(:,2)).*(2 - 16*R(:,1).^2 - R(:,2));
        c3      = (1 - R(:,1).^2 - R(:,2)).*(1.2 - 1.2*R(:,1).^2 - R(:,2));
        invalid = c1<0 | c2>0 | c3>0;
        while any(invalid)
            R(invalid,:) = R(invalid,:).*1.001;
            R(any(R>1.3,2),:) = [];
            c1      = (2 - 4*R(:,1).^2 - R(:,2)).*(2 - 8*R(:,1).^2 - R(:,2));
            c2      = (2 - 2*R(:,1).^2 - R(:,2)).*(2 - 16*R(:,1).^2 - R(:,2));
            c3      = (1 - R(:,1).^2 - R(:,2)).*(1.2 - 1.2*R(:,1).^2 - R(:,2));
            invalid = c1<0 | c2>0 | c3>0;
        end
        R = R(NDSort(R,1)==1,:);
        P = R;

    case 'MW11'
        R(:,1)  = (0:1/(N-1):1)';
        R(:,2)  = 1 - R(:,1);
        R       = R./repmat(sqrt(sum(R.^2,2)/2),1,2);
        c1      = (3 - R(:,1).^2 - R(:,2)).*(3 - 2*R(:,1).^2 - R(:,2));
        c2      = (3 - 0.625*R(:,1).^2 - R(:,2)).*(3 - 7*R(:,1).^2 - R(:,2));
        c3      = (1.62 - 0.18*R(:,1).^2 - R(:,2)).*(1.125 - 0.125*R(:,1).^2 - R(:,2));
        c4      = (2.07 - 0.23*R(:,1).^2 - R(:,2)).*(0.63 - 0.07*R(:,1).^2 - R(:,2));
        invalid = c1<0 | c2>0 | c3<0 | c4>0;
        while any(invalid)
            R(invalid,:) = R(invalid,:).*1.001;
            R(any(R>2.2,2),:) = [];
            c1      = (3 - R(:,1).^2 - R(:,2)).*(3 - 2*R(:,1).^2 - R(:,2));
            c2      = (3 - 0.625*R(:,1).^2 - R(:,2)).*(3 - 7*R(:,1).^2 - R(:,2));
            c3      = (1.62 - 0.18*R(:,1).^2 - R(:,2)).*(1.125 - 0.125*R(:,1).^2 - R(:,2));
            c4      = (2.07 - 0.23*R(:,1).^2 - R(:,2)).*(0.63 - 0.07*R(:,1).^2 - R(:,2));
            invalid = c1<0 | c2>0 | c3<0 | c4>0;
        end
        R = [R;1,1];
        R = R(NDSort(R,1)==1,:);
        P = R;

    case 'MW12'
        R(:,1)  = (0:1/(N-1):1)';
        R(:,2)  = 0.85 - 0.8*R(:,1) - 0.08*abs(sin(3.2*pi*R(:,1)));
        c1      = (1-0.8*R(:,1)-R(:,2)+0.08*sin(2*pi*(R(:,2)-R(:,1)/1.5))).*(1.8-1.125*R(:,1)-R(:,2)+0.08*sin(2*pi*(R(:,2)/1.8-R(:,1)/1.6)));
        invalid = c1>0;
        while any(invalid)
            R(invalid,:) = R(invalid,:).*1.001;
            c1      = (1-0.8*R(:,1)-R(:,2)+0.08*sin(2*pi*(R(:,2)-R(:,1)/1.5))).*(1.8-1.125*R(:,1)-R(:,2)+0.08*sin(2*pi*(R(:,2)/1.8-R(:,1)/1.6)));
            invalid = c1>0;
        end
        P = R;

    case 'MW13'
        R(:,1)  = (0:1.5/(N-1):1.5)';
        R(:,2)  = 5 - exp(R(:,1)) - 0.5*abs(sin(3*pi*R(:,1)));
        c1      = (5-exp(R(:,1))-0.5*sin(3*pi*R(:,1))-R(:,2)).*(5-(1+0.4*R(:,1))-0.5*sin(3*pi*R(:,1))-R(:,2));
        invalid = c1>0;
        while any(invalid)
            R(invalid,:) = R(invalid,:).*1.001;
            c1      = (5-exp(R(:,1))-0.5*sin(3*pi*R(:,1))-R(:,2)).*(5-(1+0.4*R(:,1))-0.5*sin(3*pi*R(:,1))-R(:,2));
            invalid = c1>0;
        end
        R = R(NDSort(R,1)==1,:);
        P = R;

    case 'MW14'
        interval     = [0,0.731000,1.331000,1.500000];
        median       = (interval(2)-interval(1))/(interval(4)-interval(3)+interval(2)-interval(1));
        X            = UniformPoint(N,M-1,'grid');
        X(X<=median) = X(X<=median)*(interval(2)-interval(1))/median+interval(1);
        X(X>median)  = (X(X>median)-median)*(interval(4)-interval(3))/(1-median)+interval(3);
        R            = [X,1/(M-1)*sum(6 - exp(X) - 1.5*sin(1.1*pi*X.^2),2)];
        P = R;

    case 'CF1'
        R(:,1) = (0:1/20:1)';
        R(:,2) = 1 - R(:,1);
        P = R;

    case 'CF2'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - sqrt(R(:,1));
        R(0<R(:,1) & R(:,1)<1/16 | 1/4<R(:,1) & R(:,1)<9/16,:) = [];
        P = R;

    case 'CF3'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - R(:,1).^2;
        R(0<R(:,1) & R(:,1)<1/2 | sqrt(1/2)<R(:,1) & R(:,1)<sqrt(3/4),:) = [];
        P = R;

    case 'CF4'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - R(:,1);
        temp1  = 0.5<R(:,1) & R(:,1)<=0.75;
        temp2  = 0.75<R(:,1);
        R(temp1,2) = -0.5*R(temp1,1) + 3/4;
        R(temp2,2) = 1 - R(temp2,1) + 0.125;
        P = R;

    case 'CF5'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - R(:,1);
        temp1  = 0.5<R(:,1) & R(:,1)<=0.75;
        temp2  = 0.75<R(:,1);
        R(temp1,2) = -0.5*R(temp1,1) + 3/4;
        R(temp2,2) = 1 - R(temp2,1) + 0.125;
        P = R;

    case 'CF6'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = (1-R(:,1)).^2;
        temp1  = 0.5<R(:,1) & R(:,1)<=0.75;
        temp2  = 0.75<R(:,1);
        R(temp1,2) = 0.5*(1-R(temp1,1));
        R(temp2,2) = 0.25*sqrt(1-R(temp2,1));
        P = R;

    case 'CF7'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = (1-R(:,1)).^2;
        temp1  = 0.5<R(:,1) & R(:,1)<=0.75;
        temp2  = 0.75<R(:,1);
        R(temp1,2) = 0.5*(1-R(temp1,1));
        R(temp2,2) = 0.25*sqrt(1-R(temp2,1));
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
        R(1e-5<R(:,1) & R(:,1)<sqrt((1-R(:,3).^2)/4) | sqrt((1-R(:,3).^2)/2)<R(:,1) & R(:,1)<sqrt(3*(1-R(:,3).^2)/4),:) = [];
        P = R;

    case 'CF10'
        R = UniformPoint(N,3);
        R = R./repmat(sqrt(sum(R.^2,2)),1,3);
        R(1e-5<R(:,1) & R(:,1)<sqrt((1-R(:,3).^2)/4) | sqrt((1-R(:,3).^2)/2)<R(:,1) & R(:,1)<sqrt(3*(1-R(:,3).^2)/4),:) = [];
        P = R;

    case 'DASCMOP1'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - R(:,1).^2;
        R      = R + 0.5;
        X1     = (sqrt(1-4*(-R(:,1)+R(:,2)-1))-1)/2;
        sum1   = R(:,1) - X1;
        C      = Constraint1(X1,sum1);
        R(any(C>0,2),:) = [];
        R = [R;1.5,0.5];
        P = R;

    case 'DASCMOP2'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - sqrt(R(:,1));
        R      = R + 0.5;
        theta_k = -0.25 * pi;
        C = 0.25 - ((R(:,1) - 1) * cos(theta_k) - (R(:,2) - 0.5) * sin(theta_k)).^2 ./ 0.3 -...
            ((R(:,1) - 1) * sin(theta_k) + (R(:,2) - 0.5) * cos(theta_k)).^2 ./ 1.2;
        invalid = C>0;
        while any(invalid)
            R(invalid,:) = (R(invalid,:)-0.5).*1.001 + 0.5;
            C = 0.25 - ((R(:,1) - 1) * cos(theta_k) - (R(:,2) - 0.5) * sin(theta_k)).^2 ./ 0.3 -...
                ((R(:,1) - 1) * sin(theta_k) + (R(:,2) - 0.5) * cos(theta_k)).^2 ./ 1.2;
            invalid = C>0;
        end
        P = R;

    case 'DASCMOP3'
        R = [0.5000,1.5000;0.5010,1.4762;0.5020,1.4710;0.5030,1.4688;0.5040,1.4681;0.6502,1.4652;0.7002,1.0541;
            0.9044,0.8986;1.1066,0.7729;1.3008,0.6114;1.5000,0.5000;0.9069,0.8951;1.1126,0.7727;0.9129,0.8950;
            1.1151,0.7690;0.9153,0.8914;1.1175,0.7653;1.1200,0.7616;0.9213,0.8913;1.1260,0.7613;1.1285,0.7576];
        P = R;

    case 'DASCMOP4'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - R(:,1).^2;
        R      = R + 0.5;
        X1     = (sqrt(1-4*(-R(:,1)+R(:,2)-1))-1)/2;
        sum1   = R(:,1) - X1;
        C      = Constraint2(X1,sum1);
        R(any(C>0,2),:) = [];
        R = [R;1.5,0.5];
        P = R;

    case 'DASCMOP5'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - sqrt(R(:,1));
        R      = R + 0.5;
        R(sin(20*pi*R(:,1))<-1e-10,:) = [];
        theta_k = -0.25 * pi;
        C = 0.25 - ((R(:,1) - 1) * cos(theta_k) - (R(:,2) - 0.5) * sin(theta_k)).^2 ./ 0.3 -...
            ((R(:,1) - 1) * sin(theta_k) + (R(:,2) - 0.5) * cos(theta_k)).^2 ./ 1.2;
        invalid = C>0;
        while any(invalid)
            R(invalid,:) = (R(invalid,:)-0.5).*1.001 + 0.5;
            C = 0.25 - ((R(:,1) - 1) * cos(theta_k) - (R(:,2) - 0.5) * sin(theta_k)).^2 ./ 0.3 -...
                ((R(:,1) - 1) * sin(theta_k) + (R(:,2) - 0.5) * cos(theta_k)).^2 ./ 1.2;
            invalid = C>0;
        end
        P = R;

    case 'DASCMOP6'
        R = [0.5000,1.5000;0.5010,1.4762;0.5020,1.4710;0.5030,1.4688;0.5040,1.4681;0.6502,1.4652;0.7002,1.0541;
            0.9044,0.8986;1.1066,0.7729;1.3008,0.6114;1.5000,0.5000;0.9069,0.8951;1.1126,0.7727;0.9129,0.8950;
            1.1151,0.7690;0.9153,0.8914;1.1175,0.7653;1.1200,0.7616;0.9213,0.8913;1.1260,0.7613;1.1285,0.7576];
        P = R;

    case 'DASCMOP7'
        R = UniformPoint(N,3);
        X(:,1) = 1./(1+R(:,2)./R(:,1));
        X(:,2) = R(:,1)./X(:,1);
        C(:,1) = -sin(20*pi*X(:,1));
        C(:,2) = -cos(20*pi*X(:,2));
        R(any(C>1e-2,2),:) = [];
        R = R + 0.5;
        P = R;

    case 'DASCMOP8'
        R = UniformPoint(N,3);
        R = R./repmat(sqrt(sum(R.^2,2)),1,3);
        X(:,2) = atan(R(:,2)./R(:,1))/0.5/pi;
        X(:,1) = acos(R(:,1)./cos(0.5*pi*X(:,2)))/0.5/pi;
        C(:,1) = -sin(20*pi*X(:,1));
        C(:,2) = -cos(20*pi*X(:,2));
        R(any(C>1e-2,2),:) = [];
        R = R + 0.5;
        P = R;

    case 'DASCMOP9'
        R = UniformPoint(N,3);
        R = R./repmat(sqrt(sum(R.^2,2)),1,3);
        X(:,2) = atan(R(:,2)./R(:,1))/0.5/pi;
        X(:,1) = acos(R(:,1)./cos(0.5*pi*X(:,2)))/0.5/pi;
        C(:,1) = -sin(20*pi*X(:,1));
        C(:,2) = -cos(20*pi*X(:,2));
        R(any(C>1e-2,2),:) = [];
        R = R + 0.5;
        P = R;

    case 'LIRCMOP1'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - R(:,1).^2;
        R      = R + 0.5;
        P = R;

    case 'LIRCMOP2'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - sqrt(R(:,1));
        R      = R + 0.5;
        P = R;

    case 'LIRCMOP3'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - R(:,1).^2;
        R(sin(20*pi*R(:,1))<0.5,:) = [];
        R      = R + 0.5;
        P = R;

    case 'LIRCMOP4'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - sqrt(R(:,1));
        R(sin(20*pi*R(:,1))<0.5,:) = [];
        R      = R + 0.5;
        P = R;

    case 'LIRCMOP5'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - sqrt(R(:,1));
        R      = R + 0.7057;
        R(any(Constraint4(R)>0,2),:) = [];
        P = R;

    case 'LIRCMOP6'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - R(:,1).^2;
        R      = R + 0.7057;
        R(any(Constraint5(R)>0,2),:) = [];
        P = R;

    case 'LIRCMOP7'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - sqrt(R(:,1));
        R      = R + 0.7057;
        theta  = -0.25*pi;
        c1     = 0.1 - ((R(:,1)-1.2)*cos(theta)-(R(:,2)-1.2)*sin(theta)).^2/(2^2) -...
            ((R(:,1)-1.2)*sin(theta)+(R(:,2)-1.2)*cos(theta)).^2/(6^2);
        invalid = c1>0;
        while any(invalid)
            R(invalid,:) = (R(invalid,:)-0.7057).*1.001 + 0.7057;
            c1 = 0.1 - ((R(:,1)-1.2)*cos(theta)-(R(:,2)-1.2)*sin(theta)).^2/(2^2) -...
                ((R(:,1)-1.2)*sin(theta)+(R(:,2)-1.2)*cos(theta)).^2/(6^2);
            invalid = c1>0;
        end
        P = R;

    case 'LIRCMOP8'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - sqrt(R(:,1));
        R      = R + 0.7057;
        theta  = -0.25*pi;
        c1     = 0.1 - ((R(:,1)-1.2)*cos(theta)-(R(:,2)-1.2)*sin(theta)).^2/(2^2) -...
            ((R(:,1)-1.2)*sin(theta)+(R(:,2)-1.2)*cos(theta)).^2/(6^2);
        invalid = c1>0;
        while any(invalid)
            R(invalid,:) = (R(invalid,:)-0.7057).*1.001 + 0.7057;
            c1 = 0.1 - ((R(:,1)-1.2)*cos(theta)-(R(:,2)-1.2)*sin(theta)).^2/(2^2) -...
                ((R(:,1)-1.2)*sin(theta)+(R(:,2)-1.2)*cos(theta)).^2/(6^2);
            invalid = c1>0;
        end
        P = R;

    case 'LIRCMOP9'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - R(:,1).^2;
        R      = R*1.7057;
        R(any(Constraint6(R)>0,2),:) = [];
        R = [R;0,2.182;1.856,0];
        P = R;

    case 'LIRCMOP10'
        R(:,1) = linspace(0,1,N)';
        R(:,2) = 1 - sqrt(R(:,1));
        R      = R*1.7057;
        R(any(Constraint7(R)>0,2),:) = [];
        R = [R;1.747,0];
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
        R = 1.7057*R./repmat(sqrt(sum(R.^2,2)),1,3);
        P = R;

    case 'LIRCMOP14'
        R = UniformPoint(N,3);
        R = sqrt(3.0625)*R./repmat(sqrt(sum(R.^2,2)),1,3);
        P = R;

    case 'FCP1'
        P = 8.5*UniformPoint(N,M);
    
    case 'FCP1Mod'
        t = 0.5*pi*(0:1/N:1)';
        P = 8.5*[cos(t),sin(t)];

    case 'Test1'
        P = 8.5*UniformPoint(N,M);

    case 'Test2'
        t = 0.5*pi*(0:1/N:1)';
        P = 8.5*[cos(t),sin(t)];

    case 'FCP2'
        t = 0:1/N:1;
        f1 = cos(0.5*pi*t);
        f2 = sin(0.5*pi*t)+0.2*sin(4*pi*t);
        P = [f1',f2'];
        P = P(find(NDSort(P,1)==1),:);
        P = 8.5*P;

    case 'FCP3'
        t = 0.5*pi*(0:1/N:1)';
        P = 8.5*[cos(t),sin(t)];

    case 'FCP4'
        t = 0:1/N:1;
        f1 = 1-t;
        f2 = t+0.2*sin(4*pi*t);
        P = [f1',f2'];
        P = P(find(NDSort(P,1)==1),:);
        P = 8.5*P;

    case 'FCP5'
        t=0.5*pi*(0:1/N:1);
        c1x1 = 0.1*[9+0.5*cos(t),9-0.5*cos(t)];
        c1g  = repmat(3-0.5*sin(t),1,2);
        c2x1 = 0.1*[6+0.95*cos(t),6-0.95*cos(t)];
        c2g  = repmat(6-0.95*sin(t),1,2);
        c3x1 = 0.1*[sqrt(2)+sqrt(2)*cos(t),sqrt(2)-sqrt(2)*cos(t)];
        c3g  = repmat(10-sqrt(2)*sin(t),1,2);

        x1 = [c1x1,c2x1,c3x1];
        g = [c1g,c2g,c3g];
        P1 = x1.*g;
        P2 = (1-x1).*g;
        P=[P1',P2'];
        P = P(find(NDSort(P,1)==1),:);

    otherwise
        P = [];

end
end


function W = ReplicatePoint(SampleNum,M)
if M > 1
    SampleNum = (ceil(SampleNum^(1/M)))^M;
    Gap       = 0:1/(SampleNum^(1/M)-1):1;
    eval(sprintf('[%s]=ndgrid(Gap);',sprintf('c%d,',1:M)))
    eval(sprintf('W=[%s];',sprintf('c%d(:),',1:M)))
else
    W = (0:1/(SampleNum-1):1)';
end
end

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
end

function Output = convex(x)
Output = fliplr(cumprod([ones(size(x,1),1),1-cos(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),1-sin(x(:,end-1:-1:1)*pi/2)];
end

function Output = mixed(x)
Output = 1-x(:,1)-cos(10*pi*x(:,1)+pi/2)/10/pi;
end

function Output = disc(x)
Output = 1-x(:,1).*(cos(5*pi*x(:,1))).^2;
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


function P = InterX(L1,varargin)
%INTERX Intersection of curves
%   P = INTERX(L1,L2) returns the intersection points of two curves L1
%   and L2. The curves L1,L2 can be either closed or open and are described
%   by two-row-matrices, where each row contains its x- and y- coordinates.
%   The intersection of groups of curves (e.g. contour lines, multiply
%   connected regions etc) can also be computed by separating them with a
%   column of NaNs as for example
%
%         L  = [x11 x12 x13 ... NaN x21 x22 x23 ...;
%               y11 y12 y13 ... NaN y21 y22 y23 ...]
%
%   P has the same structure as L1 and L2, and its rows correspond to the
%   x- and y- coordinates of the intersection points of L1 and L2. If no
%   intersections are found, the returned P is empty.
%
%   P = INTERX(L1) returns the self-intersection points of L1. To keep
%   the code simple, the points at which the curve is tangent to itself are
%   not included. P = INTERX(L1,L1) returns all the points of the curve
%   together with any self-intersection points.
%
%   Example:
%       t = linspace(0,2*pi);
%       r1 = sin(4*t)+2;  x1 = r1.*cos(t); y1 = r1.*sin(t);
%       r2 = sin(8*t)+2;  x2 = r2.*cos(t); y2 = r2.*sin(t);
%       P = InterX([x1;y1],[x2;y2]);
%       plot(x1,y1,x2,y2,P(1,:),P(2,:),'ro')

%   Author : NS
%   Version: 3.0, 21 Sept. 2010

%   Two words about the algorithm: Most of the code is self-explanatory.
%   The only trick lies in the calculation of C1 and C2. To be brief, this
%   is essentially the two-dimensional analog of the condition that needs
%   to be satisfied by a function F(x) that has a zero in the interval
%   [a,b], namely
%           F(a)*F(b) <= 0
%   C1 and C2 exactly do this for each segment of curves 1 and 2
%   respectively. If this condition is satisfied simultaneously for two
%   segments then we know that they will cross at some point.
%   Each factor of the 'C' arrays is essentially a matrix containing
%   the numerators of the signed distances between points of one curve
%   and line segments of the other.

%...Argument checks and assignment of L2
error(nargchk(1,2,nargin));
if nargin == 1,
    L2 = L1;    hF = @lt;   %...Avoid the inclusion of common points
else
    L2 = varargin{1}; hF = @le;
end

%...Preliminary stuff
x1  = L1(1,:)';  x2 = L2(1,:);
y1  = L1(2,:)';  y2 = L2(2,:);
dx1 = diff(x1); dy1 = diff(y1);
dx2 = diff(x2); dy2 = diff(y2);

%...Determine 'signed distances'
S1 = dx1.*y1(1:end-1) - dy1.*x1(1:end-1);
S2 = dx2.*y2(1:end-1) - dy2.*x2(1:end-1);

C1 = feval(hF,D(bsxfun(@times,dx1,y2)-bsxfun(@times,dy1,x2),S1),0);
C2 = feval(hF,D((bsxfun(@times,y1,dx2)-bsxfun(@times,x1,dy2))',S2'),0)';

%...Obtain the segments where an intersection is expected
[i,j] = find(C1 & C2);
if isempty(i),P = zeros(2,0);return; end;

%...Transpose and prepare for output
i=i'; dx2=dx2'; dy2=dy2'; S2 = S2';
L = dy2(j).*dx1(i) - dy1(i).*dx2(j);
i = i(L~=0); j=j(L~=0); L=L(L~=0);  %...Avoid divisions by 0

%...Solve system of eqs to get the common points
P = unique([dx2(j).*S1(i) - dx1(i).*S2(j), ...
    dy2(j).*S1(i) - dy1(i).*S2(j)]./[L L],'rows')';

    function u = D(x,y)
        u = bsxfun(@minus,x(:,1:end-1),y).*bsxfun(@minus,x(:,2:end),y);
    end
end


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


function PopCon = Constraint7(PopObj)
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