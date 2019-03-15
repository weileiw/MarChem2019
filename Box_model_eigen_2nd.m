% this is the calling function and used to load data, set up
% constants and parameters. It also sets up optimization conditions
% for Matlab routine fminunc. At the end, it calculates errorbars
% for optimal parameter values and save them to file.
load DB2_May2005.mat
global GN
p.dp   = depth;     % sampling depths.
p.Chl  = Chla;      % large sized Chl a.
p.POC  = POC;       % large sized particulate organic matter.
p.Phyo = Phyeo;     % large sized phaeopigment.
p.chl  = chla;      % small sized Chl a.
p.poc  = poc;       % small sized particulated organic matter.
p.phyo = phyeo;     % small sized phaeopigment.

% observation data vector
p.O = [p.POC;p.poc;p.Chl;p.chl;p.Phyo;p.phyo];
GN  = p.O;
% second per year
spy = 365*24*60*60;

jj = length(depth);
M2d = ones(1,jj);

grd = buildgrd(p,M2d); % build model grid;
bs  = 0.84;    % Martin Curve exponential for small particles.
bl  = 0.84;    % Martin Curve exponential for large particles.
kappa1  = 1;   % small sized Chl a to phyeopigment rate constant [y^[-1]]
kappa2  = 1;   % small sized POC remineralization rate [y^[-1]];
kappa3  = 1;   % Phyeopigment reminearlization rate [y^[-1]];
alpha   = 3;   % aggregation rate      [y^[-1]];
beta    = 150; % disaggregation rate   [y^[-1]];
xPOC = 1; % inital guess of large-size POC flux;
xpoc = 1; % inital guess of small-size POC flux;
xChl = 1; % inital guess of large-size Chl a flux;
xchl = 1; % inital guess of small-size Chl a flux;
xPhy = 1; % inital guess of large-size pheopigment flux;
xphy = 1; % inital guess of small-size pheopigment flux;

% vector of intial parameter guess 
x0 = [bl;bs;kappa1;kappa2;kappa3;alpha;beta;xPOC;xpoc;xChl;xchl;xPhy;xphy];
x0  = log(x0); % lognormal transformation to ensure positive
               % parameter values;

nip = length(x0);

lambda = linspace(2,100,50);      % hyper parameter scaling parameters;
gamma  = linspace(0.001,0.01,10); % hyper parameter scaling data;
theta  = logspace(0,3,30);        % hyper parameter scaling amplitudes 
% theta = 1e-1;

% lambda = R.LAMBDA;
% gamma  = R.GAMMA;
% theta  = R.THETA;

[X,Y,Z] = meshgrid(lambda,gamma,theta);

logZ = 0*X;
ndata = length(p.O);
for jj = 1:length(X(:))

    if mod(jj,20) == 0
        fprintf('current iteration is %d \n',jj)
    end
    p.LAMBDA = X(jj);
    p.GAMMA  = Y(jj);
    p.THETA  = Z(jj);
    
    L = @(x) neglogpost_eigen(x,p,grd,M2d);

    options = optimoptions(@fminunc,...
                           'Algorithm','trust-region',...
                           'GradObj','on',...
                           'Hessian','on',...
                           'Display','off',...
                           'MaxFunEvals',100,...
                           'MaxIter',100,...
                           'TolX',1e-16,...
                           'DerivativeCheck','off',...
                           'FinDiffType','central',...
                           'TolFun',1e-16,...
                           'PrecondBandWidth',Inf);

    % optimize model under steady-state 
    [xhat,fval,exitflag] = fminunc(L,x0,options);
    % extract the slowly decaying modes and save
    % them to a file
    [f,dfdx,d2fdx2,eigV,eigE,M] = L(xhat);
    ev = diag(eigE);
    tr = real(1./ev)*365;    % e-folding decay in days
    ti = 2*pi./imag(ev)*365; % period in days
    
    [~,ii] = sort(tr,'descend');
    tr = tr(ii);
    ti = ti(ii);

    E = eigV(:,ii);
    ikeep = find(tr>5); % keep modes that have an e-folding decay
                        % time-scale of more than 5 days
    tr = tr(ikeep);
    ti = ti(ikeep);
    E  = E(:,ikeep);

    EE = zeros(size(E,1),length(tr));
    for ii = 1:length(tr)
        if (isreal(E(:,ii)))
            EE(:,ii)  = real(E(:,ii))/max(abs(real(E(:,ii))));
        else
            EE(:,ii)   = real(E(:,ii))/max(abs(real(E(:,ii))));
            EE(:,ii+1) = imag(E(:,ii))/max(imag(E(:,ii)));
            ii = ii+1;
        end
    end

    tt = length(ikeep);
    
    x1 = [xhat; zeros(tt,1)];
    p.prior = xhat;
    % optimize model with slowly-decaying modes.
    F = @(x) cost_2nd(x,p,grd,M2d,EE,tt);
    [xhat,fval,exitflag] = fminunc(F,x1,options);
    
    [f,dfdx,d2fdx2] = cost_2nd(xhat,p,grd,M2d,EE,tt);
    HH = nearestSPD(d2fdx2);
    
    logZ(jj) = -f-0.5*log(det(HH))+...
        (nip/2)*log(p.LAMBDA)+...
        (tt/2)*log(p.THETA)+...
        (ndata/2)*log(p.GAMMA);
    
    Q(jj).HH = HH;
    Q(jj).xhat = xhat;
    Q(jj).f = f;
    
end

%%%%%%%% Comment this out to get a contour plot for alpha and beta %%%%%%%%%%%%
%figure(2)
% contourf(X,Y,logZ);colorbar
%set(gca,'XTick',[0:1:10])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6','7','8','9','10'})
%set(gca,'YTick',[20:5:60])
%set(gca,'YTickLabel',{'20','25','30','35','40','45','50','55','60'})
%text('Interpreter','latex','String','$$\Lambda = 21.89$$','Position',[50 35],'fontsize',16)
%text('Interpreter','latex','String','$$\Gamma = 27.76$$','Position',[50 37],'fontsize',16)
%text(5,55,'MedFlux')
% xlabel('\Lambda-scaling factor for parameter')
% ylabel('\Gamma-scaling factor for data')
% set(gca,'fontsize',16)
%%%%%%%%  Comment this out to get a contour plot for alpha and beta %%%%%%%%%%%

% finding optimal alpha and beta, recalculate  Hessian matrix based
% on them, and calculate error bars for each parameter values.
imax = find(logZ(:)==max(logZ(:)));
p.LAMBDA = X(imax);
p.GAMMA  = Y(imax);
p.THETA  = Z(imax);
HH = Q(imax).HH;
xhat = Q(imax).xhat;

N_p = nip-p.LAMBDA*(trace(HH)).^(-1);
error = sqrt(diag(inv(HH)));

R.X = X;
R.Y = Y;
R.logZ = logZ;

% calculate parameters errorbars
R.upbar = [(exp(xhat(1:13)+error(1:13))-exp(xhat(1:13))); error(14:end)];
R.lowbar = [(exp(xhat(1:13))-exp(xhat(1:13)-error(1:13))); error(14:end)];
R.xhat = [exp(xhat(1:13));xhat(14:end)];
R.LAMBDA = p.LAMBDA;
R.GAMMA = p.GAMMA;
R.THETA = p.THETA;

fname = sprintf('xhat_theta_%1.0e',p.THETA);
% fname = sprintf('xhat_mode_%s',date);
save(fname,'R');

fprintf('----done----- \n \n')