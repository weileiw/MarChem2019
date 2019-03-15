%+++++ this script contains two functions,
% function [f,dfdx,d2fdx2,eigV,eigE,M] = cost(x,p,grd,M2d,EE,tt) and
% function [M,D] = Pcycle(p,PFD,dPFDdb,dPFDdbeta,M2d);
% this first function is used to calculated non steady-state
% objective function
% value(f), along with gradient(dfdx) and hessian(d2fdx2)
% towards parameters with the output of the second function.
% x is a vector of parameters, p is a
% structure that contains data and constants.M2d is a mask. EE is
% the eigen vector based on steady-state model, tt is the number of
% eigen vector with e-folding decaying time longer than 5 days.

function [f,dfdx,d2fdx2,eigV,eigE,M] = cost_2nd(x,p,grd,M2d,EE,tt)

    LAMBDA = p.LAMBDA;
    GAMMA  = p.GAMMA;
    THETA = p.THETA;
    
    nip    = length(x);
    d2fdx2 = zeros(nip,nip);
    dx     = sqrt(-1)*eps.^3*eye(nip);
    O = p.O;
    % prior = p.prior;
    prior =    [-0.11; -0.11;  0.32; 0.32; 0.32; 1.70; 5.00; ...
                2; -6.00; -6.00; 0.46; -6.00; -6.00];

    U   = d0([1/0.16; 1/0.16; 1/2.98; 1/2.98; 1/2.98; 1/4.01; 1/6.64; ...
              1/100; 1/100; 1/100; 1/100; 1/100; 1/100]);
    
    W = d0([1./(p.POC*0.02).^2;1./(p.poc*0.02).^2; 1./(p.Chl*0.30).^2;...
            1./(p.chl*0.30).^2;1./(p.Phyo*0.30).^2;1./(p.phyo*0.30).^2]);
    
    ap = 0*ones(tt,1);
    aU = d0(ones(tt,1));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % set parameter values back to its prior once the computer
    % suggests a weird numbers, such as Martin curve exponential b
    % of >600.
    xii_old = exp(p.prior);
    a = 5;
    xii_new = exp(x(1:13));
    tmp = find(abs(xii_new)>=abs(a*xii_old) | ...
               abs(xii_new)<= abs((1/a)*xii_old));

    for ii = 1:length(tmp)
        xii_new(tmp(ii)) = xii_old(tmp(ii));
    end
    
    % restrict b value below 1.5 and .5.
    if (xii_new(1)>1.5 | xii_new(1)<0.5);
        xii_new(1) = xii_old(1);
    end
    if (xii_new(2)>1.5 | xii_new(2)<0.5);
        xii_new(2) = xii_old(2);
    end
    x = [log(xii_new);x(14:end)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % the gradient is calculated analyticbally.
    % the for loop is used to calculated Hessian matrix by
    % usin{g complex step method.
  for ii = 1:nip
      x    = real(x)+dx(:,ii);
      p.bl  = exp(x(1));
      p.bs  = exp(x(2));
      p.kappa1 = exp(x(3));
      p.kappa2 = exp(x(4));
      p.kappa3 = exp(x(5));
      p.alpha  = exp(x(6));
      p.beta   = exp(x(7));
      p.xPOC   = exp(x(8));
      p.xChl   = exp(x(9));
      p.xPhy   = exp(x(10));
      p.xpoc   = exp(x(11));
      p.xchl   = exp(x(12));
      p.xphy   = exp(x(13));
      amp = x(14:end);

      [PFDl,dPFDldb,dPFDldbeta]  = buildPFD_l(M2d,p,grd);
      [PFDs,dPFDsdb,dPFDsdalpha] = buildPFD_s(M2d,p,grd);

      p.PFDl = PFDl; p.dPFDldb = dPFDldb; p.dPFDldbeta  = dPFDldbeta;
      p.PFDs = PFDs; p.dPFDsdb = dPFDsdb; p.dPFDsdalpha = dPFDsdalpha;
      
      [M,G] = Pcycle(p,M2d,EE);
      
      % ep = (x-prior);
      ep = x(1:13)-prior;
      ea = amp-ap;
      
      EM = EE*amp;  dem = EE;
      
      ed = M+EM-O;

      % f = 0.5*GAMMA*(ed.'*W*ed) + 0.5*LAMBDA*(ep.'*U*ep);
      f = 0.5*GAMMA*ed.'*W*ed + 0.5*LAMBDA*(ep.'*U*ep) + 0.5*THETA*(ea.'*aU*ea);
      
      dpdx  = diag(exp(x(1:13)));

      dedx = G*dpdx;
      dfdx = GAMMA*[ed.'*W*dedx, ed.'*W*dem] + [LAMBDA*ep.'*U, ...
                          THETA*ea.'*aU];
      % dfdx = GAMMA*[ed.'*W*dedx, ed.'*W*dem] + LAMBDA*ep.'*U;
      
      dfdx_test(ii) = imag(f)./eps.^3;
      d2fdx2(ii,:)  = imag(dfdx)./eps.^3;
      
  end
  
  f = real(f);
  dfdx = real(dfdx);
  
  POC_O = [p.POC;  p.poc];
  Chl_O = [p.Chl;  p.chl];
  Phy_O = [p.Phyo; p.phyo];
  jj = length(p.POC);
  M = [M+EM];
  POC_M = M(0*jj+1:2*jj);
  Chl_M = M(2*jj+1:4*jj);
  Phy_M = M(4*jj+1:6*jj);
  
%%%%%%%% Comment this out of model versus observation plot %%%%%%%%%%
% figure(1)
% loglog(POC_O,real(POC_M),'rp',[1e-8:10],[1e-8:10],'MarkerSize', 12,'MarkerFaceColor','r')
% hold on
% loglog(Chl_O,real(Chl_M),'co',[1e-8:10],[1e-8:10],'MarkerSize', 12,'MarkerFaceColor','c')
% hold on
% loglog(Phy_O,real(Phy_M),'b^',[1e-8:10],[1e-8:10],'MarkerSize', 12,'MarkerFaceColor','b')
% hold off
% xlim([1e-8,10]);
% ylim([1e-8,10]);
% O = [POC_O;Chl_O;Phy_O];
% r2 = rsquare(real(M),O);
% txt = sprintf('R^2 = %.2f',r2);
% text(0.01,1,txt)
% xlabel('Observation (\mumol L^-^1)','FontSize',12)
% ylabel('Model prediction (\mumol L^-^1)','FontSize',12)

% set(gca,'fontsize',12)
% set(findall(gca, 'Type', 'Line'),'LineWidth',2);
% exportfig(gcf,'MvsO_2nd_eig','fontsize',20,'fontmode','fixed','renderer','painters')
% keyboard
%%++++++++++++++++++++++++++++++++++++++++++++++
%%%%%%%% Comment this out of model versus observation plot
%%%%%%%% %%%%%%%%%%

%%%%%%++++++function [M,G] = Pcycle(p,PFD,dPFDdb,dPFDdbeta,M2d)
% function Pcycle output model field (M) and first derivative(D)
function [M,D] = Pcycle(p,M2d,EE)
    
    global GN
    iwet = find(M2d(:));
    nwet = length(iwet);
    I = speye(nwet);
    
    PFDl = p.PFDl;
    PFDs = p.PFDs;
    dPFDldb = p.dPFDldb;
    dPFDsdb = p.dPFDsdb;
    dPFDldbeta  = p.dPFDldbeta;
    dPFDsdalpha = p.dPFDsdalpha; 
    kappa1 = p.kappa1;
    kappa2 = p.kappa2;
    kappa3 = p.kappa3;
    alpha  = p.alpha;
    beta   = p.beta;
    xPOC   = p.xPOC;
    xChl   = p.xChl;
    xPhy   = p.xPhy;
    xpoc   = p.xpoc;
    xchl   = p.xchl;
    xphy   = p.xphy;
    
    POC = GN(0*nwet+1:1*nwet);
    poc = GN(1*nwet+1:2*nwet);
    Chl = GN(2*nwet+1:3*nwet);
    chl = GN(3*nwet+1:4*nwet);
    Phy = GN(4*nwet+1:5*nwet);
    phy = GN(5*nwet+1:6*nwet);
    
    X0 = [POC;poc;Chl;chl;Phy;phy];
    options.atol = 1e-12;
    options.rtol = 1e-12;
    options.iprint = 0;
      
    [M,ierr] = nsnew(X0,@(X) eqs(X,p,M2d),options);
    GN = 0.999*real(M);
    
    if (ierr ~=0)
        fprintf('cost eqNcycle did not converge.\n');
        fprintf('current alpha, beta, gamma %f,%f,%f.\n', p.LAMBDA, ...
                p.GAMMA, p.THETA);
        keyboard
    end
    
    if nargout>1
        POC = M(0*nwet+1:1*nwet);
        poc = M(1*nwet+1:2*nwet);
        Chl = M(2*nwet+1:3*nwet);
        chl = M(3*nwet+1:4*nwet);
        Phy = M(4*nwet+1:5*nwet);
        phy = M(5*nwet+1:6*nwet);
        
        Z = sparse(nwet,1);
        O = [1;Z(1:end-1)];
        dFdbl  = [-dPFDldb*POC;...
                  Z;...
                  -dPFDldb*Chl;...
                  Z;...
                  -dPFDldb*Phy;...
                  Z];
        
        dFdbs = [Z;...
                 -dPFDsdb*poc;...
                 Z;...
                 -dPFDsdb*chl;...
                 Z;...
                 -dPFDsdb*phy];
        
        dFdkappa1 = [Z;...
                     Z;...
                     Z;...
                     -chl;...
                     Z;...
                     chl];   
        
        dFdkappa2 = [Z;...
                     -poc;...
                     Z;...
                     Z;...
                     Z;...
                     Z];        
        
        dFdkappa3 = [Z;...
                     Z;...
                     Z;...
                     Z;...
                     Z;...
                     -phy];
        
        dFdalpha = [poc.*poc;...
                    -poc.*poc-dPFDsdalpha*poc;...
                    chl.*chl;...
                    -chl.*chl-dPFDsdalpha*chl;...
                    phy.*phy;...
                    -phy.*phy-dPFDsdalpha*phy];
        
        dFdbeta  = [-POC-dPFDldbeta*POC;...
                    POC;...
                    -Chl-dPFDldbeta*Chl;...
                    Chl;...
                    -Phy-dPFDldbeta*Phy;...
                    Phy];
        
        dFdxPOC = [O;...
                   Z;...
                   Z;...
                   Z;...
                   Z;...
                   Z];
        
        dFdxChl = [Z;...
                   Z;...
                   O;...
                   Z;...
                   Z;...
                   Z];
        
        dFdxPhy = [Z;...
                   Z;...
                   Z;...
                   Z;...
                   O;...
                   Z];

        dFdxpoc = [Z;...
                   O;...
                   Z;...
                   Z;...
                   Z;...
                   Z];
        
        dFdxchl = [Z;...
                   Z;...
                   Z;...
                   O;...
                   Z;...
                   Z];  

        dFdxphy = [Z;...
                   Z;...
                   Z;...
                   Z;...
                   Z;...
                   O];
        
        Fx = [dFdbl,...
              dFdbs,...
              dFdkappa1,...
              dFdkappa2,...
              dFdkappa3,...
              dFdalpha,...
              dFdbeta,...
              dFdxPOC,...
              dFdxChl,...
              dFdxPhy,...
              dFdxpoc,...
              dFdxchl,...
              dFdxphy];
        
        % dFdb,...
        Jac = [[-(beta*I+PFDl),      2*alpha*d0(poc), 0*I, 0*I, 0*I, 0*I];...
               [beta*I, -(2*alpha*d0(poc)+kappa2*I+PFDs), 0*I, 0*I, 0*I, 0*I];...
               [0*I, 0*I, -(beta*I+PFDl),      2*alpha*d0(chl), 0*I, 0*I];...
               [0*I, 0*I, beta*I, -(2*alpha*d0(chl)+kappa1*I+PFDs), 0*I, 0*I];...
               [0*I, 0*I, 0*I, 0*I,   -(beta*I+PFDl),    2*alpha*d0(phy)];...
               [0*I, 0*I, 0*I, kappa1*I, beta*I, -(2*alpha*d0(phy)+kappa3*I+PFDs)]];
        
        [eigV,eigE,FLAG] = eig(full(real(-Jac)));
        
        DF = mfactor(Jac);
        
        D = -mfactor(DF,Fx);
    end
    
function [F, DF] = eqs(X,p,M2d)
    
    iwet = find(M2d(:));
    nwet = length(iwet);
    I = speye(nwet);

    PFDs = p.PFDs; PFDl = p.PFDl;
    kappa1 = p.kappa1;
    kappa2 = p.kappa2;
    kappa3 = p.kappa3;
    alpha  = p.alpha;
    beta   = p.beta;
    xPOC   = p.xPOC;
    xChl   = p.xChl;
    xPhy   = p.xPhy;
    xpoc   = p.xpoc;
    xchl   = p.xchl;
    xphy   = p.xphy;

    POC = X(0*nwet+1:1*nwet);
    poc = X(1*nwet+1:2*nwet);
    Chl = X(2*nwet+1:3*nwet);
    chl = X(3*nwet+1:4*nwet);
    Phy = X(4*nwet+1:5*nwet);
    phy = X(5*nwet+1:6*nwet);

    fPOC = [xPOC;zeros(nwet-1,1)];
    fChl = [xChl;zeros(nwet-1,1)];
    fPhy = [xPhy;zeros(nwet-1,1)];

    fpoc = [xpoc;zeros(nwet-1,1)];
    fchl = [xchl;zeros(nwet-1,1)];
    fphy = [xphy;zeros(nwet-1,1)];
    
    F = [[alpha*poc.*poc-beta*POC-PFDl*POC+fPOC];...
         [beta*POC-alpha*poc.*poc-kappa2*poc-PFDs*poc+fpoc];...
         [alpha*chl.*chl-beta*Chl-PFDl*Chl+fChl];...
         [beta*Chl-alpha*chl.*chl-kappa1*chl-PFDs*chl+fchl];...
         [alpha*phy.*phy-beta*Phy-PFDl*Phy+fPhy];...
         [beta*Phy+kappa1*chl-alpha*phy.*phy-kappa3*phy-PFDs*phy+fphy]];

    if nargout > 1
        
        Jac = [[-(beta*I+PFDl),      2*alpha*d0(poc), 0*I, 0*I, 0*I, 0*I];...
               [beta*I, -(2*alpha*d0(poc)+kappa2*I+PFDs), 0*I, 0*I, 0*I, 0*I];...
               [0*I, 0*I, -(beta*I+PFDl),      2*alpha*d0(chl), 0*I, 0*I];...
               [0*I, 0*I, beta*I, -(2*alpha*d0(chl)+kappa1*I+PFDs), 0*I, 0*I];...
               [0*I, 0*I, 0*I, 0*I,   -(beta*I+PFDl),    2*alpha*d0(phy)];...
               [0*I, 0*I, 0*I, kappa1*I, beta*I, -(2*alpha*d0(phy)+kappa3*I+PFDs)]];
        
        DF = mfactor(Jac);
        
    end
    