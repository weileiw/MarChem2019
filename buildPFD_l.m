function [PFdiv,dPFDdb,dPFDdd] = buildPFD_l(M2d,p,grd);
% this code is used to build Particle Flux Diverngence (PFD)
%%%
%                      ________________________
%                     /         A             /|  POP sinking
%                 top/_______________________/ |       |
%                    |  |                    | |       | -w
%                    | /        V            | /       |
%                 bot|/______________________|/        V
%
% PFD = (A*w(top)*POP(top)-A*w(bot)*POP(bot))/V;
% add an exra layer of zeros at the bottom to ensure there is no
% flux in or out of the domain when using periodic shift operators
% for finite differences and averaging
    [nx,nz] = size(M2d);
    M2D = zeros(nx,nz+1);
    M2D(:,1:end-1) = M2d;
    % add the zw coordinate at the top of the extra layer
    ZW2d = grd.ZW2d;
    ZW2d = ZW2d(:,[1:end,end]);
    ZW2d(:,end) = grd.ZW2d(:,end)+grd.dzt(end);
    % areas of the top of the grid box
    %dAt = (grd.DXT3d.*grd.DYT3d).*M2d;
    % volume of the grid boxes
    %dVt = (grd.DXT3d.*grd.DZT3d).*M2d;
    dVt = grd.dzt;
    n  = nx*(nz+1);
    I0 = speye(n);
    i0 = zeros(nx,nz+1);
    i0(:) = 1:n;
    % periodic shifts OK because M3D has a layer of zeros on the bottom
    iu = i0(:,[nz+1,1:nz]); %use a periodic upward shift
    ib = i0(:,[2:nz+1,1]); % use a periodic downward shift
    IU = I0(iu,:);
    IB = I0(ib,:);
    % keep only wet boxes
    iocn = find(M2D(:));
    I0 = I0(iocn,:); I0 = I0(:,iocn);
    IU = IU(:,iocn); IU = IU(iocn,:);
    IB = IB(:,iocn); IB = IB(iocn,:);

    % (averages POP onto the top of the grid boxes)
    %AVG = (I0+IU)\(I0+IU);
    % (compute the divergence in the center of the grid boxes)
    DIV = d0(dVt)\(I0-IB);
    % (compute the flux at the top of the grid cells)
    % mimics a Martin curve flux attenuation profile
    %(see Kriest and Oschelies 2008 in Biogeosciences)
    r = p.beta;
    b = p.bl;
    a = r/b;
    % particle sinking velocity at the top of the grid cells.
    MSK = M2D.*M2D(:,[nz+1,1:nz]);
    M = MSK.*ZW2d;
    w = -a*M;
    
    dadb = -r/(b^2);
    dadr = 1/b;
    dwdb = -dadb*M;
    dwdr = -dadr*M;
    %FLUX = d0(w(iocn))*AVG;
    FLUX = d0(w(iocn))*IU;
    dFLUXdb = d0(dwdb(iocn))*IU;
    dFLUXdr = d0(dwdr(iocn))*IU;
    % particle flux divergence operator
    PFdiv = DIV*FLUX;
    dPFDdb = DIV*dFLUXdb;
    dPFDdr = DIV*dFLUXdr;
    dPFDdd = dPFDdr;
end
