function grd = buildgrd(p,M2d)

depth = p.dp;
tmp   = zeros(length(depth)+1,1);
tmp   = [depth;1000];% add a bottom layer,
                    %which represents bottom depth.
dzt   = tmp(2:end)-depth;
ZW2d  = M2d;

for ii = 1:length(depth)
    
  ZW2d(:,ii)  = depth(ii)-p.dp(1);

end
grd.dzt   = dzt;
grd.ZW2d  = ZW2d;


