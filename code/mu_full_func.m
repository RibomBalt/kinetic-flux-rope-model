function mu = mu_full_func(vr,vp,vz,xx,Bzfunc,Bphifunc,params)
    % vc = ExB + gradB, ignoring curvB and use local-Bz in mu of gradB
    [rc, tc] = rcenter3(vr,vp,vz,xx,Bzfunc,Bphifunc,params);

    Bz = Bzfunc(rc);
    Bp = Bphifunc(rc);
    ele = atan2(Bz, Bp);
    Bt = hypot(Bz,Bp);
    
%     [Bx,By] = pol2cart(tc+pi/2, Bp);
%     
%     vpara_co = (vr*Bx + vp*By + vz*Bz)/Bt^2;
%     vparax = vpara_co*Bx;
%     vparay = vpara_co*By;
%     vparaz = vpara_co*Bz;
% 
%     mu = ((vr - vparax).^2 + (vp - vparay).^2 + (vz - vparaz).^2)./Bt*params.Bst;
    
    [bx,by,bz] = sph2cart(tc + pi/2, ele, 1);
    vb = vr.*bx + vp.*by + vz.*bz;
    mu = ((vr - vb.*bx).^2 + (vp - vb.*by).^2 + (vz - vb.*bz).^2)./Bt*params.Bst;
end

function [rc, tc] = rcenter3(vr,vp,vz,xx,Bzfunc,Bphifunc,params)
    % [~,vperp] = cart2pol(vr, vp);
    Bz = Bzfunc(xx);
    Bp = Bphifunc(xx);
    Bt = hypot(Bz,Bp);
    co = params.me*params.v0e/params.e/(Bt*Bt*1e-9)/(params.Lst*1e-2);
    [tc, rc] = cart2pol(xx - co*Bz*vp + co*Bp*vz, +co*Bz*vr);
%     [tc, rc, zc] = cart2pol(xx - co*Bz*vp + co*Bp*vz, +co*Bz*vr, -co*Bp*vr);
end