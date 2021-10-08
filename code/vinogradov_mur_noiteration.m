%% 1. params
% MODEL PARAMETERS
params.alpha_rot = 0.75;
params.alpha_z = 0.35;
% ga: uz/(Omegae*Lst)
params.ga = 2;
% epsilon: Omegae/omega_ge
params.epsilon = 0.2e-2;
% bz nT at center
params.bz = 22;
% n0 cm-3
params.n0i = 40;
% bebg: background mu
params.bebg = -0.07;
params.berot = -0.;
params.bez = 0.53;
% taub Teb/Te >1
params.teb = 34.8;
params.taub = 1.05;
% Te eV
params.te = params.teb/params.taub;
% tau Ti/Te
params.tau = 3.7;
params.ti = params.te*params.tau;
% tauz Tez/Te
params.tauz = 0.4;
params.tez = params.te*params.tauz;


% PHYSICAL CONSTANT
% mp, kg
params.mp=1.6726231e-27;
% me/mp
params.m0=5.196e-4;
% c m/s
params.c=2.99792458e8;
% e Coulomb
params.e=1.60218e-19;
params.egau = params.e*params.c*1e1;
% mp,me kg
params.mp=1.6726231e-27;
params.me=9.1093897e-31;
% mu0 
params.mu0 = 4*pi*1e-7;
% OUTPUT PARAMETERS
% Bst nT
params.Bst = sqrt(4*pi*params.n0i*(params.te*params.e*1e7))*1e5;
% cs: ion sound speed cm/s
params.cs = sqrt(params.te*params.e/params.mp)*1e2;
% Omegae: s-1
params.Omegae = params.epsilon*(params.egau*(params.Bst*1e-5)/(params.me*1e3)/(params.c*1e2));
% di: ion inertial length cm
params.di = (params.c*1e2)/sqrt(4*pi*params.n0i*params.egau^2/(params.mp*1e3));
% Lst: spatial scale: cm
params.Lst = sqrt(params.di*params.cs/params.Omegae);
% E0 mV/m
params.E0 = (params.te*params.e)/params.e/(params.Lst*1e-2)*1e3;
% n0i cm-3
params.n0e = params.n0i / ((1 + params.alpha_rot - params.alpha_z) +...
    params.alpha_z - params.alpha_rot);

%% 2. solver & model moments
inputs = 0:0.05:8;
[y,yp] = solver(params,inputs);

x = inputs;
Bz = (params.epsilon + y(2,:))*params.Bst;
Bphi = -yp(3,:)*params.Bst;
E = -yp(5,:)*params.E0;

xie = (1 - params.bebg./(y(2,:) + params.epsilon))/(1 - params.bebg*params.Bst/params.bz);
xier = (1 - params.berot./(y(2,:) + params.epsilon))/(1 - params.berot*params.Bst/params.bz);
xiez = (1 - params.bez./(y(2,:) + params.epsilon))/(1 - params.bez*params.Bst/params.bz);
Neb = (1+params.alpha_rot-params.alpha_z)./xie.*exp(y(5,:)/params.taub);
Ner = - params.alpha_rot./xier.*exp(y(5,:) - y(1,:));
Nez = params.alpha_z./xiez.*exp((y(5,:)-params.ga*y(3,:))/params.tauz);
Ne = params.n0e*(Neb + Ner + Nez);
% Ne = params.n0e*( (1+params.alpha_rot-params.alpha_z)./xie.*exp(y(5,:)/params.taub) ...
%     - params.alpha_rot*exp(y(5,:) - y(1,:)) + params.alpha_z*exp((y(5,:)-params.ga*y(3,:))/params.tauz));
Ni = params.n0i*(exp(-y(5,:)/params.tau));
% vz,vphi km/s
vphi = (params.Omegae*x*params.Lst*1e-5).*Ner*params.n0e./Ne;
vz = (params.ga*params.Lst*1e-5*params.Omegae)*Nez*params.n0e./Ne;


% Pb Pe Pt nPa
% Pb = (Bz*1e-9).^2/(2*params.mu0)*1e9;
% Pe = (params.e*params.te)*(1e6*params.n0e)*(params.taub*(1+params.alpha_rot)*exp(y(5,:)/params.taub) ...
%     - params.alpha_rot*exp(y(5,:) - y(1,:)) )*1e9;
% 

%% 3. distribution
% TODO
params.v0e = sqrt(2*params.te*params.e/params.me);
% fe_dist = @(vr,vp,vz,ii)(params.n0e*1e6)*(...
%     (1+params.alpha_rot-params.alpha_z)*(1/pi/params.taub)^1.5*exp(-1/params.taub*(vr.^2+vp.^2+vz.^2)).*exp(y(5,ii)/params.taub)...
%     + params.alpha_z*(1/pi/params.tauz)^1.5*exp((-(vr.^2+vp.^2+(vz - params.ga*params.Lst*1e-2*params.Omegae/params.v0e).^2)...
%         +y(5,ii)-params.ga*y(3,ii))/params.tauz)...
%     - params.alpha_rot*(1/pi)^1.5*exp(-(vr.^2+(vp - params.Omegae*x(ii)*params.Lst*1e-2/params.v0e).^2+vz.^2)+y(5,ii)-y(1,ii)));

psi_interp = griddedInterpolant(x, y(1,:));
a_interp = griddedInterpolant(x, y(3,:));
phi_interp = griddedInterpolant(x, y(5,:));
xie_interp = griddedInterpolant(x, xie);

Bzfunc = griddedInterpolant(x, Bz);
Bzfunc.ExtrapolationMethod = 'nearest';
E_interp = griddedInterpolant(x, E);
Erfunc = @(xx)(xx<=x(end)).*E_interp(xx);
gradBz_interp = griddedInterpolant([0 x(1:end-1) + diff(x)/2], [0 diff(Bz)./diff(x)]);
gradBzfunc = @(xx)(xx<=x(end)).*gradBz_interp(xx);

% fe_dist_interp = @(vr,vp,vz,x)(params.n0e*1e6)*(...
% (1+params.alpha_rot-params.alpha_z)*(1/pi/params.taub)^1.5*exp(-1/params.taub*((vr.^2+vp.^2).*xie_interp(x)+vz.^2)).*exp(phi_interp(x)/params.taub)...
% + params.alpha_z*(1/pi/params.tauz)^1.5*exp((-(vr.^2+vp.^2+(vz - params.ga*params.Lst*1e-2*params.Omegae/params.v0e).^2)...
%     +phi_interp(x)-params.ga*a_interp(x))/params.tauz)...
% - params.alpha_rot*(1/pi)^1.5*exp(-(vr.^2+(vp - params.Omegae*x*params.Lst*1e-2/params.v0e).^2+vz.^2)+phi_interp(x)-psi_interp(x)));

%% 4. plot
% h = irf_plot(5,'newfigure');
% ax = h(1);
% irf_plot(ax, [x;Bz;-Bphi;hypot(Bz,Bphi)]')
% ax = h(2);
% irf_plot(ax, [x;Ni;Ne]')
% ax = h(3);
% irf_plot(ax, [x;vphi;vz]')
% ax = h(4);
% irf_plot(ax, [x; (Neb*params.taub./(1 - params.bebg*params.Bst./Bz)...
%     + Ner./(1 - params.berot*params.Bst./Bz)...
%     + Nez*params.tauz)*(params.te * params.n0e)./Ne; ...
%     (Neb*params.taub + Ner + Nez*params.tauz)*(params.te * params.n0e)./Ne]')
% irf_legend(ax, {'perp','para'},[0.99 0.99])
% ax = h(5);
% irf_plot(ax, [x; Neb; Ner; Nez; Neb+Ner+Nez]');
fprintf('Bzpeak:%.2f; Bzinf:%.2f; Bphipeak:%.2f\n', Bz(1),Bz(end),min(Bphi))
fprintf('Nemax:%.2f; Nemin:%.2f; Necenter:%.2f\n', Ne(end), min(Ne) ,Ne(1))
fprintf('scale:%.2f;\n', (params.Lst*1e-5))

%% function
function [y,yp] = solver(params,inputs)
    opts = odeset('Mass',diag([1 1 1 1 0]),'RelTol',1e-4,'AbsTol',1e-10);
    sol = ode23t(@(x,y)odefun(x,y,params),...
        inputs, [0;params.bz/params.Bst - params.epsilon;0;0;0],opts);
    [y,yp] = deval(sol,inputs);
end

function dydx = odefun(x,y,params)
    xie = (1 - params.bebg/(y(2) + params.epsilon));
    xier = (1 - params.berot/(y(2) + params.epsilon));
    xiez = (1 - params.bez/(y(2) + params.epsilon));
    % DAE
    dydx = zeros([5 1]);

    % psi
    dydx(1) = y(2)*x;
    % dydx(2) = -x*params.alpha_rot.*exp(y(5) - y(1));
    dydx(2) = -params.alpha_rot*exp(y(5) - y(1))*(1 - params.berot*params.Bst/params.bz)/xier^2*x*exp(params.me/2/params.te*(params.Omegae*x*params.Lst*1e-2)^2*(1 - 1/xier));
    % z
    if x~=0
        dydx(3) = y(4)/x;
    else
        dydx(3) = 0;
    end
    dydx(4) = x*params.ga*params.alpha_z*(1 - params.bez*params.Bst/params.bz)/xiez*exp((y(5)-params.ga*y(3))/params.tauz);
    % phi: e - i
    dydx(5) = (1+params.alpha_rot-params.alpha_z)*(1 - params.bebg*params.Bst/params.bz)/xie*exp(y(5)/params.taub)...
        - params.alpha_rot*(1 - params.berot*params.Bst/params.bz)/xier*exp(y(5) - y(1))...
        + params.alpha_z*(1 - params.bez*params.Bst/params.bz)/xiez*exp((y(5) - params.ga*y(3))/params.tauz)...
         - params.n0i/params.n0e*exp(-y(5)/params.tau);
    if abs(dydx(5)) < 2*eps
        dydx(5) = 0;
    end
end

function y=dif2(x,t)
    x_diff=diff(x)./diff(t);
    y = interp1((t(1:end-1)+t(2:end))/2,x_diff,t,'spline',0);
end
