%% 0. get initial conditions
if ~exist('yp','var') % continue
    vinogradov_mur_noiteration;
end

%% 1. iterative process
% during iteration, assume Bphi=0. mu = m(vr^2+vphi^2)/Bz

% xie_interp = griddedInterpolant(x, xie);
% fe_dist_interp = @(vr,vp,vz,x)(params.n0e*1e6)*(...
% (1+params.alpha_rot-params.alpha_z)*(1/pi/params.taub)^1.5*exp(-1/params.taub*((vr.^2+vp.^2).*xie_interp(x)+vz.^2)).*exp(phi_interp(x)/params.taub)...
% + params.alpha_z*(1/pi/params.tauz)^1.5*exp((-(vr.^2+vp.^2+(vz - params.ga*params.Lst*1e-2*params.Omegae/params.v0e).^2)...
%     +phi_interp(x)-params.ga*a_interp(x))/params.tauz)...
% - params.alpha_rot*(1/pi)^1.5*exp(-(vr.^2+(vp - params.Omegae*x*params.Lst*1e-2/params.v0e).^2+vz.^2)+phi_interp(x)-psi_interp(x)));

psi_interp = griddedInterpolant(x, y(1,:));
a_interp = griddedInterpolant(x, y(3,:));
phi_interp = griddedInterpolant(x, y(5,:));

yfunc = cell(1,5);
% extrapolant: 1/2Bzr^2
yfunc{1} = @(xx)(xx<=x(end)).*psi_interp(xx) +...
    (xx>x(end)).*(y(1,end) + y(2,end)/2*(xx.^2 - x(end)^2));
% extrapolant: Bphi = mu0I/(2pi r)
yfunc{3} = @(xx)y3func(xx,x,y,a_interp);
phi_interp.ExtrapolationMethod = 'nearest';
yfunc{5} = phi_interp;
    
Bzfunc = griddedInterpolant(x, Bz);
Bzfunc.ExtrapolationMethod = 'nearest';
Bphi_func = griddedInterpolant(x, Bphi);
Bphifunc = @(xx)bpfunc(xx,x,Bphi_func);

E_interp = griddedInterpolant(x, E);
Erfunc = @(xx)(xx<=x(end)).*E_interp(xx);
gradBz_interp = griddedInterpolant([0 x(1:end-1) + diff(x)/2], [0 diff(Bz)./diff(x)]);
gradBzfunc = @(xx)(xx<=x(end)).*gradBz_interp(xx);

niteration = 3;
% x,y,yp
iterate_data = cell(4, 1+niteration);
iterate_data{1,1} = x;
iterate_data{2,1} = y;
iterate_data{3,1} = yp;
% f => EB
% feb = @(vr,vp,xx)exp(-1/params.taub*((vr.^2+vp.^2) - params.bebg*mufunc(vr,vp,xx,Bzfunc,Erfunc,gradBzfunc,params)));
% fer = @(vr,vp,xx)exp(-(vr.^2+(vp - params.Omegae*xx*params.Lst*1e-2/params.v0e).^2 - params.berot*mufunc(vr,vp,xx,Bzfunc,Erfunc,gradBzfunc,params)));
% % fez has same vr-vp distribution as feb
% fez = @(vr,vp,xx)exp(-1/params.tauz*((vr.^2+vp.^2) - params.bez*mufunc(vr,vp,xx,Bzfunc,Erfunc,gradBzfunc,params)));
fe3 = @(vr,vp,vz,xx,yy)(...
        (1+params.alpha_rot-params.alpha_z)*(1 - params.bebg*params.Bst/params.bz)*(1/pi/params.taub)^1.5*exp(-1/params.taub*((vr.^2+vp.^2+vz.^2) - params.bebg*mu_full_func(vr,vp,vz,xx,Bzfunc,Bphifunc,params))).*exp(yy(5)/params.taub)...
        + params.alpha_z*(1 - params.bez*params.Bst/params.bz)*(1/pi/params.tauz)^1.5*exp(-1/params.tauz*(vr.^2+vp.^2+(vz - params.ga*params.Lst*1e-2*params.Omegae/params.v0e).^2 - params.bez*mu_full_func(vr,vp,vz,xx,Bzfunc,Bphifunc,params))).*exp((yy(5) - params.ga*yy(3))/params.tauz)...
        - params.alpha_rot*(1 - params.berot*params.Bst/params.bz)*(1/pi)^1.5*exp(-(vr.^2+(vp - params.Omegae*xx*params.Lst*1e-2/params.v0e).^2 + vz.^2 - params.berot*mu_full_func(vr,vp,vz,xx,Bzfunc,Bphifunc,params))).*exp(yy(5) - yy(1)) );
% iterate_data{4,1} = feb;
% iterate_data{5,1} = fer;
% iterate_data{6,1} = fez;
iterate_data{4,1} = fe3;
for ii = 1:niteration
    
    % solve initial condition consistency
    integral_span = 7.2*max([params.taub params.taub/(1 - params.bebg*params.Bst/15) ...
        params.tauz params.tauz/(1 - params.bez*params.Bst/15) ...
        1 1/(1 - params.berot*params.Bst/15)]);
    % neb0_part = integral2(@(vr,vp)feb(vr,vp,0), -integral_span,integral_span,-integral_span,integral_span,'AbsTol',1e-7,'RelTol',1e-4);
    % ner0_part = integral2(@(vr,vp)fer(vr,vp,0), -integral_span,integral_span,-integral_span,integral_span,'AbsTol',1e-7,'RelTol',1e-4);
    % nez0_part = integral2(@(vr,vp)fez(vr,vp,0), -integral_span,integral_span,-integral_span,integral_span,'AbsTol',1e-7,'RelTol',1e-4);
    params.n0e = params.n0i / (integral3(@(vr,vp,vz)fe3(vr,vp,vz,0,y(:,1)),-integral_span,integral_span,0,integral_span,-integral_span,integral_span,'AbsTol',1e-7,'RelTol',1e-4)+...
        integral3(@(vr,vp,vz)fe3(vr,vp,vz,0,y(:,1)),-integral_span,integral_span,-integral_span,0,-integral_span,integral_span,'AbsTol',1e-7,'RelTol',1e-4) );
    
    fprintf('- n0i: %.2f; n0e: %.2f\n',params.n0i, params.n0e);
    opts = odeset('Mass',diag([1 1 1 1 0]),'AbsTol',1e-7,'RelTol',1e-4);
    % solve problem
    sol = ode15s(@(x,y)iterate_odefun(x,y,yfunc,fe3,params,integral_span),...
        inputs, [0;params.bz/params.Bst - params.epsilon;0;0;0],opts);
    [y,yp] = deval(sol,inputs);
    % cache y
    x = inputs;
    iterate_data{1,1+ii} = x;
    iterate_data{2,1+ii} = y;
    iterate_data{3,1+ii} = yp;
    % calculate EB, interps
    Bz = (params.epsilon + y(2,:))*params.Bst;
    Bphi = -yp(3,:)*params.Bst;
    E = -yp(5,:)*params.E0;

    psi_interp = griddedInterpolant(x, y(1,:));
    a_interp = griddedInterpolant(x, y(3,:));
    phi_interp = griddedInterpolant(x, y(5,:));

    % extrapolant: 1/2Bzr^2
    yfunc{1} = @(xx)(xx<=x(end)).*psi_interp(xx) +...
    (xx>x(end)).*(y(1,end) + y(2,end)/2*(xx.^2 - x(end)^2));
    % extrapolant: Bphi = mu0I/(2pi r)
    yfunc{3} = @(xx)y3func(xx,x,y,a_interp);
    phi_interp.ExtrapolationMethod = 'nearest';
    yfunc{5} = phi_interp;

    Bzfunc = griddedInterpolant(x, Bz);
    Bzfunc.ExtrapolationMethod = 'nearest';
    Bphi_func = griddedInterpolant(x, Bphi);
    Bphifunc = @(xx)bpfunc(xx,x,Bphi_func);

    E_interp = griddedInterpolant(x, E);
    Erfunc = @(xx)(xx<=x(end)).*E_interp(xx);
    gradBz_interp = griddedInterpolant([0 x(1:end-1) + diff(x)/2], [0 diff(Bz)./diff(x)]);
    gradBzfunc = @(xx)(xx<=x(end)).*gradBz_interp(xx);
    
    % feb = @(vr,vp,xx)exp(-1/params.taub*((vr.^2+vp.^2) - params.bebg*mufunc(vr,vp,xx,Bzfunc,Erfunc,gradBzfunc,params)));
    % iterate_data{4,ii+1} = feb;
    % fer = @(vr,vp,xx)exp(-(vr.^2+(vp - params.Omegae*xx*params.Lst*1e-2/params.v0e).^2 - params.berot*mufunc(vr,vp,xx,Bzfunc,Erfunc,gradBzfunc,params)));
    % iterate_data{5,ii+1} = fer;
    % fez = @(vr,vp,xx)exp(-1/params.tauz*((vr.^2+vp.^2) - params.bez*mufunc(vr,vp,xx,Bzfunc,Erfunc,gradBzfunc,params)));
    % iterate_data{6,ii+1} = fez;
    fe3 = @(vr,vp,vz,xx,yy)(...
        (1+params.alpha_rot-params.alpha_z)*(1 - params.bebg*params.Bst/params.bz)*(1/pi/params.taub)^1.5*exp(-1/params.taub*((vr.^2+vp.^2+vz.^2) - params.bebg*mu_full_func(vr,vp,vz,xx,Bzfunc,Bphifunc,params))).*exp(yy(5)/params.taub)...
        + params.alpha_z*(1 - params.bez*params.Bst/params.bz)*(1/pi/params.tauz)^1.5*exp(-1/params.tauz*(vr.^2+vp.^2+(vz - params.ga*params.Lst*1e-2*params.Omegae/params.v0e).^2 - params.bez*mu_full_func(vr,vp,vz,xx,Bzfunc,Bphifunc,params))).*exp((yy(5) - params.ga*yy(3))/params.tauz)...
        - params.alpha_rot*(1 - params.berot*params.Bst/params.bz)*(1/pi)^1.5*exp(-(vr.^2+(vp - params.Omegae*xx*params.Lst*1e-2/params.v0e).^2 + vz.^2 - params.berot*mu_full_func(vr,vp,vz,xx,Bzfunc,Bphifunc,params))).*exp(yy(5) - yy(1)) );
    iterate_data{4,ii+1} = fe3;
    fprintf('- iteration #%d complete\n', ii)
    disp(datestr(now));
end
disp('iteration finished')

%% 2. finalize quantities
disp('Start calc integrals...')
int_struct = calc_integral_parts(fe3,x,y,integral_span);
disp('calc integrals complete')
disp(datestr(now));
%%
calc_linemodel_func = @calc_linemodel;
[final_line, final_integrals] = calc_linemodel_func(x,y,yp,params,'intparts',int_struct);
disp('finalize complete')

%% 
fe3b = @(vr,vp,vz,xx,yy)(1+params.alpha_rot-params.alpha_z)*(1 - params.bebg*params.Bst/params.bz)*(1/pi/params.taub)^1.5*exp(-1/params.taub*((vr.^2+vp.^2+vz.^2) - params.bebg*mu_full_func(vr,vp,vz,xx,Bzfunc,Bphifunc,params))).*exp(yy(5)/params.taub);
fe3r = @(vr,vp,vz,xx,yy)-params.alpha_rot*(1 - params.berot*params.Bst/params.bz)*(1/pi)^1.5*exp(-(vr.^2+(vp - params.Omegae*xx*params.Lst*1e-2/params.v0e).^2 + vz.^2 - params.berot*mu_full_func(vr,vp,vz,xx,Bzfunc,Bphifunc,params))).*exp(yy(5) - yy(1));
fe3z = @(vr,vp,vz,xx,yy)params.alpha_z*(1 - params.bez*params.Bst/params.bz)*(1/pi/params.tauz)^1.5*exp(-1/params.tauz*(vr.^2+vp.^2+(vz - params.ga*params.Lst*1e-2*params.Omegae/params.v0e).^2 - params.bez*mu_full_func(vr,vp,vz,xx,Bzfunc,Bphifunc,params))).*exp((yy(5) - params.ga*yy(3))/params.tauz);
Neb = zeros(size(x));
Ner = zeros(size(x));
Nez = zeros(size(x));
parfor ii = 1:length(x)
    Neb(ii) = integral3(@(vr,vp,vz)fe3b(vr,vp,vz,x(ii),y(:,ii)), -integral_span,integral_span,-integral_span,integral_span,-integral_span,integral_span,'AbsTol',1e-7,'RelTol',1e-4);
    Ner(ii) = integral3(@(vr,vp,vz)fe3r(vr,vp,vz,x(ii),y(:,ii)), -integral_span,integral_span,-integral_span,integral_span,-integral_span,integral_span,'AbsTol',1e-7,'RelTol',1e-4);
    Nez(ii) = integral3(@(vr,vp,vz)fe3z(vr,vp,vz,x(ii),y(:,ii)), -integral_span,integral_span,-integral_span,integral_span,-integral_span,integral_span,'AbsTol',1e-7,'RelTol',1e-4);
end
final_line.Neb = Neb;
final_line.Ner = Ner;
final_line.Nez = Nez;
disp('calc parts complete')
% plot(x,final_line.vphi,x,final_line.vExBphi,x,final_line.vgradPphi,x,final_line.vExBphi + final_line.vgradPphi);
% legend

%% function
function dydx = iterate_odefun(x, y, yfunc, fe3, params, integral_span)
    % feb = @(vr,vp,xx)exp(-1/params.taub*((vr.^2+vp.^2) - params.bebg*mufunc(vr,vp,xx)));
    % integral
    % neb_part = integral2(@(vr,vp)feb(vr,vp,x), -integral_span,integral_span,-integral_span,integral_span,'RelTol',1e-8,'AbsTol',1e-11);
    % ueb_part = integral2(@(vr,vp)vp.*feb(vr,vp,x), -integral_span,integral_span,-integral_span,integral_span,'AbsTol',1e-7,'RelTol',1e-4);
    % ner_part = integral2(@(vr,vp)fer(vr,vp,x), -integral_span,integral_span,-integral_span,integral_span,'RelTol',1e-8,'AbsTol',1e-11);
    % uer_part = integral2(@(vr,vp)vp.*fer(vr,vp,x), -integral_span,integral_span,-integral_span,integral_span,'AbsTol',1e-7,'RelTol',1e-4);
    % nez_part = integral2(@(vr,vp)fez(vr,vp,x), -integral_span,integral_span,-integral_span,integral_span,'RelTol',1e-8,'AbsTol',1e-11);
    % uez_part = integral2(@(vr,vp)vp.*fez(vr,vp,x), -integral_span,integral_span,-integral_span,integral_span,'AbsTol',1e-7,'RelTol',1e-4);
    ne = integral3(@(vr,vp,vz)fe3(vr,vp,vz,x,y), -integral_span,integral_span, -integral_span,integral_span, -integral_span,integral_span,'AbsTol',1e-7,'RelTol',1e-4);
    if x>0.1
        ue = integral3(@(vr,vp,vz)vp.*fe3(vr,vp,vz,x,y), -integral_span,integral_span, -integral_span,integral_span, -integral_span,integral_span,'AbsTol',1e-7,'RelTol',1e-4);
    else
        ue = integral3(@(vr,vp,vz)vp.*fe3(vr,vp,vz,x,y), -integral_span,integral_span, 0,integral_span, -integral_span,integral_span,'AbsTol',1e-7,'RelTol',1e-4) + ...
            integral3(@(vr,vp,vz)vp.*fe3(vr,vp,vz,x,y), -integral_span,integral_span, -integral_span, 0, -integral_span,integral_span,'AbsTol',1e-7,'RelTol',1e-4);
    end
    uz = integral3(@(vr,vp,vz)vz.*fe3(vr,vp,vz,x,y), -integral_span,integral_span, -integral_span,integral_span, -integral_span,integral_span,'AbsTol',1e-7,'RelTol',1e-4);
    dydx = zeros(5,1);
    % psi
    dydx(1) = y(2)*x;
    if x==0
        dydx(2) = 0; % prevent error in integrals
    else
        dydx(2) = ue/(params.Omegae*params.Lst*1e-2/params.v0e);
    end
    % z
    if x~=0
        dydx(3) = y(4)/x;
    else
        dydx(3) = 0;
    end
    % TODO: mu on Pz should effect this term
    dydx(4) = x*uz/(params.Omegae*params.Lst*1e-2/params.v0e);
    % dydx(4) = params.alpha_z*(1/pi/params.tauz)*exp((y(5)-params.ga*y(3))/params.tauz)*nez_part*(x*params.ga);
    % phi: e - i
    % xie = (1 - params.bebg/(y(2) + params.epsilon));
    % feb = @(vr,vp,xx)exp(-1/params.taub*((vr.^2+vp.^2) - params.bebg*mufunc(vr,vp,xx)));
    dydx(5) = params.n0e*ne - params.n0i*exp(-y(5)/params.tau);

    

%     if (abs(dydx(5)) < eps) && (x<0.1)
%         dydx(5) = 0;
%     end
%     dydx(5) = dydx(5)*exp(-x/5);
end

function intparts = calc_integral_parts(fe3_,x,y,integral_span)
    neb_part = zeros(size(x));
    ueb_part = zeros(size(x));
    uez_part = zeros(size(x));
    Ppp_part = zeros(size(x));
    Prr_part = zeros(size(x));
    Pzz_part = zeros(size(x));
    parfor ii = 1:length(x)
        fe3 = fe3_;
        neb_part(ii) = integral3(@(vr,vp,vz)fe3(vr,vp,vz,x(ii),y(:,ii)), -integral_span,integral_span,-integral_span,integral_span,-integral_span,integral_span,'AbsTol',1e-7,'RelTol',1e-4);
        ueb_part(ii) = integral3(@(vr,vp,vz)vp.*fe3(vr,vp,vz,x(ii),y(:,ii)), -integral_span,integral_span,-integral_span,integral_span,-integral_span,integral_span,'AbsTol',1e-7,'RelTol',1e-4);
        uez_part(ii) = integral3(@(vr,vp,vz)vz.*fe3(vr,vp,vz,x(ii),y(:,ii)), -integral_span,integral_span,-integral_span,integral_span,-integral_span,integral_span,'AbsTol',1e-7,'RelTol',1e-4);
        Ppp_part(ii) = integral3(@(vr,vp,vz)vp.^2.*fe3(vr,vp,vz,x(ii),y(:,ii)), -integral_span,integral_span,-integral_span,integral_span,-integral_span,integral_span,'AbsTol',1e-7,'RelTol',1e-4);
        Prr_part(ii) = integral3(@(vr,vp,vz)vr.^2.*fe3(vr,vp,vz,x(ii),y(:,ii)), -integral_span,integral_span,-integral_span,integral_span,-integral_span,integral_span,'AbsTol',1e-7,'RelTol',1e-4);
        Pzz_part(ii) = integral3(@(vr,vp,vz)vz.^2.*fe3(vr,vp,vz,x(ii),y(:,ii)), -integral_span,integral_span,-integral_span,integral_span,-integral_span,integral_span,'AbsTol',1e-7,'RelTol',1e-4);
    end
    intparts = struct('n',neb_part,'u',ueb_part,'uz', uez_part,'pp',Ppp_part,'rr',Prr_part,'zz',Pzz_part);
end

function [linedata, intparts] = calc_linemodel(x,y,yp,params,varargin)
    opts = struct(varargin{:});
    if isfield(opts, 'intparts')
        useint = true;
    elseif isfield(opts, 'fe3')
        useint = false;
    else
        error('cannot calc integral3')
    end
    linedata = struct();
    linedata.x = x;
    linedata.Bz = (params.epsilon + y(2,:))*params.Bst;
    linedata.Bphi = -yp(3,:)*params.Bst;
    linedata.E = -yp(5,:)*params.E0;
    % integrals
    % if ~useint
    %     warning('No given integrals. Recalculating may take some time.')
    %     opts.intbparts = calc_integral_parts()
    % end
    if useint
        intparts = opts.intparts;
        ne_part = opts.intparts.n;
        ue_part = opts.intparts.u;
        uz_part = opts.intparts.uz;
        Ppp_part = opts.intparts.pp;
        Prr_part = opts.intparts.rr;
        Pzz_part = opts.intparts.zz;
    else
        % fe3 = opts.fe3;
        error('Use Calc_integral instead');
    end
    % Ni: note that params.n0i may change during iteration.
    linedata.Ni = params.n0i*(exp(-y(5,:)/params.tau));
    % linedata.Ner = - params.alpha_rot*(1/pi)*exp(y(5,:) - y(1,:)).*ner_part;
    % linedata.Nez = params.alpha_z*(1/pi/params.tauz)*exp((y(5,:)-params.ga*y(3,:))/params.tauz).*nez_part;
    % linedata.Neb = (1+params.alpha_rot-params.alpha_z)*(1/pi/params.taub)*exp(y(5,:)/params.taub).*neb_part;
    linedata.Ne = params.n0e*ne_part;
    linedata.Ni = linedata.Ni * (linedata.Ne(1)/linedata.Ni(1));
    fprintf('Ne0/Ni0: %.2f\n', linedata.Ne(1)/linedata.Ni(1));
    % Ve
    % linedata.nvphi = ((1+params.alpha_rot-params.alpha_z)*(1/pi/params.taub)*ueb_part.*exp(y(5,:)/params.taub)*params.v0e...
    %     -params.alpha_rot*(1/pi)*uer_part.*exp(y(5,:) - y(1,:))*params.v0e...
    %     +params.alpha_z*(1/pi/params.tauz)*uez_part.*exp((y(5,:)-params.ga*y(3,:))/params.tauz)*params.v0e )*(params.n0e*1e6); % nv = j/-e, /m2 s
    linedata.nvphi = ue_part *params.n0e*1e6 * params.v0e;
    linedata.vphi = linedata.nvphi./(linedata.Ne*1e6)/1e3;
    % linedata.vz = (params.ga*params.Lst*1e-2*params.Omegae)*...
    %     linedata.Nez*params.n0e./linedata.Ne/1e3;
    linedata.vz = uz_part *params.n0e* params.v0e./linedata.Ne/1e3;
    % Pe
    linedata.Prr = Prr_part* (2*params.te*params.e * params.n0e*1e6);
    linedata.Ppp = Ppp_part* (2*params.te*params.e * params.n0e*1e6) - linedata.nvphi.*linedata.vphi*1e3*params.me; % 2 means mv0e^2=2Te
    linedata.Pzz = Pzz_part* (2*params.te*params.e * params.n0e*1e6) - (linedata.vz*1e3).^2.*(linedata.Ne * 1e6)*params.me;
    % grad Pe
    linedata.gradPe = (dif2(linedata.Prr, x) + (linedata.Prr - linedata.Ppp)./x)/(params.Lst*1e-2);
    linedata.gradPe(1) = interp1(x(2:end), linedata.gradPe(2:end), 0, 'linear', 'extrap'); % handle 0/0 at x=0
    % Trr, Tzz
    linedata.Trr = linedata.Prr./(linedata.Ne*1e6)/params.e;
    linedata.Tpp = linedata.Ppp./(linedata.Ne*1e6)/params.e;
    linedata.Tzz = linedata.Pzz./(linedata.Ne*1e6)/params.e;
    % Generalized Ohm's Law
    Btotal = vecnorm([linedata.Bz; linedata.Bphi], 2, 1);
    linedata.vExBphi = -(linedata.E*1e-3).*(linedata.Bz*1e-9)./(Btotal*1e-9).^2*1e-3;
    linedata.vExBz = (linedata.E*1e-3).*(linedata.Bphi*1e-9)./(Btotal*1e-9).^2*1e-3;
    linedata.vgradPphi = -linedata.gradPe.*(linedata.Bz*1e-9)./(Btotal*1e-9).^2./(linedata.Ne*1e6)/params.e*1e-3;
    linedata.vgradPz = linedata.gradPe.*(linedata.Bphi*1e-9)./(Btotal*1e-9).^2./(linedata.Ne*1e6)/params.e*1e-3;
    linedata.vExB = (linedata.E*1e-3)./(Btotal*1e-9)*1e-3;
    linedata.vgradP = linedata.gradPe./(Btotal*1e-9)./(linedata.Ne*1e6)/params.e*1e-3;
    linedata.vinert =  params.me/params.e*(linedata.vphi*1e3).^2./(Btotal*1e-9)./(linedata.x*params.Lst*1e-2)*1e-3;

    linedata.Econv = ((linedata.vz.*linedata.Bphi) - (linedata.vphi.*linedata.Bz))*(1e3*1e-9*1e3);
    linedata.Egradp =  - linedata.gradPe./(linedata.Ne*1e6)/params.e*1e3;

    % interps
    psi_interp = griddedInterpolant(x, y(1,:));
    a_interp = griddedInterpolant(x, y(3,:));
    phi_interp = griddedInterpolant(x, y(5,:));

    yfunc = cell(1,5);
    % extrapolant: 1/2Bzr^2
    yfunc{1} = @(xx)(xx<=x(end)).*psi_interp(xx) +...
        (xx>x(end)).*(y(1,end) + y(2,end)/2*(xx.^2 - x(end)^2));
    % extrapolant: Bphi = mu0I/(2pi r)
    yfunc{3} = @(xx)y3func(xx,x,y,a_interp);
    phi_interp.ExtrapolationMethod = 'nearest';
    yfunc{5} = phi_interp;

    Bzfunc = griddedInterpolant(x, linedata.Bz);
    Bzfunc.ExtrapolationMethod = 'nearest';
    Bphi_func = griddedInterpolant(x, linedata.Bphi);
    Bphifunc = @(xx)bpfunc(xx,x,Bphi_func);
%     E_interp = griddedInterpolant(x, linedata.E);
%     Erfunc = @(xx)(xx<=x(end)).*E_interp(xx);
%     gradBz_interp = griddedInterpolant([0 x(1:end-1) + diff(x)/2], [0 diff(linedata.Bz)./diff(x)]);
%     gradBzfunc = @(xx)(xx<=x(end)).*gradBz_interp(xx);
    % fe
    linedata.fe = @(vr,vp,vz,xx)(params.n0e*1e6)*(...
        (1+params.alpha_rot-params.alpha_z)*(1 - params.bebg*params.Bst/params.bz)*(1/pi/params.taub)^1.5*exp(-1/params.taub*((vr.^2+vp.^2+vz.^2) - params.bebg*mu_full_func(vr,vp,vz,xx,Bzfunc,Bphifunc,params))).*exp(yfunc{5}(xx)/params.taub)...
        + params.alpha_z*(1 - params.bez*params.Bst/params.bz)*(1/pi/params.tauz)^1.5*exp(-1/params.tauz*(vr.^2+vp.^2+(vz - params.ga*params.Lst*1e-2*params.Omegae/params.v0e).^2 - params.bez*mu_full_func(vr,vp,vz,xx,Bzfunc,Bphifunc,params))).*exp((yfunc{5}(xx) - params.ga*yfunc{3}(xx))/params.tauz)...
        - params.alpha_rot*(1 - params.berot*params.Bst/params.bz)*(1/pi)^1.5*exp(-(vr.^2+(vp - params.Omegae*xx*params.Lst*1e-2/params.v0e).^2 + vz.^2 - params.berot*mu_full_func(vr,vp,vz,xx,Bzfunc,Bphifunc,params))).*exp(yfunc{5}(xx) - yfunc{1}(xx)) );
end


function a = y3func(xx,x,y,a_interp)
    a = zeros(size(xx));
    inbound = xx <= x(end);
    a(inbound) = a_interp(xx(inbound));
    a(~inbound) = y(3,end) + y(4,end)*log(xx(~inbound)/x(end));
    % for ii = 1:numel(xx)
    %     if xx(ii) <= x(end)
    %         a(ii) = a_interp(xx(ii));
    %     else
    %         a(ii) = y(3,end) + y(4,end)*log(xx(ii)/x(end));
    %     end
    % end
end

function bphi = bpfunc(xx,x,bp_interp)
    bphi = zeros(size(xx));
    inbound = xx <= x(end);
    bphi(inbound) = bp_interp(xx(inbound));
    bphi(~inbound) = bp_interp(x(end))*x(end)./xx(~inbound);
%     for ii = 1:numel(xx)
%         if xx(ii) <= x(end)
%             bphi(ii) = bp_interp(xx(ii));
%         else
%             bphi(ii) = bp_interp(x(end))*x(end)/xx(ii);
%         end
%     end
end