function s = symbSetting(par,cfg)
%SYMBSETTING Authoritative five-state plant and symbolic HJ derivation.
% z=[theta;v;pressure;hysteresis;RAW command u], optimized input nu=du/dt.
% The pressure equation contains APPLIED command sigma(u), not raw u.
% Saturation is part of the plant BEFORE constructing the pre-Hamiltonian.
% Edit the model here; the main solver consumes functions derived from here.
% Symbolic Math Toolbox is required. Cost has NO factor 1/2.

if exist('sym','file')~=2
    error('HJ:SymbolicRequired','This formulation requires MATLAB Symbolic Math Toolbox.');
end
n=5; y=sym('y',[n 1],'real'); mu=sym('mu',[n 1],'real');
nu=sym('nu','real'); R_sym=sym('Rnu','positive');
qw=sym('qweight',[n 1],'positive'); Q_sym=diag(qw);
S=sym(cfg.scale(:),'f');
J=sym(par.J,'f'); k=sym(par.k,'f'); c=sym(par.c,'f');
bp=sym(par.bp,'f'); kh=sym(par.kh,'f');
taup=sym(par.taup,'f'); ku=sym(par.ku,'f');
Ah=sym(par.Ah,'f'); beta=sym(par.beta,'f'); gamma=sym(par.gamma,'f');
lambda_h=sym(par.lambda_h,'f'); hr=sym(par.hRelaxEq,'f');
eh=sym(par.eps_abs_h,'f'); er=sym(par.eps_abs_rate,'f');
umin=sym(cfg.uMin,'f'); umax=sym(cfg.uMax,'f');
mid=(umin+umax)/2; half=(umax-umin)/2;

% Smooth bounded amplitude saturation, with 0 <= sigma'(u) <= 1.
% It is a smooth actuator model, not an exact hard clip.
uraw=sym('uraw','real'); ua=sym('uapplied','real');
sat_sym=mid+half*tanh((uraw-mid)/half);
dsat_sym=diff(sat_sym,uraw);
inverse_sat_sym=mid+half*atanh((ua-mid)/half);

thetaStar=sym(cfg.thetaStar,'f'); pStar=(k*thetaStar-kh*hr)/bp;
uAppliedStar=pStar/ku;
target=double(uAppliedStar);
if ~(target>cfg.uMin && target<cfg.uMax)
    error('HJ:Equilibrium','The target APPLIED command must lie strictly inside saturation limits.');
end
uRawStar=subs(inverse_sat_sym,ua,uAppliedStar);
zStar=[thetaStar;sym(0);pStar;hr;uRawStar];
q=S.*y; raw=uRawStar+q(5); h=hr+q(4);
uApplied=subs(sat_sym,uraw,raw);

% Equilibrium-shifted version of r=(-p+ku*sigma(u))/taup.
r=(-q(3)+ku*(uApplied-uAppliedStar))/taup;
% Rationalized smooth absolute value avoids subtraction near zero.
sh=h^2/(sqrt(h^2+eh^2)+eh);
sr=r^2/(sqrt(r^2+er^2)+er);
hd=(Ah-gamma*sh)*r-beta*h*sr-lambda_h*q(4);
Fphysical=[q(2);(-k*q(1)-c*q(2)+bp*q(3)+kh*q(4))/J;r;hd;sym(0)];
f_sym=Fphysical./S; g_sym=sym([0;0;0;0;1])./S;

% 1. PRE-HAMILTONIAN: saturation is already present in f_sym.
preH_sym=mu.'*(f_sym+g_sym*nu)+y.'*Q_sym*y+R_sym*nu^2;
% 2. DIFFERENTIATE with respect to the optimized RATE nu.
dpreH_dnu=diff(preH_sym,nu);
d2preH_dnu2=diff(dpreH_dnu,nu);
% 3. SOLVE stationarity. u is a state, so do not minimize over u.
nu_star_sym=solve(dpreH_dnu==0,nu);
if numel(nu_star_sym)~=1
    error('HJ:Stationarity','Expected one strictly convex quadratic rate minimizer.');
end
nu_star_sym=simplify(nu_star_sym);
stationarity_identity=simplify(subs(dpreH_dnu,nu,nu_star_sym));
if ~isAlways(stationarity_identity==0)
    error('HJ:Stationarity','Substituted symbolic rate does not satisfy stationarity.');
end
% 4. SUBSTITUTE into the PRE-Hamiltonian; do not insert a memorized Hbar.
Hamilnian_sym=subs(preH_sym,nu,nu_star_sym); % reference spelling retained
% 5. DIFFERENTIATE the resulting minimized Hamiltonian.
dHdy=jacobian(Hamilnian_sym,y).';
dHdmu=jacobian(Hamilnian_sym,mu).';
Ham_vecField=[dHdmu;-dHdy];
Ajac=jacobian(f_sym,y);
A=simplify(subs(Ajac,y,zeros(n,1)));
B=g_sym;
linearMat_HamVec=simplify(subs(jacobian(Ham_vecField,[y;mu]), ...
    [y;mu],zeros(2*n,1)));
nonlinear=Ham_vecField-linearMat_HamVec*[y;mu];

s=struct('dim_x',n,'dim_u',1,'x',y,'p',mu,'u',nu, ...
    'Q_sym',Q_sym,'qw',qw,'R_sym',R_sym,'scale',S,'zStar',zStar, ...
    'uAppliedStar',uAppliedStar,'uraw',uraw,'uapplied',ua, ...
    'sat_sym',sat_sym,'dsat_sym',dsat_sym,'inverse_sat_sym',inverse_sat_sym, ...
    'f_sym',f_sym,'g_sym',g_sym,'Ajac',Ajac,'A',A,'B',B, ...
    'preH_sym',preH_sym,'dpreH_dnu',dpreH_dnu, ...
    'd2preH_dnu2',d2preH_dnu2,'nu_star_sym',nu_star_sym, ...
    'stationarity_identity',stationarity_identity, ...
    'Hamilnian_sym',Hamilnian_sym,'Ham_vecField',Ham_vecField, ...
    'linearMat_HamVec',linearMat_HamVec,'HamF_sym',nonlinear(1:n), ...
    'HamG_sym',nonlinear(n+1:end));
end
