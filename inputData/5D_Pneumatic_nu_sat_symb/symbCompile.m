function problem = symbCompile(par,cfg)
%SYMBCOMPILE Derive and compile the symbolic problem used by the solver.
%   problem = symbCompile(par,cfg)   (called by buildProblem.m when symbFlag = 1)
%
% Same roles as in the reference _sym examples: symbSetting.m defines the model
% and the Hamiltonian; this function substitutes the weights Q_y, R_nu and
% compiles f/g/HamF/HamG, the substituted Hamiltonian, the stationary rate nu*
% and the canonical field into vectorised function handles. Runtime handles
% close over their own numeric expressions, so the relaxation sweep can build
% one problem per lambda_h without file-name clashes.
%
% All generated files and the derivation are exported for inspection to
%   result/symbolic_generated/lambda_<lambda_h>_Rnu_<Rnu>/
% The files f.m, g.m, HamF.m, HamG.m kept in this folder were generated in the
% same way for the nominal parameters (lambda_h = 0.02, Rnu = 0.02).
%
% (Originally systemData.m of the Revision-5 code that produced the paper results.)

s=symbSetting(par,cfg);
Qy=cfg.Qdiag(:).*cfg.scale(:).^2;
weights=[s.qw;s.R_sym]; numbers=sym([Qy;cfg.Rnu],'f');
numeric=@(expr) subs(expr,weights,numbers);
f=s.f_sym; g=s.g_sym;
pre=numeric(s.preH_sym); dpre=numeric(s.dpreH_dnu);
nustar=numeric(s.nu_star_sym); H=numeric(s.Hamilnian_sym);
canonical=numeric(s.Ham_vecField);
linear=numeric(s.linearMat_HamVec);
hf=numeric(s.HamF_sym); hg=numeric(s.HamG_sym);

problem=struct('n',s.dim_x,'nInput',s.dim_u, ...
    'zStar',double(s.zStar),'uAppliedStar',double(s.uAppliedStar), ...
    'A',double(s.A),'B',double(s.B), ...
    'canonicalLinear',double(linear), ...
    'f',compile_vector(f,{s.x}),'g',compile_vector(g,{s.x}), ...
    'preHamiltonian',compile_vector(pre,{s.x,s.p,s.u}), ...
    'dPreDnu',compile_vector(dpre,{s.x,s.p,s.u}), ...
    'optimalRate',compile_vector(nustar,{s.x,s.p}), ...
    'hamiltonian',compile_vector(H,{s.x,s.p}), ...
    'canonical',compile_vector(canonical,{s.x,s.p}), ...
    'HamF',compile_vector(hf,{s.x,s.p}), ...
    'HamG',compile_vector(hg,{s.x,s.p}), ...
    'driftJacobian',matlabFunction(s.Ajac,'Vars',{s.x}), ...
    'saturation',matlabFunction(s.sat_sym,'Vars',{s.uraw}), ...
    'saturationDerivative',matlabFunction(s.dsat_sym,'Vars',{s.uraw}), ...
    'inverseSaturation',matlabFunction(s.inverse_sat_sym,'Vars',{s.uapplied}), ...
    'inputCurvature',double(numeric(s.d2preH_dnu2)), ...
    'sourceFolder',fileparts(mfilename('fullpath')), ...
    'status','Symbolic saturation/preH/derivative/substitution/canonical pipeline active');
assert(problem.inputCurvature>0,'The rate minimization must be strictly convex.');

% Export inspectable reference-style files and the derivation. These are
% generated outputs, not a separate hand-written numerical model.
label=strrep(sprintf('lambda_%g_Rnu_%g',par.lambda_h,cfg.Rnu),'.','p');
folder=fullfile(cfg.resultFolder,'symbolic_generated',label);
if ~exist(folder,'dir'), mkdir(folder); end
entries={ ...
    'f',f,{s.x}; 'g',g,{s.x}; ...
    'pre_hamiltonian',pre,{s.x,s.p,s.u}; ...
    'dpre_dnu',dpre,{s.x,s.p,s.u}; ...
    'nu_star',nustar,{s.x,s.p}; ...
    'hamiltonian',H,{s.x,s.p}; ...
    'canonicalField',canonical,{s.x,s.p}; ...
    'HamF',hf,{s.x,s.p}; 'HamG',hg,{s.x,s.p}; ...
    'driftJacobian',s.Ajac,{s.x}; ...
    'sat_u',s.sat_sym,{s.uraw}; ...
    'dsat_u',s.dsat_sym,{s.uraw}; ...
    'sat_inverse',s.inverse_sat_sym,{s.uapplied}};
for j=1:size(entries,1)
    matlabFunction(entries{j,2},'File',fullfile(folder,[entries{j,1} '.m']), ...
        'Vars',entries{j,3},'Optimize',false);
end
save(fullfile(folder,'symbolic_derivation.mat'),'s','par','cfg','-v7');
fid=fopen(fullfile(folder,'hamiltonian_derivation.txt'),'w');
if fid<0, error('HJ:SymbolicReport','Cannot write symbolic derivation.'); end
cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'FIVE-STATE SATURATED-COMMAND SYMBOLIC PROBLEM\n');
fprintf(fid,'Raw state u; applied command sigma(u); optimized input nu=du/dt.\n');
fprintf(fid,'sigma(u) = %s\n\n',char(s.sat_sym));
fprintf(fid,'PRE-HAMILTONIAN = %s\n\n',char(pre));
fprintf(fid,'DERIVATIVE d(preH)/dnu = %s\n\n',char(dpre));
fprintf(fid,'STATIONARY RATE nu_star = %s\n\n',char(nustar));
fprintf(fid,'SUBSTITUTED HAMILTONIAN = %s\n\n',char(H));
fprintf(fid,'CANONICAL VECTOR FIELD = %s\n',char(canonical));
problem.generatedFolder=folder;
problem.generated=true;
problem.stationarityIdentityVerified=true;
end

function fun = compile_vector(expr,vars)
% Compile each scalar component so constant zero rows broadcast correctly
% for the 5-by-N and 10-by-N batches used by the numerical manifold solver.
parts=cell(numel(expr),1);
for j=1:numel(expr), parts{j}=matlabFunction(expr(j),'Vars',vars); end
fun=@(varargin) evaluate_components(parts,varargin{:});
end

function values = evaluate_components(parts,varargin)
np=1;
for k=1:numel(varargin), np=max(np,size(varargin{k},2)); end
values=zeros(numel(parts),np);
for j=1:numel(parts)
    value=parts{j}(varargin{:});
    if isscalar(value)
        values(j,:)=value;
    else
        assert(numel(value)==np,'Generated symbolic function has inconsistent batch dimensions.');
        values(j,:)=reshape(value,1,np);
    end
end
end
