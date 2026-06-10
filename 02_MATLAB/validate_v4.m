function validate_v4()
%VALIDATE_V4  Cross-check bioprinting_algorithm_v4 against the legacy solvers.
% Uses C15 / 21G at Vp = 0.01 mm/s, horizontal (hydrostatic = 0) so the
% legacy totals are directly comparable.

here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(here, 'archive', 'scripts'));   % deprecated PL/Cross solvers (audit-retained)
fprintf('==================== VALIDATE v4 ====================\n');

% ---- inputs (C15 Ramp1) ----
name = 'C15';
K = 92.5224; n = 0.397475;
eta0 = 578.26; etaInf = 0.0162966; lambda = 3.71132; mC = 0.810869;
rho = 898; Rrec = 55;
Rs = 14.3e-3/2; Rn = 0.515e-3/2; Ln = 31.75e-3; Ls = 90e-3;
Vp = 0.01e-3;   % m/s

sample = struct('name',name,'K_PL',K,'n_PL',n,'eta0',eta0,'etaInf',etaInf, ...
    'lambda',lambda,'m_Cross',mC,'rho',rho,'Rrec_pct',Rrec);
geom = struct('Rs',Rs,'R_n',Rn,'L_n',Ln,'Ls',Ls,'label','21G','h_factor',0.7);

tol_dir = tempname; mkdir(tol_dir);
res = bioprinting_algorithm_v4(sample, geom, Vp, ...
    'OutputFolder',tol_dir,'Orientation','horizontal','IncludeHydrostatic',true, ...
    'SaveData',true,'SaveFigures',false,'SaveCSV',true,'NumPoints',200,'PlotResults',false);

% ---- legacy PL ----
pl = bioprinting_algorithm_3(Rs, Rn, Ln, Ls, Vp, K, n, rho, ...
    'Name',name,'PlotResults',false,'SaveData',false,'SaveFigures',false, ...
    'Orientation','horizontal','IncludeHydrostatic',true,'OutputFolder',tol_dir,'NumPoints',200);

% ---- legacy Cross (straight: Rn_in = Rn_out) ----
cr = bioprinting_algorithm_cross_v2(Rs, Rn, Rn, Ln, Ls, Vp, eta0, etaInf, lambda, mC, rho, ...
    'Name',name,'PlotResults',false,'SaveData',false,'SaveFigures',false, ...
    'Orientation','horizontal','IncludeHydrostatic',true,'OutputFolder',tol_dir,'NumPoints',200);

% ---- comparisons ----
fprintf('\n--- POWER-LAW (v4 vs algorithm_3) ---\n');
chk('Q',                res.PL.Q(1),            pl.flow.Q);
chk('needle dP',        res.PL.dP_needle(1),    pl.pressure.drop_needle);
chk('syringe dP',       res.PL.dP_syr(1),       pl.pressure.drop_syringe);
chk('total dP',         res.PL.dP_total(1),     pl.pressure.drop_total_system);
chk('needle tau_w',     res.PL.tau_w_needle(1), pl.shear.tau_wall_needle);
chk('needle gamma_true',res.PL.gamma_w_needle(1),pl.shear.gamma_dot_true_needle);
chk('needle Re',        res.PL.Re_needle(1),    pl.reynolds.Re_needle);
chk('needle u_max',     res.PL.u_max_needle(1), pl.flow.u_max_needle);

fprintf('\n--- CROSS (v4 vs cross_v2 straight) ---\n');
chk('needle dP',        res.Cross.dP_needle(1),  cr.pressure.drop_needle);
chk('syringe dP',       res.Cross.dP_syr(1),     cr.pressure.drop_syringe);
chk('needle tau_w',     res.Cross.tau_w_needle(1),cr.shear.tau_wall_needle);
chk('needle gamma_w',   res.Cross.gamma_w_needle(1),cr.shear.gamma_dot_wall_needle);
chk('needle Re',        res.Cross.Re_needle(1),  cr.reynolds.Re_needle);

fprintf('\n--- SLICER (v4 vs v3 formulas) ---\n');
beta_PL = max(0,0.30*(1-n)); w_PL = 2*Rn*(1+beta_PL); h = 0.7*2*Rn;
vprint_PL = (Rs^2*Vp)*pi /(w_PL*h);   % Q/(w*h), Q=pi Rs^2 Vp
kflow_PL = (1+beta_PL)^2*1.0*sqrt(Rrec/100);
chk('PL v_print',  res.slicer.PL.v_print_mm_s(1), vprint_PL*1e3);
chk('PL k_flow',   res.slicer.PL.k_flow,          kflow_PL);
chk('PL w_line',   res.slicer.PL.w_line_mm,       w_PL*1e3);

% ---- data file check ----
vp_tag = sprintf('Vp%smmps', strrep(sprintf('%g',Vp*1e3),'.','p'));
df = fullfile(tol_dir, vp_tag, 'C15_21G_data.txt');
if isfile(df)
    txt = fileread(df);
    has_pl  = contains(txt,'POWER-LAW MODEL');
    has_cr  = contains(txt,'CROSS MODEL');
    has_np  = contains(txt,'Radial Profiles - Needle');
    has_sp  = contains(txt,'Radial Profiles - Syringe');
    fprintf('\n--- DATA FILE %s ---\n', df);
    fprintf('  has PL section: %d | has Cross section: %d | needle profile: %d | syringe profile: %d\n', ...
        has_pl, has_cr, has_np, has_sp);
else
    fprintf('\n[FAIL] data file not found: %s\n', df);
end

fprintf('\n==================== DONE ====================\n');
end

function chk(label, a, b)
    if b == 0
        rel = abs(a-b);
    else
        rel = abs(a-b)/abs(b);
    end
    if rel < 1e-3
        flag = 'OK  ';
    elseif rel < 2e-2
        flag = 'ok~ ';
    else
        flag = 'FAIL';
    end
    fprintf('  [%s] %-18s v4=%.6e legacy=%.6e  rel=%.2e\n', flag, label, a, b, rel);
end
