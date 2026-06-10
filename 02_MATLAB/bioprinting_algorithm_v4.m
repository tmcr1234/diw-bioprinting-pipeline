function results = bioprinting_algorithm_v4(sample, geom, Vp_vec, varargin)
%BIOPRINTING_ALGORITHM_V4  Unified superset extrusion solver (straight needle).
%
% A single algorithm that fully supersedes:
%   - bioprinting_algorithm_3      (Power-Law, straight cylinder)
%   - bioprinting_algorithm_cross_v2 (Cross, straight-needle degenerate case)
%   - bioprinting_algorithm_v3     (slicer-parameter generator)
%
% For each piston velocity in Vp_vec it computes, for BOTH the Power-Law and
% the Cross constitutive models, the complete syringe+needle flow analysis:
%   * volumetric flow rate, mean/max velocities (syringe + needle)
%   * full pressure decomposition: syringe drop, needle drop, hydrostatic
%     head, total system drop, absolute exit/needle-inlet/syringe-inlet
%   * wall shear rate / stress (syringe + needle; PL also Rabinowitsch-true)
%   * generalised Reynolds numbers + critical-Re estimate + laminar verdict
%   * radial profiles r, u(r), gamma_dot(r), tau(r) (syringe + needle)
%   * FFF-slicer-equivalent parameters (v_print, w_line, h_layer,
%     beta_swell, k_flow)
%
% It writes, per operating point, a combined human-readable _data.txt
% (Power-Law section + Cross section, legacy format) and a set of figures
% including the flow-regime estimate chart and the viscosity-across-needle-
% exit chart. Across the Vp sweep it writes a slicer-lookup CSV and a
% trend-summary figure.
%
% Geometry: STRAIGHT cylindrical needle (matches 21G/22G blunt hardware).
%
% --------------------------------------------------------------------
% INPUTS
%   sample - struct with rheology + recovery + density fields:
%       .name       sample identifier (string)
%       .K_PL       Power-Law consistency (Pa*s^n)
%       .n_PL       Power-Law flow index (-)
%       .eta0       Cross zero-shear viscosity (Pa*s)
%       .etaInf     Cross infinite-shear viscosity (Pa*s)
%       .lambda     Cross time constant (s)
%       .m_Cross    Cross exponent (-)
%       .rho        density (kg/m^3)
%       .Rrec_pct   structural recovery (%) at deposition shear rate
%
%   geom - struct:
%       .Rs         syringe inner radius (m)
%       .R_n        needle inner radius (m)
%       .L_n        needle length (m)
%       .Ls         syringe length (m)
%       .label      label, e.g. '21G' or '22G'
%       .h_factor   layer height as fraction of needle ID (default 0.7)
%
%   Vp_vec - vector of piston velocities (m/s)
%
% NAME-VALUE OPTIONAL
%   'OutputFolder'       (default 'output_v4')
%   'Orientation'        'horizontal'|'upward'|'downward' (default 'horizontal')
%   'IncludeHydrostatic' (default true)
%   'SaveData'           (default true)   per-Vp combined _data.txt
%   'SaveFigures'        (default true)   per-Vp + summary PNGs
%   'SaveCSV'            (default true)   across-Vp slicer-lookup CSV
%   'PlotResults'        (default true)   master switch for any figure
%   'NumPoints'          (default 200)    radial discretisation
%   'f_slip'             (default 1.0)    syringe-piston slip factor
%
% OUTPUT
%   results - struct with .PL and .Cross sub-structs (each field is a column
%       vector over Vp), .profiles{iVp} (per-Vp radial data for both models),
%       .slicer (PL/Cross slicer params), plus .sample, .geom, .Vp_mm_s.
%
% Author: T.M.C. Rodrigues (PEMM/COPPE/UFRJ) -- 2026-06-01
% --------------------------------------------------------------------

%% -------------------- INPUT PARSING --------------------
p = inputParser;
addRequired(p, 'sample', @isstruct);
addRequired(p, 'geom',   @isstruct);
addRequired(p, 'Vp_vec', @(x) isnumeric(x) && all(x > 0));
addParameter(p, 'OutputFolder',       'output_v4', @(x) ischar(x)||isstring(x));
addParameter(p, 'Orientation',        'horizontal', @(x) any(validatestring(x,{'horizontal','upward','downward'})));
addParameter(p, 'IncludeHydrostatic', true,  @islogical);
addParameter(p, 'SaveData',           true,  @islogical);
addParameter(p, 'SaveFigures',        true,  @islogical);
addParameter(p, 'SaveCSV',            true,  @islogical);
addParameter(p, 'PlotResults',        true,  @islogical);
addParameter(p, 'NumPoints',          200,   @(x) isnumeric(x) && x >= 20);
addParameter(p, 'f_slip',             1.0,   @(x) isnumeric(x) && x > 0);
parse(p, sample, geom, Vp_vec, varargin{:});

output_folder = char(p.Results.OutputFolder);
orientation   = validatestring(p.Results.Orientation, {'horizontal','upward','downward'});
inc_hydro     = p.Results.IncludeHydrostatic;
save_data     = p.Results.SaveData;
save_figs     = p.Results.SaveFigures;
save_csv      = p.Results.SaveCSV;
plot_results  = p.Results.PlotResults;
num_points    = p.Results.NumPoints;
f_slip        = p.Results.f_slip;

if ~exist(output_folder, 'dir'), mkdir(output_folder); end
safe_name = regexprep(char(sample.name), '[<>:"/\\|?*\.\s'']', '_');

%% -------------------- GEOMETRY / CONSTANTS --------------------
if ~isfield(geom, 'h_factor'), geom.h_factor = 0.7; end
if ~isfield(geom, 'Ls'),       geom.Ls = 90e-3;      end
Rs  = geom.Rs;   Rn = geom.R_n;   Ln = geom.L_n;   Ls = geom.Ls;
rho = sample.rho;
g   = 9.81;
P_atm = 101325;
h_layer = geom.h_factor * 2 * Rn;            % m

switch orientation
    case 'horizontal', orient_factor = 0;
    case 'upward',     orient_factor = +1;   % flow against gravity
    case 'downward',   orient_factor = -1;   % gravity assists flow
end
if inc_hydro
    hydro = orient_factor * rho * g * (Ls + Ln);
else
    hydro = 0;
end

% Rheology shorthands
K = sample.K_PL;   n = sample.n_PL;
eta0 = sample.eta0; etaInf = sample.etaInf; lambda = sample.lambda; mC = sample.m_Cross;
Rrec_pct = sample.Rrec_pct;

n_v = numel(Vp_vec);
A_syringe = pi * Rs^2;

% Pre-allocate per-Vp result vectors (PL + Cross)
[PL.Q, PL.dP_syr, PL.dP_needle, PL.dP_total, PL.absP_syr_in, PL.absP_needle_in, ...
 PL.tau_w_needle, PL.gamma_w_needle, PL.tau_w_syr, PL.gamma_w_syr, ...
 PL.Re_needle, PL.Re_syr, PL.u_avg_needle, PL.u_max_needle] = deal(zeros(n_v,1));
[CR.Q, CR.dP_syr, CR.dP_needle, CR.dP_total, CR.absP_syr_in, CR.absP_needle_in, ...
 CR.tau_w_needle, CR.gamma_w_needle, CR.tau_w_syr, CR.gamma_w_syr, ...
 CR.Re_needle, CR.Re_syr, CR.u_avg_needle, CR.u_max_needle] = deal(zeros(n_v,1));
PL.Re_crit = 2100 * n^0.75;
profiles = cell(n_v,1);

%% ==================================================================
%% MAIN Vp LOOP
%% ==================================================================
for iv = 1:n_v
    Vp = Vp_vec(iv);
    Q  = A_syringe * Vp;                      % m^3/s (same for both models)

    % ---------- POWER-LAW (analytical, straight cylinder) ----------
    pl_needle = pl_section(Q, K, n, Rn, Ln, rho, num_points);
    pl_syr    = pl_section(Q, K, n, Rs, Ls, rho, num_points);

    pl_dP_total = pl_needle.drop + pl_syr.drop + hydro;
    pl_absP_needle_in = P_atm + pl_needle.drop;
    pl_absP_syr_in    = pl_absP_needle_in + pl_syr.drop;

    PL.Q(iv)            = Q;
    PL.dP_syr(iv)       = pl_syr.drop;
    PL.dP_needle(iv)    = pl_needle.drop;
    PL.dP_total(iv)     = pl_dP_total;
    PL.absP_needle_in(iv)= pl_absP_needle_in;
    PL.absP_syr_in(iv)  = pl_absP_syr_in;
    PL.tau_w_needle(iv) = pl_needle.tau_w;
    PL.gamma_w_needle(iv)= pl_needle.gamma_true;
    PL.tau_w_syr(iv)    = pl_syr.tau_w;
    PL.gamma_w_syr(iv)  = pl_syr.gamma_true;
    PL.Re_needle(iv)    = pl_needle.Re;
    PL.Re_syr(iv)       = pl_syr.Re;
    PL.u_avg_needle(iv) = pl_needle.u_avg;
    PL.u_max_needle(iv) = pl_needle.u_max;

    % ---------- CROSS (dpdz-from-Q, straight cylinder) ----------
    cr_needle = cross_section(Q, eta0, etaInf, lambda, mC, Rn, Ln, rho, num_points);
    cr_syr    = cross_section(Q, eta0, etaInf, lambda, mC, Rs, Ls, rho, num_points);

    cr_dP_total = cr_needle.drop + cr_syr.drop + hydro;
    cr_absP_needle_in = P_atm + cr_needle.drop;
    cr_absP_syr_in    = cr_absP_needle_in + cr_syr.drop;

    CR.Q(iv)            = Q;
    CR.dP_syr(iv)       = cr_syr.drop;
    CR.dP_needle(iv)    = cr_needle.drop;
    CR.dP_total(iv)     = cr_dP_total;
    CR.absP_needle_in(iv)= cr_absP_needle_in;
    CR.absP_syr_in(iv)  = cr_absP_syr_in;
    CR.tau_w_needle(iv) = cr_needle.tau_w;
    CR.gamma_w_needle(iv)= cr_needle.gamma_w;
    CR.tau_w_syr(iv)    = cr_syr.tau_w;
    CR.gamma_w_syr(iv)  = cr_syr.gamma_w;
    CR.Re_needle(iv)    = cr_needle.Re;
    CR.Re_syr(iv)       = cr_syr.Re;
    CR.u_avg_needle(iv) = cr_needle.u_avg;
    CR.u_max_needle(iv) = cr_needle.u_max;

    % ---------- store per-Vp radial profiles (both models) ----------
    prof.Vp = Vp; prof.Q = Q;
    prof.PL = struct('needle',pl_needle,'syr',pl_syr, ...
        'dP_total',pl_dP_total,'absP_needle_in',pl_absP_needle_in,'absP_syr_in',pl_absP_syr_in);
    prof.CR = struct('needle',cr_needle,'syr',cr_syr, ...
        'dP_total',cr_dP_total,'absP_needle_in',cr_absP_needle_in,'absP_syr_in',cr_absP_syr_in);
    profiles{iv} = prof;

    % ---------- per-Vp DATA FILE ----------
    if save_data
        write_combined_data(output_folder, safe_name, geom, sample, ...
            Vp, orientation, inc_hydro, hydro, prof);
    end

    % ---------- per-Vp FIGURES ----------
    if plot_results && save_figs
        make_perVp_figures(output_folder, safe_name, geom, sample, Vp, prof, PL.Re_crit);
    end
end

%% ==================================================================
%% SLICER LAYER (across Vp)
%% ==================================================================
beta_PL = max(0, 0.30*(1 - n));
beta_CR = max(0, 0.30*(1 - mC));
w_line_PL = 2*Rn*(1 + beta_PL);
w_line_CR = 2*Rn*(1 + beta_CR);
v_print_PL = PL.Q ./ (w_line_PL * h_layer);   % m/s
v_print_CR = CR.Q ./ (w_line_CR * h_layer);
k_flow_PL = (1+beta_PL)^2 * f_slip * sqrt(Rrec_pct/100);
k_flow_CR = (1+beta_CR)^2 * f_slip * sqrt(Rrec_pct/100);

slicer.PL = struct('beta',beta_PL,'w_line_mm',w_line_PL*1e3,'h_layer_mm',h_layer*1e3, ...
    'v_print_mm_s',v_print_PL*1e3,'k_flow',k_flow_PL);
slicer.CR = struct('beta',beta_CR,'w_line_mm',w_line_CR*1e3,'h_layer_mm',h_layer*1e3, ...
    'v_print_mm_s',v_print_CR*1e3,'k_flow',k_flow_CR);

if save_csv
    write_slicer_csv(output_folder, safe_name, geom, Vp_vec, PL, CR, ...
        v_print_PL, w_line_PL, beta_PL, k_flow_PL, ...
        v_print_CR, w_line_CR, beta_CR, k_flow_CR, h_layer);
end

if plot_results && save_figs && n_v > 1
    make_summary_figure(output_folder, safe_name, geom, Vp_vec, PL, CR, ...
        v_print_PL, v_print_CR, k_flow_PL, k_flow_CR, Rrec_pct);
end

%% -------------------- ASSEMBLE RESULTS --------------------
results.sample  = sample;
results.geom    = geom;
results.Vp_mm_s = Vp_vec(:)*1e3;
results.PL      = PL;
results.Cross   = CR;
results.profiles= profiles;
results.slicer  = slicer;
results.orientation = orientation;
results.hydrostatic_Pa = hydro;

end  % ===================== END MAIN =====================


%% ====================================================================
%% LOCAL FUNCTIONS
%% ====================================================================

function s = pl_section(Q, K, n, R, L, rho, npts)
% Power-Law fully-developed pipe flow (analytical) for one constant-radius tube.
    abs_dpdz = 2*K * ((Q*(3*n+1)) / (pi*n*R^((3*n+1)/n)))^n;
    s.abs_dpdz   = abs_dpdz;
    s.drop       = abs_dpdz * L;
    s.gamma_app  = 4*Q / (pi*R^3);
    s.gamma_true = s.gamma_app * (3*n+1)/(4*n);
    s.tau_w      = abs_dpdz * R / 2;
    s.u_avg      = Q / (pi*R^2);

    r = linspace(0, R, npts);
    u = (n/(n+1)) * (abs_dpdz/(2*K))^(1/n) .* (R^((n+1)/n) - r.^((n+1)/n));
    gdot = (abs_dpdz/(2*K))^(1/n) .* r.^(1/n);   gdot(1) = 0;
    tau  = K .* gdot.^n;
    s.r = r; s.u = u; s.gamma_dot = gdot; s.tau = tau;
    s.u_max = u(1);

    D = 2*R;
    s.Re = (rho * s.u_avg^(2-n) * D^n) / (K * 8^(n-1) * ((3*n+1)/(4*n))^n);
end

function s = cross_section(Q, eta0, etaInf, lambda, m, R, L, rho, npts)
% Cross model fully-developed pipe flow for one constant-radius tube.
% dp/dz is solved from the target Q; radial profiles follow from dp/dz.
    abs_dpdz = solve_dpdz_cross(Q, R, eta0, etaInf, lambda, m, min(npts,100));
    s.abs_dpdz = abs_dpdz;
    s.drop     = abs_dpdz * L;

    [r, u, gdot, tau] = calc_radial_profiles_cross(abs_dpdz, R, eta0, etaInf, lambda, m, npts);
    s.r = r; s.u = u; s.gamma_dot = gdot; s.tau = tau;
    s.gamma_w = max(gdot);
    s.tau_w   = max(tau);
    s.u_avg   = Q / (pi*R^2);
    s.u_max   = u(1);

    eta_wall = s.tau_w / max(s.gamma_w, 1e-12);
    s.Re = (rho * s.u_avg * (2*R)) / eta_wall;
end

function dpdz = solve_dpdz_cross(Q_target, R, eta0, etaInf, lambda, m, N)
% Bracketed fzero: dpdz_lo (eta~etaInf, flows too much) ... dpdz_hi (eta~eta0).
% Robustness guard (2026-06-01): in a wide, low-shear conduit (e.g. the
% syringe barrel) a strongly shear-thinning ink with very small etaInf stays
% at its eta0 plateau, so the root sits on/beyond the Newtonian-eta0 estimate
% dpdz_hi and the naive bracket can fail to change sign. Widen the bracket
% until a sign change is found before calling fzero.
    f = @(dp) calc_Q_discrete(dp, R, eta0, etaInf, lambda, m, N) - Q_target;
    dpdz_lo = 0.5 * (8 * etaInf * Q_target) / (pi * R^4);
    dpdz_hi =       (8 * eta0   * Q_target) / (pi * R^4);
    if dpdz_lo <= 0, dpdz_lo = dpdz_hi * 1e-9; end
    flo = f(dpdz_lo);  fhi = f(dpdz_hi);
    nexp = 0;
    while flo*fhi > 0 && nexp < 80
        dpdz_hi = dpdz_hi * 1.5;   fhi = f(dpdz_hi);   nexp = nexp + 1;
    end
    if flo*fhi > 0
        error('solve_dpdz_cross:noBracket', ...
              'Could not bracket dp/dz (Q=%.3e, R=%.3e).', Q_target, R);
    end
    opts = optimset('Display','off','TolX',1e-4);
    dpdz = fzero(f, [dpdz_lo, dpdz_hi], opts);
end

function Q = calc_Q_discrete(dpdz, R, eta0, etaInf, lambda, m, N)
    r = linspace(0, R, N);
    tau_r = (dpdz * r) / 2;
    g_r = zeros(size(r));
    for i = 2:N
        g_r(i) = solve_gamma_cross(tau_r(i), eta0, etaInf, lambda, m);
    end
    Q = pi * trapz(r, r.^2 .* g_r);
end

function g = solve_gamma_cross(tau, eta0, etaInf, lambda, m)
% Invert the Cross stress-shear relation tau = eta(g)*g for g >= 0.
% The relation is monotone, so the only fragile endpoint is the low-shear
% side, where eta -> eta0 (Newtonian plateau). Guard it with the analytic
% eta0 limit so the bracket always changes sign.
    if tau <= 0, g = 0; return; end
    obj = @(gg) (etaInf + (eta0 - etaInf)./(1 + (lambda*abs(gg)).^m)).*gg - tau;
    g_lo = 1e-12;
    g_hi = tau / etaInf;             % high-shear (etaInf) limit, obj(g_hi) >= 0
    if obj(g_lo) >= 0                % tau below the eta0-plateau resolution
        g = tau / eta0;             % exact low-shear Newtonian-eta0 limit
    else
        g = fzero(obj, [g_lo, g_hi]);
    end
end

function [r, u, g_r, tau_r] = calc_radial_profiles_cross(dpdz, R, eta0, etaInf, lambda, m, N)
    r = linspace(0, R, N);
    tau_r = (dpdz * r) / 2;
    g_r = zeros(size(r));
    for i = 2:N
        g_r(i) = solve_gamma_cross(tau_r(i), eta0, etaInf, lambda, m);
    end
    u = zeros(size(r));
    for i = 1:N-1
        u(i) = trapz(r(i:end), g_r(i:end));   % u(R)=0, du/dr = -gamma_dot
    end
end

%% -------------------- DATA WRITER --------------------
function write_combined_data(out, safe_name, geom, sample, Vp, orientation, inc_hydro, hydro, prof)
    vp_tag = sprintf('Vp%smmps', strrep(sprintf('%g', Vp*1e3),'.','p'));
    folder = fullfile(out, vp_tag);
    if ~exist(folder,'dir'), mkdir(folder); end
    fname = fullfile(folder, sprintf('%s_%s_data.txt', safe_name, geom.label));
    fid = fopen(fname, 'w');
    if fid == -1, warning('Could not open %s', fname); return; end

    P_atm = 101325;
    fprintf(fid, '=== BIOPRINTING ANALYSIS RESULTS (v4 unified, straight needle) ===\n');
    fprintf(fid, 'Sample: %s\n', sample.name);
    fprintf(fid, 'Needle: %s\n', geom.label);
    fprintf(fid, 'Orientation: %s | Hydrostatic included: %d (contribution %.6e Pa)\n\n', ...
        orientation, inc_hydro, hydro);

    fprintf(fid, 'System Parameters:\n');
    fprintf(fid, '  Syringe radius (mm): %.6f\n', geom.Rs*1000);
    fprintf(fid, '  Needle radius (mm): %.6f\n', geom.R_n*1000);
    fprintf(fid, '  Needle length (mm): %.6f\n', geom.L_n*1000);
    fprintf(fid, '  Syringe length (mm): %.6f\n', geom.Ls*1000);
    fprintf(fid, '  Piston velocity (mm/s): %.6f\n', Vp*1000);
    fprintf(fid, '  Density (kg/m^3): %.6f\n\n', sample.rho);

    % -------- POWER-LAW SECTION --------
    pl = prof.PL;
    fprintf(fid, '############### POWER-LAW MODEL ###############\n');
    fprintf(fid, 'Rheology Parameters:\n');
    fprintf(fid, '  Power-law index n: %.6f\n', sample.n_PL);
    fprintf(fid, '  Consistency index K (Pa.s^n): %.6f\n\n', sample.K_PL);
    fprintf(fid, 'Flow Analysis:\n');
    fprintf(fid, '  Volumetric flow rate Q (m^3/s): %.9e\n', prof.Q);
    fprintf(fid, '  Syringe avg velocity (m/s): %.9e\n', pl.syr.u_avg);
    fprintf(fid, '  Syringe max velocity (m/s): %.9e\n', pl.syr.u_max);
    fprintf(fid, '  Needle avg velocity (m/s): %.9e\n', pl.needle.u_avg);
    fprintf(fid, '  Needle max velocity (m/s): %.9e\n\n', pl.needle.u_max);
    fprintf(fid, 'Pressure Analysis:\n');
    fprintf(fid, '  Syringe |dp/dz| (Pa/m): %.9e\n', pl.syr.abs_dpdz);
    fprintf(fid, '  Needle  |dp/dz| (Pa/m): %.9e\n', pl.needle.abs_dpdz);
    fprintf(fid, '  Syringe pressure drop (Pa): %.9e\n', pl.syr.drop);
    fprintf(fid, '  Needle pressure drop (Pa): %.9e\n', pl.needle.drop);
    fprintf(fid, '  Total viscous pressure drop (Pa): %.9e\n', pl.syr.drop + pl.needle.drop);
    fprintf(fid, '  Hydrostatic contribution (Pa): %.9e\n', hydro);
    fprintf(fid, '  Total estimated system pressure drop (Pa): %.9e\n', pl.dP_total);
    fprintf(fid, '  Absolute exit pressure (Pa): %.9e\n', P_atm);
    fprintf(fid, '  Absolute needle inlet pressure (Pa): %.9e\n', pl.absP_needle_in);
    fprintf(fid, '  Absolute syringe inlet pressure (Pa): %.9e\n\n', pl.absP_syr_in);
    fprintf(fid, 'Shear Analysis:\n');
    fprintf(fid, '  Syringe apparent wall shear rate (s^-1): %.9e\n', pl.syr.gamma_app);
    fprintf(fid, '  Syringe true wall shear rate (s^-1): %.9e\n', pl.syr.gamma_true);
    fprintf(fid, '  Syringe wall shear stress (Pa): %.9e\n', pl.syr.tau_w);
    fprintf(fid, '  Needle apparent wall shear rate (s^-1): %.9e\n', pl.needle.gamma_app);
    fprintf(fid, '  Needle true wall shear rate (s^-1): %.9e\n', pl.needle.gamma_true);
    fprintf(fid, '  Needle wall shear stress (Pa): %.9e\n\n', pl.needle.tau_w);
    fprintf(fid, 'Reynolds Analysis:\n');
    fprintf(fid, '  Syringe generalized Reynolds number: %.9e\n', pl.syr.Re);
    fprintf(fid, '  Needle generalized Reynolds number: %.9e\n', pl.needle.Re);
    Re_crit = 2100 * sample.n_PL^0.75;
    fprintf(fid, '  Critical Reynolds estimate: %.9e\n', Re_crit);
    if pl.needle.Re < Re_crit && pl.syr.Re < Re_crit
        verdict = 'Likely laminar based on generalized Reynolds-number estimate';
    else
        verdict = 'Flow may be transitional; check assumptions and correlations carefully';
    end
    fprintf(fid, '  Comment: %s\n\n', verdict);
    write_profile_block(fid, 'Radial Profiles - Needle', 'needle', pl.needle);
    write_profile_block(fid, 'Radial Profiles - Syringe', 'syringe', pl.syr);

    % -------- CROSS SECTION --------
    cr = prof.CR;
    fprintf(fid, '\n############### CROSS MODEL (straight needle) ###############\n');
    fprintf(fid, 'Rheology Parameters (Cross Model):\n');
    fprintf(fid, '  Zero-shear viscosity eta0 (Pa.s): %.6f\n', sample.eta0);
    fprintf(fid, '  Infinite-shear viscosity etainf (Pa.s): %.6f\n', sample.etaInf);
    fprintf(fid, '  Cross time constant lambda (s): %.6f\n', sample.lambda);
    fprintf(fid, '  Cross rate constant m: %.6f\n\n', sample.m_Cross);
    fprintf(fid, 'Flow Analysis:\n');
    fprintf(fid, '  Volumetric flow rate Q (m^3/s): %.9e\n', prof.Q);
    fprintf(fid, '  Syringe avg velocity (m/s): %.9e\n', cr.syr.u_avg);
    fprintf(fid, '  Needle avg velocity (m/s): %.9e\n\n', cr.needle.u_avg);
    fprintf(fid, 'Pressure Analysis (LOWER BOUNDS -- fully-developed flow assumed):\n');
    beta_contr = (geom.Rs / geom.R_n)^2;
    fprintf(fid, '  *** Entrance losses not modelled. Contraction beta = (Rs/Rn)^2 = %.0f ***\n', beta_contr);
    fprintf(fid, '  *** For beta >> 1, Bagley correction may add 20-50%% to needle pressure. ***\n');
    fprintf(fid, '  Syringe pressure drop (Pa): %.9e\n', cr.syr.drop);
    fprintf(fid, '  Needle pressure drop (Pa): %.9e\n', cr.needle.drop);
    fprintf(fid, '  Hydrostatic contribution (Pa): %.9e\n', hydro);
    fprintf(fid, '  Total estimated system pressure drop (Pa): %.9e\n', cr.dP_total);
    fprintf(fid, '  Absolute exit pressure (Pa): %.9e\n', P_atm);
    fprintf(fid, '  Absolute needle inlet pressure (Pa): %.9e\n', cr.absP_needle_in);
    fprintf(fid, '  Absolute syringe inlet pressure (Pa): %.9e\n\n', cr.absP_syr_in);
    fprintf(fid, 'Shear Analysis (Max Wall Values):\n');
    fprintf(fid, '  Syringe wall shear rate (s^-1): %.9e\n', cr.syr.gamma_w);
    fprintf(fid, '  Syringe wall shear stress (Pa): %.9e\n', cr.syr.tau_w);
    fprintf(fid, '  Needle exit wall shear rate (s^-1): %.9e\n', cr.needle.gamma_w);
    fprintf(fid, '  Needle exit wall shear stress (Pa): %.9e\n\n', cr.needle.tau_w);
    fprintf(fid, 'Reynolds Analysis:\n');
    fprintf(fid, '  Syringe Reynolds number: %.9e\n', cr.syr.Re);
    fprintf(fid, '  Needle exit Reynolds number: %.9e\n\n', cr.needle.Re);
    write_profile_block(fid, 'Radial Profiles - Needle Exit', 'needle', cr.needle);
    write_profile_block(fid, 'Radial Profiles - Syringe', 'syringe', cr.syr);

    fclose(fid);
    fprintf('[Saved data] %s\n', fname);
end

function write_profile_block(fid, title, tag, s)
    fprintf(fid, '%s:\n', title);
    fprintf(fid, 'r_%s(m),u_%s(m/s),gamma_dot_%s(s^-1),tau_%s(Pa)\n', tag, tag, tag, tag);
    for i = 1:numel(s.r)
        fprintf(fid, '%.9e,%.9e,%.9e,%.9e\n', s.r(i), s.u(i), s.gamma_dot(i), s.tau(i));
    end
    fprintf(fid, '\n');
end

%% -------------------- SLICER CSV --------------------
function write_slicer_csv(out, safe_name, geom, Vp_vec, PL, CR, ...
        v_print_PL, w_line_PL, beta_PL, k_flow_PL, ...
        v_print_CR, w_line_CR, beta_CR, k_flow_CR, h_layer)
    nrows = numel(Vp_vec);
    T = table( ...
        repmat({safe_name}, nrows, 1), repmat({geom.label}, nrows, 1), ...
        Vp_vec(:)*1e3, PL.Q(:)*1e9, PL.u_avg_needle(:)*1e3, ...
        PL.tau_w_needle(:), PL.gamma_w_needle(:), PL.dP_total(:)/1e3, v_print_PL(:)*1e3, ...
        repmat(w_line_PL*1e3,nrows,1), repmat(beta_PL,nrows,1), repmat(k_flow_PL,nrows,1), ...
        CR.tau_w_needle(:), CR.gamma_w_needle(:), CR.dP_total(:)/1e3, v_print_CR(:)*1e3, ...
        repmat(w_line_CR*1e3,nrows,1), repmat(beta_CR,nrows,1), repmat(k_flow_CR,nrows,1), ...
        repmat(h_layer*1e3,nrows,1), ...
        'VariableNames', { ...
            'sample','needle','Vp_mm_s','Q_mm3_s','v_mean_needle_mm_s', ...
            'PL_tau_w_Pa','PL_gamma_w_invs','PL_dP_total_kPa','PL_v_print_mm_s','PL_w_line_mm','PL_beta_swell','PL_k_flow', ...
            'Cross_tau_w_Pa','Cross_gamma_w_invs','Cross_dP_total_kPa','Cross_v_print_mm_s','Cross_w_line_mm','Cross_beta_swell','Cross_k_flow', ...
            'h_layer_mm'});
    csv_name = fullfile(out, sprintf('slicer_lookup_%s_%s.csv', safe_name, geom.label));
    writetable(T, csv_name);
    fprintf('[Saved CSV] %s\n', csv_name);
end

%% -------------------- PER-Vp FIGURES --------------------
function make_perVp_figures(out, safe_name, geom, sample, Vp, prof, Re_crit)
    vp_tag = sprintf('Vp%smmps', strrep(sprintf('%g', Vp*1e3),'.','p'));
    folder = fullfile(out, vp_tag);
    if ~exist(folder,'dir'), mkdir(folder); end
    pl = prof.PL; cr = prof.CR;
    ttl = sprintf('%s | %s | Vp=%.4g mm/s', sample.name, geom.label, Vp*1e3);

    % ---- Figure 1: radial profiles (3x2), PL solid + Cross dashed ----
    f1 = figure('Visible','off','Position',[80 80 1100 1000]);
    % needle velocity
    subplot(3,2,1);
    plot(pl.needle.r*1e3, pl.needle.u*1e3,'r-','LineWidth',1.8); hold on;
    plot(cr.needle.r*1e3, cr.needle.u*1e3,'b--','LineWidth',1.8);
    xlabel('r (mm)'); ylabel('u (mm/s)'); title('Needle velocity');
    legend({'PL','Cross'},'Location','best'); grid on;
    % syringe velocity
    subplot(3,2,2);
    plot(pl.syr.r*1e3, pl.syr.u*1e3,'r-','LineWidth',1.8); hold on;
    plot(cr.syr.r*1e3, cr.syr.u*1e3,'b--','LineWidth',1.8);
    xlabel('r (mm)'); ylabel('u (mm/s)'); title('Syringe velocity'); grid on;
    % needle shear rate
    subplot(3,2,3);
    plot(pl.needle.r*1e3, pl.needle.gamma_dot,'r-','LineWidth',1.8); hold on;
    plot(cr.needle.r*1e3, cr.needle.gamma_dot,'b--','LineWidth',1.8);
    xlabel('r (mm)'); ylabel('\gamma_{dot} (1/s)'); title('Needle shear rate'); grid on;
    % syringe shear rate
    subplot(3,2,4);
    plot(pl.syr.r*1e3, pl.syr.gamma_dot,'r-','LineWidth',1.8); hold on;
    plot(cr.syr.r*1e3, cr.syr.gamma_dot,'b--','LineWidth',1.8);
    xlabel('r (mm)'); ylabel('\gamma_{dot} (1/s)'); title('Syringe shear rate'); grid on;
    % needle shear stress
    subplot(3,2,5);
    plot(pl.needle.r*1e3, pl.needle.tau,'r-','LineWidth',1.8); hold on;
    plot(cr.needle.r*1e3, cr.needle.tau,'b--','LineWidth',1.8);
    xlabel('r (mm)'); ylabel('\tau (Pa)'); title('Needle shear stress'); grid on;
    % syringe shear stress
    subplot(3,2,6);
    plot(pl.syr.r*1e3, pl.syr.tau,'r-','LineWidth',1.8); hold on;
    plot(cr.syr.r*1e3, cr.syr.tau,'b--','LineWidth',1.8);
    xlabel('r (mm)'); ylabel('\tau (Pa)'); title('Syringe shear stress'); grid on;
    sgtitle(['Radial profiles -- ' ttl],'Interpreter','none');
    exportgraphics(f1, fullfile(folder, sprintf('%s_%s_profiles.png', safe_name, geom.label)), 'Resolution', 200);
    close(f1);

    % ---- Figure 2: system (2x2) incl. flow-regime + viscosity-across-exit ----
    f2 = figure('Visible','off','Position',[100 100 1100 850]);
    % (1) axial absolute pressure profile (both models)
    subplot(2,2,1);
    z1 = linspace(0, geom.Ls, 60);
    P1_PL = pl.absP_syr_in - pl.syr.abs_dpdz*z1;
    P2_PL = pl.absP_needle_in - pl.needle.abs_dpdz*linspace(0,geom.L_n,60);
    z2 = linspace(0, geom.L_n, 60) + geom.Ls;
    P1_CR = cr.absP_syr_in - cr.syr.abs_dpdz*z1;
    P2_CR = cr.absP_needle_in - cr.needle.abs_dpdz*linspace(0,geom.L_n,60);
    plot(z1*1e3, P1_PL/1e3,'r-','LineWidth',1.8); hold on;
    plot(z2*1e3, P2_PL/1e3,'r-','LineWidth',1.8);
    plot(z1*1e3, P1_CR/1e3,'b--','LineWidth',1.8);
    plot(z2*1e3, P2_CR/1e3,'b--','LineWidth',1.8);
    xlabel('Axial position (mm)'); ylabel('Absolute pressure (kPa)');
    title('Axial pressure profile'); legend({'PL','','Cross',''},'Location','best'); grid on;
    % (2) FLOW-REGIME estimate (Reynolds) -- both models
    subplot(2,2,2);
    bar([pl.syr.Re pl.needle.Re; cr.syr.Re cr.needle.Re]');
    set(gca,'XTickLabel',{'Syringe','Needle'},'YScale','log');
    yline(Re_crit,'k--','Re_{crit}','LineWidth',1.5);
    ylabel('Generalized Reynolds number'); title('Flow-regime estimate');
    legend({'PL','Cross'},'Location','best'); grid on;
    % (3) VISCOSITY across needle exit -- both models
    subplot(2,2,3);
    eta_pl = pl.needle.tau ./ max(pl.needle.gamma_dot, 1e-12);
    eta_cr = cr.needle.tau ./ max(cr.needle.gamma_dot, 1e-12);
    plot(pl.needle.r*1e3, eta_pl,'r-','LineWidth',1.8); hold on;
    plot(cr.needle.r*1e3, eta_cr,'b--','LineWidth',1.8);
    set(gca,'YScale','log');
    xlabel('r (mm)'); ylabel('\eta_{app} (Pa.s)'); title('Viscosity across needle exit');
    legend({'PL','Cross'},'Location','best'); grid on;
    % (4) wall shear stress comparison
    subplot(2,2,4);
    bar([pl.syr.tau_w pl.needle.tau_w; cr.syr.tau_w cr.needle.tau_w]');
    set(gca,'XTickLabel',{'Syringe','Needle'});
    ylabel('Wall shear stress (Pa)'); title('Wall shear stress');
    legend({'PL','Cross'},'Location','best'); grid on;
    sgtitle(['System analysis -- ' ttl],'Interpreter','none');
    exportgraphics(f2, fullfile(folder, sprintf('%s_%s_system.png', safe_name, geom.label)), 'Resolution', 200);
    close(f2);
end

%% -------------------- ACROSS-Vp SUMMARY FIGURE --------------------
function make_summary_figure(out, safe_name, geom, Vp_vec, PL, CR, ...
        v_print_PL, v_print_CR, k_flow_PL, k_flow_CR, Rrec_pct)
    vmm = Vp_vec(:)*1e3;
    f = figure('Visible','off','Position',[120 120 1100 800]);
    subplot(2,2,1);
    plot(vmm, PL.dP_total/1e3,'-o','LineWidth',1.4); hold on;
    plot(vmm, CR.dP_total/1e3,'-s','LineWidth',1.4);
    xlabel('V_p (mm/s)'); ylabel('\Delta P_{total} (kPa)'); title('Total extrusion pressure');
    legend({'PL','Cross'},'Location','best'); grid on;
    subplot(2,2,2);
    plot(vmm, PL.gamma_w_needle,'-o','LineWidth',1.4); hold on;
    plot(vmm, CR.gamma_w_needle,'-s','LineWidth',1.4);
    set(gca,'YScale','log'); xlabel('V_p (mm/s)'); ylabel('\gamma_w needle (1/s)');
    title('Needle wall shear rate'); legend({'PL','Cross'},'Location','best'); grid on;
    subplot(2,2,3);
    plot(vmm, v_print_PL*1e3,'-o','LineWidth',1.4); hold on;
    plot(vmm, v_print_CR*1e3,'-s','LineWidth',1.4);
    xlabel('V_p (mm/s)'); ylabel('v_{print} (mm/s)'); title('Predicted head speed');
    legend({'PL','Cross'},'Location','best'); grid on;
    subplot(2,2,4);
    bar([k_flow_PL k_flow_CR]); set(gca,'XTickLabel',{'PL','Cross'});
    ylabel('k_{flow}'); title(sprintf('Flow factor (R_{rec}=%.0f%%)', Rrec_pct)); grid on;
    sgtitle(sprintf('Vp-sweep summary -- %s | %s', safe_name, geom.label),'Interpreter','none');
    exportgraphics(f, fullfile(out, sprintf('%s_%s_summary.png', safe_name, geom.label)), 'Resolution', 200);
    close(f);
end
