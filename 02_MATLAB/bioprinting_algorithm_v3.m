function results = bioprinting_algorithm_v3(sample, geom, Vp_vec, varargin)
%BIOPRINTING_ALGORITHM_V3  Unified Power-Law + Cross extrusion solver with slicer-parameter output.
%
% Replaces v2 / v3 / cross_v2 with a single algorithm that:
%   (1) accepts a vector of piston velocities,
%   (2) computes Q, tau_w, gamma_w, dP for BOTH Power-Law and Cross models,
%   (3) derives the FFF-slicer-equivalent parameters
%       (v_print, w_line, h_layer, beta_swell, k_flow), and
%   (4) writes a per-sample CSV lookup table.
%
% Geometry assumed: straight cylindrical needle.
%
% --------------------------------------------------------------------
% INPUTS
%   sample - struct with rheology + recovery fields:
%       .name        sample identifier (string)
%       .K_PL        Power-Law consistency (Pa*s^n)
%       .n_PL        Power-Law flow index (-)
%       .eta0        Cross zero-shear viscosity (Pa*s)
%       .etaInf      Cross infinite-shear viscosity (Pa*s)
%       .lambda      Cross time constant (s)
%       .m_Cross     Cross exponent (-)
%       .Rrec_pct    structural recovery (%) at deposition shear rate
%
%   geom - struct:
%       .Rs          syringe inner radius (m)
%       .R_n         needle inner radius (m)
%       .L_n         needle length (m)
%       .label       label, e.g. '21G' or '22G'
%       .h_factor    layer height as fraction of needle ID (default 0.7)
%
%   Vp_vec - vector of piston velocities (m/s)
%
% NAME-VALUE OPTIONAL
%   'OutputFolder' (default 'output_v3')
%   'SaveCSV'      (default true)
%   'PlotResults'  (default true)
%   'f_slip'       (default 1.0)   syringe-piston slip factor
%   'NumPoints'    (default 400)   for Cross numerical integration
%
% OUTPUT
%   results - struct with .PL, .Cross sub-structs, each containing
%       Q, tau_w, gamma_w, dP, v_mean, v_print, w_line, h_layer,
%       beta_swell, k_flow, gamma_w_dep_estimate
%   Plus results.sample, results.geom, results.Vp
%
% Author: T.M.C. Rodrigues (PEMM/COPPE/UFRJ) -- 2026-05-14
% --------------------------------------------------------------------

%% -------------------- INPUT PARSING --------------------
p = inputParser;
addRequired(p, 'sample', @isstruct);
addRequired(p, 'geom',   @isstruct);
addRequired(p, 'Vp_vec', @(x) isnumeric(x) && all(x > 0));
addParameter(p, 'OutputFolder', 'output_v3', @(x) ischar(x)||isstring(x));
addParameter(p, 'SaveCSV',      true, @islogical);
addParameter(p, 'PlotResults',  true, @islogical);
addParameter(p, 'f_slip',       1.0,  @(x) isnumeric(x) && x > 0);
addParameter(p, 'NumPoints',    400,  @(x) isnumeric(x) && x > 10);
parse(p, sample, geom, Vp_vec, varargin{:});

output_folder = char(p.Results.OutputFolder);
save_csv      = p.Results.SaveCSV;
plot_results  = p.Results.PlotResults;
f_slip        = p.Results.f_slip;
num_points    = p.Results.NumPoints;

if ~exist(output_folder, 'dir'), mkdir(output_folder); end
safe_name = regexprep(char(sample.name), '[<>:"/\\|?*\.\s]', '_');

%% -------------------- GEOMETRY DEFAULTS --------------------
if ~isfield(geom, 'h_factor'),  geom.h_factor = 0.7;  end
R   = geom.R_n;
L   = geom.L_n;
Rs  = geom.Rs;
h_layer = geom.h_factor * 2 * R;   % m

%% -------------------- VOLUMETRIC FLOW --------------------
A_syringe = pi * Rs^2;
A_needle  = pi * R^2;
Q_vec     = A_syringe * Vp_vec;          % m^3/s
v_mean    = Q_vec / A_needle;            % m/s in needle

%% -------------------- POWER-LAW SOLUTION --------------------
K  = sample.K_PL;
n  = sample.n_PL;

%   Q = (pi R^3 n / (3n+1)) * (tau_w / K)^(1/n)
%   => tau_w = K * ((3n+1) Q / (pi R^3 n))^n
tau_w_PL  = K * ( (3*n + 1) .* Q_vec ./ (pi * R^3 * n) ).^n;
gamma_w_PL = (tau_w_PL ./ K).^(1/n);
dP_PL     = 2 * L * tau_w_PL / R;        % Pa

%% -------------------- CROSS SOLUTION --------------------
eta0   = sample.eta0;
etaInf = sample.etaInf;
lambda = sample.lambda;
mC     = sample.m_Cross;

n_v = numel(Vp_vec);
tau_w_C   = zeros(size(Q_vec));
gamma_w_C = zeros(size(Q_vec));
dP_C      = zeros(size(Q_vec));
for ii = 1:n_v
    Qi = Q_vec(ii);
    % tau_w from inversion of Weissenberg-Rabinowitsch integral
    tau_max = 1e8;   % upper bracket (Pa)
    fQ = @(tw) Q_of_tau_w_cross(tw, R, eta0, etaInf, lambda, mC, num_points) - Qi;
    % bracket
    Q_max = Q_of_tau_w_cross(tau_max, R, eta0, etaInf, lambda, mC, num_points);
    if Q_max < Qi
        tau_w_C(ii)   = NaN;
        gamma_w_C(ii) = NaN;
        dP_C(ii)      = NaN;
        continue
    end
    tw = fzero(fQ, [1e-3, tau_max]);
    tau_w_C(ii)   = tw;
    gamma_w_C(ii) = invert_cross_tau(tw, eta0, etaInf, lambda, mC);
    dP_C(ii)      = 2 * L * tw / R;
end

%% -------------------- SLICER LOOKUP --------------------
% Die-swell proxy: linear in (1 - flow_index_proxy)
%   PL:    beta = 0.30 * (1 - n_PL)        (n_PL near 1 -> low swell)
%   Cross: beta = 0.30 * (1 - m_Cross)     (m near 1 -> low swell)
beta_PL    = max(0, 0.30 * (1 - n));
beta_Cross = max(0, 0.30 * (1 - mC));

w_line_PL    = 2 * R * (1 + beta_PL);
w_line_Cross = 2 * R * (1 + beta_Cross);

v_print_PL    = Q_vec ./ (w_line_PL    * h_layer);   % m/s
v_print_Cross = Q_vec ./ (w_line_Cross * h_layer);   % m/s

Rrec_pct = sample.Rrec_pct;
k_flow_PL    = (1 + beta_PL).^2    * f_slip * sqrt(Rrec_pct/100);
k_flow_Cross = (1 + beta_Cross).^2 * f_slip * sqrt(Rrec_pct/100);

%% -------------------- ASSEMBLE RESULTS --------------------
results.sample   = sample;
results.geom     = geom;
results.Vp_mm_s  = Vp_vec(:) * 1e3;
results.Q_mm3_s  = Q_vec(:)  * 1e9;

results.PL = struct( ...
    'tau_w_Pa',     tau_w_PL(:), ...
    'gamma_w_invs', gamma_w_PL(:), ...
    'dP_kPa',       dP_PL(:)/1e3, ...
    'v_mean_mm_s',  v_mean(:)*1e3, ...
    'v_print_mm_s', v_print_PL(:)*1e3, ...
    'w_line_mm',    w_line_PL*1e3, ...
    'h_layer_mm',   h_layer*1e3, ...
    'beta_swell',   beta_PL, ...
    'k_flow',       k_flow_PL);

results.Cross = struct( ...
    'tau_w_Pa',     tau_w_C(:), ...
    'gamma_w_invs', gamma_w_C(:), ...
    'dP_kPa',       dP_C(:)/1e3, ...
    'v_mean_mm_s',  v_mean(:)*1e3, ...
    'v_print_mm_s', v_print_Cross(:)*1e3, ...
    'w_line_mm',    w_line_Cross*1e3, ...
    'h_layer_mm',   h_layer*1e3, ...
    'beta_swell',   beta_Cross, ...
    'k_flow',       k_flow_Cross);

%% -------------------- CSV EXPORT --------------------
if save_csv
    nrows = length(Vp_vec);
    T = table( ...
        repmat({safe_name}, nrows, 1), ...
        repmat({geom.label}, nrows, 1), ...
        Vp_vec(:)*1e3, ...
        Q_vec(:)*1e9, ...
        v_mean(:)*1e3, ...
        tau_w_PL(:),    gamma_w_PL(:),    dP_PL(:)/1e3,    v_print_PL(:)*1e3, ...
        repmat(w_line_PL*1e3,    nrows, 1), ...
        repmat(beta_PL,          nrows, 1), ...
        repmat(k_flow_PL,        nrows, 1), ...
        tau_w_C(:),     gamma_w_C(:),     dP_C(:)/1e3,     v_print_Cross(:)*1e3, ...
        repmat(w_line_Cross*1e3, nrows, 1), ...
        repmat(beta_Cross,       nrows, 1), ...
        repmat(k_flow_Cross,     nrows, 1), ...
        repmat(h_layer*1e3,      nrows, 1), ...
        'VariableNames', { ...
            'sample','needle', ...
            'Vp_mm_s','Q_mm3_s','v_mean_mm_s', ...
            'PL_tau_w_Pa','PL_gamma_w_invs','PL_dP_kPa','PL_v_print_mm_s','PL_w_line_mm','PL_beta_swell','PL_k_flow', ...
            'Cross_tau_w_Pa','Cross_gamma_w_invs','Cross_dP_kPa','Cross_v_print_mm_s','Cross_w_line_mm','Cross_beta_swell','Cross_k_flow', ...
            'h_layer_mm'});
    csv_name = fullfile(output_folder, ...
                        sprintf('slicer_lookup_%s_%s.csv', safe_name, geom.label));
    writetable(T, csv_name);
    fprintf('[Saved CSV] %s\n', csv_name);
end

%% -------------------- PLOTS --------------------
if plot_results
    fig = figure('Visible','off','Position',[100 100 1100 800]);

    subplot(2,2,1);
    plot(Vp_vec*1e3, dP_PL/1e3, '-o', 'LineWidth', 1.4); hold on
    plot(Vp_vec*1e3, dP_C /1e3, '-s', 'LineWidth', 1.4);
    xlabel('Piston velocity V_p (mm/s)');
    ylabel('\Delta P (kPa)');
    title(sprintf('%s -- %s: extrusion pressure', safe_name, geom.label), 'Interpreter','none');
    legend({'Power Law','Cross'},'Location','best');
    grid on

    subplot(2,2,2);
    plot(Vp_vec*1e3, gamma_w_PL, '-o', 'LineWidth', 1.4); hold on
    plot(Vp_vec*1e3, gamma_w_C , '-s', 'LineWidth', 1.4);
    set(gca,'YScale','log');
    xlabel('Piston velocity V_p (mm/s)');
    ylabel('Wall shear rate \dot\gamma_w (1/s)');
    title('Wall shear rate'); legend({'PL','Cross'},'Location','best'); grid on

    subplot(2,2,3);
    plot(Vp_vec*1e3, results.PL.v_print_mm_s,    '-o', 'LineWidth', 1.4); hold on
    plot(Vp_vec*1e3, results.Cross.v_print_mm_s, '-s', 'LineWidth', 1.4);
    xlabel('Piston velocity V_p (mm/s)');
    ylabel('Head speed v_{print} (mm/s)');
    title(sprintf('Predicted print speed (h_{layer}=%.2f mm)', h_layer*1e3));
    legend({'PL','Cross'},'Location','best'); grid on

    subplot(2,2,4);
    bar([k_flow_PL k_flow_Cross]);
    set(gca,'XTickLabel',{'PL','Cross'});
    ylim([0 max(1, 1.2*max(k_flow_PL,k_flow_Cross))]);
    ylabel('k_{flow} (predicted)');
    title(sprintf('Flow factor (R_{rec}=%.0f%%)', Rrec_pct));
    text(0.55, 0.92, sprintf('PL: %.3f\nCross: %.3f', k_flow_PL, k_flow_Cross), ...
         'Units','normalized','FontWeight','bold');

    fig_name = fullfile(output_folder, ...
                        sprintf('plots_%s_%s.png', safe_name, geom.label));
    exportgraphics(fig, fig_name, 'Resolution', 200);
    close(fig);
    fprintf('[Saved PNG] %s\n', fig_name);
end

end  % main function

%% ========================================================
%% LOCAL FUNCTIONS
%% ========================================================
function gdot = invert_cross_tau(tau, eta0, etaInf, lambda, m)
% Solve eta(gdot)*gdot = tau for shear rate.
    if tau < 1e-12, gdot = 0; return; end
    f = @(g) (etaInf + (eta0-etaInf)./(1 + (lambda*g).^m)) .* g - tau;
    % Robust bracket
    g_lo = 1e-12; g_hi = 1e8;
    f_lo = f(g_lo); f_hi = f(g_hi);
    if f_lo * f_hi > 0
        % Try expanding bracket
        for k = 1:6
            g_hi = g_hi * 10;
            f_hi = f(g_hi);
            if f_lo * f_hi < 0, break, end
        end
    end
    gdot = fzero(f, [g_lo, g_hi]);
end

function Q = Q_of_tau_w_cross(tau_w, R, eta0, etaInf, lambda, m, npts)
% Weissenberg-Rabinowitsch: Q = (pi R^3 / tau_w^3) integral_0^{tau_w} tau^2 gdot(tau) dtau
    tau_grid = linspace(0, tau_w, npts);
    gdot_grid = zeros(size(tau_grid));
    for k = 2:npts
        gdot_grid(k) = invert_cross_tau(tau_grid(k), eta0, etaInf, lambda, m);
    end
    integrand = (tau_grid.^2) .* gdot_grid;
    Q = (pi * R^3 / tau_w^3) * trapz(tau_grid, integrand);
end
