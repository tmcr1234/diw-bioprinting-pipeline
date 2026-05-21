function results = bioprinting_algorithm_conical(Rs, Rn, Ln, Ls, Vp, K, n, rho, varargin)
%BIOPRINTING_ALGORITHM_4
% Steady, laminar syringe-transition-needle flow analysis for a Power Law fluid
% in extrusion-based bioprinting. 
% *NOTE: Transition length (Lt) is set to equal the needle length (Ln).*
%
% INPUTS
%   Rs   - Syringe radius (m)
%   Rn   - Needle radius (m)
%   Ln   - Needle length (m) AND Transition length (m)
%   Ls   - Syringe length (m)
%   Vp   - Piston velocity (m/s)
%   K    - Power-law consistency index (Pa·s^n)
%   n    - Power-law index (dimensionless)
%   rho  - Fluid density (kg/m^3)
%
% OPTIONAL NAME-VALUE PAIRS
%   'PlotResults'       - true/false (default: true)
%   'NumPoints'         - number of radial points (default: 200)
%   'Name'              - sample name (default: "Unknown")
%   'SaveData'          - true/false for TXT export (default: true)
%   'SaveFigures'       - true/false for figure export (default: true)
%   'Orientation'       - 'horizontal', 'upward', or 'downward'
%                         default: 'horizontal'
%   'OutputFolder'      - folder for exported files (default: pwd)
%   'IncludeHydrostatic'- true/false (default: true)

%% -------------------- INPUT PARSING --------------------
p = inputParser;
addRequired(p, 'Rs',  @(x) validateattributes(x, {'numeric'}, {'scalar','real','positive'}));
addRequired(p, 'Rn',  @(x) validateattributes(x, {'numeric'}, {'scalar','real','positive'}));
addRequired(p, 'Ln',  @(x) validateattributes(x, {'numeric'}, {'scalar','real','>=', 0}));
addRequired(p, 'Ls',  @(x) validateattributes(x, {'numeric'}, {'scalar','real','>=', 0}));
addRequired(p, 'Vp',  @(x) validateattributes(x, {'numeric'}, {'scalar','real','positive'}));
addRequired(p, 'K',   @(x) validateattributes(x, {'numeric'}, {'scalar','real','positive'}));
addRequired(p, 'n',   @(x) validateattributes(x, {'numeric'}, {'scalar','real','positive'}));
addRequired(p, 'rho', @(x) validateattributes(x, {'numeric'}, {'scalar','real','positive'}));

addParameter(p, 'PlotResults', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'NumPoints', 200, @(x) isnumeric(x) && isscalar(x) && x >= 20);
addParameter(p, 'Name', 'Unknown', @(x) ischar(x) || isstring(x));
addParameter(p, 'SaveData', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'SaveFigures', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Orientation', 'horizontal', @(x) any(validatestring(x, {'horizontal','upward','downward'})));
addParameter(p, 'OutputFolder', pwd, @(x) ischar(x) || isstring(x));
addParameter(p, 'IncludeHydrostatic', true, @(x) islogical(x) && isscalar(x));

parse(p, Rs, Rn, Ln, Ls, Vp, K, n, rho, varargin{:});

% Assign transition length to equal needle length
Lt = Ln; 

plot_results       = p.Results.PlotResults;
num_points         = p.Results.NumPoints;
sample_name        = char(p.Results.Name);
save_data          = p.Results.SaveData;
save_figures       = p.Results.SaveFigures;
orientation        = validatestring(p.Results.Orientation, {'horizontal','upward','downward'});
output_folder      = char(p.Results.OutputFolder);
include_hydrostatic = p.Results.IncludeHydrostatic;

if Rn >= Rs
    error('Needle radius Rn must be smaller than syringe radius Rs.');
end

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
safe_name = regexprep(sample_name, '[<>:"/\\|?*]', '_');

%% -------------------- HEADER --------------------
fprintf('=== BIOPRINTING FLOW ANALYSIS - POWER LAW FLUID ===\n\n');
fprintf('Sample: %s\n\n', sample_name);
fprintf('System Parameters:\n');
fprintf('  Syringe radius: %.4f mm\n', Rs*1000);
fprintf('  Needle radius:  %.4f mm\n', Rn*1000);
fprintf('  Syringe length: %.4f mm\n', Ls*1000);
fprintf('  Transition len: %.4f mm (Set equal to Ln)\n', Lt*1000);
fprintf('  Needle length:  %.4f mm\n', Ln*1000);
fprintf('  Piston vel: %.4f mm/s\n', Vp*1000);
fprintf('  Density: %.1f kg/m^3\n', rho);
fprintf('  Power-law index n: %.4f\n', n);
fprintf('  Consistency index K: %.6g Pa·s^n\n', K);
fprintf('  Orientation: %s\n\n', orientation);

%% -------------------- STEP 1: FLOW RATE --------------------
Q = pi * Rs^2 * Vp;

%% -------------------- STEP 2: PRESSURE ANALYSIS --------------------
abs_dpdz_needle  = 2*K * ((Q * (3*n + 1)) / (pi * n * Rn^((3*n + 1)/n)))^n;
abs_dpdz_syringe = 2*K * ((Q * (3*n + 1)) / (pi * n * Rs^((3*n + 1)/n)))^n;

pressure_drop_needle  = abs_dpdz_needle  * Ln;
pressure_drop_syringe = abs_dpdz_syringe * Ls;

% Transition Zone Pressure Drop (Lubrication approximation for conical taper)
if Lt > 0
    term1 = (2 * K * Lt) / (3 * n * (Rs - Rn));
    term2 = ((Q * (3 * n + 1)) / (pi * n))^n;
    term3 = (1 / Rn^(3 * n)) - (1 / Rs^(3 * n));
    pressure_drop_transition = term1 * term2 * term3;
else
    pressure_drop_transition = 0;
end

total_viscous_pressure_drop = pressure_drop_syringe + pressure_drop_transition + pressure_drop_needle;

g = 9.81;
switch orientation
    case 'horizontal'
        orientation_factor = 0;
    case 'upward'
        orientation_factor = +1;  
    case 'downward'
        orientation_factor = -1;  
end

if include_hydrostatic
    hydrostatic_pressure_drop = orientation_factor * rho * g * (Ls + Lt + Ln);
else
    hydrostatic_pressure_drop = 0;
end

total_system_pressure_drop = total_viscous_pressure_drop + hydrostatic_pressure_drop;

% Absolute Pressure Calculations 
P_atm = 101325; 
abs_pressure_exit = P_atm; 
abs_pressure_needle_inlet = abs_pressure_exit + pressure_drop_needle;
abs_pressure_transition_inlet = abs_pressure_needle_inlet + pressure_drop_transition;
abs_pressure_syringe_inlet = abs_pressure_transition_inlet + pressure_drop_syringe;

fprintf('STEP 2 - Pressure Analysis:\n');
fprintf('  Pressure drop across needle:     %.6f Pa\n', pressure_drop_needle);
fprintf('  Pressure drop across transition: %.6f Pa\n', pressure_drop_transition);
fprintf('  Pressure drop across syringe:    %.6f Pa\n', pressure_drop_syringe);
fprintf('  Total viscous pressure drop:     %.6f Pa\n', total_viscous_pressure_drop);

%% -------------------- STEP 3: SHEAR RATES AND SHEAR STRESSES --------------------
gamma_dot_app_needle = 4 * Q / (pi * Rn^3);
gamma_dot_app_syringe = 4 * Q / (pi * Rs^3);

gamma_dot_true_needle = gamma_dot_app_needle * (3*n + 1) / (4*n);
gamma_dot_true_syringe = gamma_dot_app_syringe * (3*n + 1) / (4*n);

tau_wall_needle = abs_dpdz_needle * Rn / 2;
tau_wall_syringe = abs_dpdz_syringe * Rs / 2;

%% -------------------- STEP 4: VELOCITY PROFILES --------------------
r_needle = linspace(0, Rn, num_points);
r_syringe = linspace(0, Rs, num_points);

u_needle = (n/(n+1)) * (abs_dpdz_needle/(2*K))^(1/n) .* (Rn^((n+1)/n) - r_needle.^((n+1)/n));
u_syringe = (n/(n+1)) * (abs_dpdz_syringe/(2*K))^(1/n) .* (Rs^((n+1)/n) - r_syringe.^((n+1)/n));

u_max_needle = u_needle(1);
u_max_syringe = u_syringe(1);
u_avg_needle = Q / (pi * Rn^2);
u_avg_syringe = Q / (pi * Rs^2);

%% -------------------- STEP 5: REYNOLDS NUMBER --------------------
D_needle = 2*Rn;
D_syringe = 2*Rs;
Re_needle = (rho * u_avg_needle^(2-n) * D_needle^n) / (K * 8^(n-1) * ((3*n+1)/(4*n))^n);
Re_syringe = (rho * u_avg_syringe^(2-n) * D_syringe^n) / (K * 8^(n-1) * ((3*n+1)/(4*n))^n);
Re_crit = 2100 * (n^0.75);

if Re_needle < Re_crit && Re_syringe < Re_crit
    flow_regime_message = 'Likely laminar based on generalized Reynolds-number estimate';
else
    flow_regime_message = 'Flow may be transitional';
end

%% -------------------- STEP 6 & 7: RADIAL SHEAR DISTRIBUTIONS --------------------
gamma_dot_needle = (abs_dpdz_needle/(2*K))^(1/n) .* r_needle.^(1/n);
gamma_dot_syringe = (abs_dpdz_syringe/(2*K))^(1/n) .* r_syringe.^(1/n);
gamma_dot_needle(1) = 0;
gamma_dot_syringe(1) = 0;

tau_needle = K .* gamma_dot_needle.^n;
tau_syringe = K .* gamma_dot_syringe.^n;

%% -------------------- STORE RESULTS --------------------
results.sample_name = sample_name;
results.safe_name = safe_name;
results.output_folder = output_folder;
results.model.type = 'Power Law';
results.model.assumptions = {'steady flow', 'laminar flow', 'incompressible fluid', 'fully developed flow', 'lubrication approx for transition'};
results.notes.hydrostatic_included = include_hydrostatic;
results.notes.orientation = orientation;
results.notes.reynolds_comment = flow_regime_message;

results.system_params.Rs = Rs;
results.system_params.Rn = Rn;
results.system_params.Ln = Ln;
results.system_params.Ls = Ls;
results.system_params.Lt = Lt;
results.system_params.Vp = Vp;
results.system_params.K = K;
results.system_params.n = n;
results.system_params.rho = rho;

results.flow.Q = Q;
results.flow.u_avg_needle = u_avg_needle;
results.flow.u_avg_syringe = u_avg_syringe;
results.flow.u_max_needle = u_max_needle;
results.flow.u_max_syringe = u_max_syringe;

results.pressure.abs_dpdz_needle = abs_dpdz_needle;
results.pressure.abs_dpdz_syringe = abs_dpdz_syringe;
results.pressure.drop_needle = pressure_drop_needle;
results.pressure.drop_transition = pressure_drop_transition;
results.pressure.drop_syringe = pressure_drop_syringe;
results.pressure.drop_total_viscous = total_viscous_pressure_drop;
results.pressure.drop_hydrostatic = hydrostatic_pressure_drop;
results.pressure.drop_total_system = total_system_pressure_drop;

results.pressure.abs_exit = abs_pressure_exit;
results.pressure.abs_needle_inlet = abs_pressure_needle_inlet;
results.pressure.abs_transition_inlet = abs_pressure_transition_inlet;
results.pressure.abs_syringe_inlet = abs_pressure_syringe_inlet;

results.shear.tau_wall_needle = tau_wall_needle;
results.shear.tau_wall_syringe = tau_wall_syringe;
results.reynolds.Re_needle = Re_needle;
results.reynolds.Re_syringe = Re_syringe;
results.reynolds.Re_crit_estimated = Re_crit;

results.profiles.r_needle = r_needle;
results.profiles.r_syringe = r_syringe;
results.profiles.u_needle = u_needle;
results.profiles.u_syringe = u_syringe;
results.profiles.gamma_dot_needle = gamma_dot_needle;
results.profiles.gamma_dot_syringe = gamma_dot_syringe;
results.profiles.tau_needle = tau_needle;
results.profiles.tau_syringe = tau_syringe;

%% -------------------- PLOTS & SAVING --------------------
if plot_results
    create_bioprinting_plots_improved(results, save_figures);
end
if save_data
    save_bioprinting_data_improved(results);
end
fprintf('=== ANALYSIS COMPLETE ===\n');
end

%% ========================================================================
function create_bioprinting_plots_improved(results, save_figs)
if nargin < 2; save_figs = false; end

Rs = results.system_params.Rs;
Rn = results.system_params.Rn;
Ls = results.system_params.Ls;
Lt = results.system_params.Lt;
Ln = results.system_params.Ln;
K = results.system_params.K;
n = results.system_params.n;
Q = results.flow.Q;

sample_name = results.sample_name;
safe_name = results.safe_name;
output_folder = results.output_folder;

%% Figure 4: Summary (Showing Pressure Plot adjustments)
f4 = figure('Position', [250, 250, 1100, 650]);

subplot(2,2,1)
% 1. Syringe part
z1 = linspace(0, Ls, 101);
P1 = results.pressure.abs_syringe_inlet - results.pressure.abs_dpdz_syringe * z1;
plot(z1*1000, P1, 'b-', 'LineWidth', 2); hold on;

% 2. Transition part (Non-linear pressure drop)
if Lt > 0
    z_trans = linspace(Ls, Ls + Lt, 101);
    R_trans = Rs - (Rs - Rn) * (z_trans - Ls) / Lt;
    dpdz_trans = 2*K * ((Q * (3*n + 1)) ./ (pi * n * R_trans.^((3*n + 1)/n))).^n;
    P_trans = results.pressure.abs_transition_inlet - cumtrapz(z_trans, dpdz_trans);
    plot(z_trans*1000, P_trans, 'g-', 'LineWidth', 2);
end

% 3. Needle part
z2 = linspace(Ls + Lt, Ls + Lt + Ln, 101);
P2 = results.pressure.abs_needle_inlet - results.pressure.abs_dpdz_needle * (z2 - (Ls + Lt));
plot(z2*1000, P2, 'r-', 'LineWidth', 2);

xlabel('Axial Position (mm)');
ylabel('Absolute Pressure (Pa)');
title(['Absolute Pressure Distribution - ', sample_name]);
legend('Syringe', 'Transition', 'Needle', 'Location', 'best');
grid on;

subplot(2,2,2)
bar([results.reynolds.Re_syringe, results.reynolds.Re_needle]); hold on;
yline(results.reynolds.Re_crit_estimated, 'k--', 'LineWidth', 2);
set(gca, 'XTickLabel', {'Syringe','Needle'});
ylabel('Generalized Reynolds Number'); title('Flow Regime Estimate'); grid on;

subplot(2,2,3)
velocities = [results.flow.u_avg_syringe*1000, results.flow.u_max_syringe*1000, ...
              results.flow.u_avg_needle*1000, results.flow.u_max_needle*1000];
bar(velocities);
set(gca, 'XTickLabel', {'Syr Avg','Syr Max','Need Avg','Need Max'});
ylabel('Velocity (mm/s)'); title('Velocity Comparison'); grid on;

subplot(2,2,4)
bar([results.shear.tau_wall_syringe, results.shear.tau_wall_needle]);
set(gca, 'XTickLabel', {'Syringe Wall','Needle Wall'});
ylabel('Wall Shear Stress (Pa)'); title('Wall Shear Stress'); grid on;

sgtitle(['Bioprinting System Summary - ', sample_name]);
if save_figs
    exportgraphics(f4, fullfile(output_folder, sprintf('%s - System Overview.png', safe_name)), 'Resolution', 1200);
end
end

%% ========================================================================
function save_bioprinting_data_improved(results)
filename = fullfile(results.output_folder, [results.safe_name, '_data.txt']);
fid = fopen(filename, 'w');
if fid == -1; return; end

fprintf(fid, '=== BIOPRINTING ANALYSIS RESULTS ===\n');
fprintf(fid, 'Sample: %s\n\n', results.sample_name);
fprintf(fid, 'System Parameters:\n');
fprintf(fid, '  Syringe radius (mm): %.6f\n', results.system_params.Rs*1000);
fprintf(fid, '  Needle radius (mm): %.6f\n', results.system_params.Rn*1000);
fprintf(fid, '  Syringe length (mm): %.6f\n', results.system_params.Ls*1000);
fprintf(fid, '  Transition len(mm): %.6f (Equal to Ln)\n', results.system_params.Lt*1000);
fprintf(fid, '  Needle length (mm): %.6f\n', results.system_params.Ln*1000);
fprintf(fid, '  Piston vel (mm/s): %.6f\n', results.system_params.Vp*1000);

fprintf(fid, '\nPressure Analysis:\n');
fprintf(fid, '  Syringe pressure drop (Pa): %.9e\n', results.pressure.drop_syringe);
fprintf(fid, '  Transition pressure drop (Pa): %.9e\n', results.pressure.drop_transition);
fprintf(fid, '  Needle pressure drop (Pa): %.9e\n', results.pressure.drop_needle);
fprintf(fid, '  Total viscous pressure drop (Pa): %.9e\n', results.pressure.drop_total_viscous);

fclose(fid);
fprintf('Data saved to file: %s\n', filename);
end