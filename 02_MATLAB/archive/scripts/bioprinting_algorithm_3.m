function results = bioprinting_algorithm_3(Rs, Rn, Ln, Ls, Vp, K, n, rho, varargin)
%BIOPRINTING_ALGORITHM_3
% Steady, laminar syringe-needle flow analysis for a Power Law fluid
% in extrusion-based bioprinting.
%
% INPUTS
%   Rs   - Syringe radius (m)
%   Rn   - Needle radius (m)
%   Ln   - Needle length (m)
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
%
% OUTPUT
%   results - structure containing computed quantities
%
% MODEL NOTES
%   This model assumes:
%   - incompressible flow
%   - steady flow
%   - laminar flow
%   - fully developed flow in syringe and needle
%   - no wall slip
%   - no entrance/contraction losses
%   - no yield stress
%   - no viscoelasticity or thixotropy
%   - Exit pressure is 1 atm (101325 Pa)

%% -------------------- INPUT PARSING --------------------
p = inputParser;
addRequired(p, 'Rs',  @(x) validateattributes(x, {'numeric'}, {'scalar','real','positive'}));
addRequired(p, 'Rn',  @(x) validateattributes(x, {'numeric'}, {'scalar','real','positive'}));
addRequired(p, 'Ln',  @(x) validateattributes(x, {'numeric'}, {'scalar','real','positive'}));
addRequired(p, 'Ls',  @(x) validateattributes(x, {'numeric'}, {'scalar','real','positive'}));
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
if n >= 1
    warning('n >= 1. The fluid is not shear-thinning. Confirm whether the Power Law model is appropriate.');
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
fprintf('  Needle length:  %.4f mm\n', Ln*1000);
fprintf('  Syringe length: %.4f mm\n', Ls*1000);
fprintf('  Piston velocity: %.4f mm/s\n', Vp*1000);
fprintf('  Density: %.1f kg/m^3\n', rho);
fprintf('  Power-law index n: %.4f\n', n);
fprintf('  Consistency index K: %.6g Pa·s^n\n', K);
fprintf('  Orientation: %s\n\n', orientation);

%% -------------------- STEP 1: FLOW RATE --------------------
Q = pi * Rs^2 * Vp;
fprintf('STEP 1 - Flow Rate Calculation:\n');
fprintf('  Q = pi*Rs^2*Vp = %.6e m^3/s\n', Q);
fprintf('  Q = %.6f mm^3/s\n\n', Q*1e9);

%% -------------------- STEP 2: PRESSURE ANALYSIS --------------------
% Power-law pipe flow relation:
% Q = pi * R^((3n+1)/n) * (n/(3n+1)) * (abs(dp/dz)/(2K))^(1/n)
% Rearranged:
% abs(dp/dz) = 2K * [ Q(3n+1)/(pi*n*R^((3n+1)/n)) ]^n

abs_dpdz_needle  = 2*K * ((Q * (3*n + 1)) / (pi * n * Rn^((3*n + 1)/n)))^n;
abs_dpdz_syringe = 2*K * ((Q * (3*n + 1)) / (pi * n * Rs^((3*n + 1)/n)))^n;

pressure_drop_needle  = abs_dpdz_needle  * Ln;
pressure_drop_syringe = abs_dpdz_syringe * Ls;

total_viscous_pressure_drop = pressure_drop_syringe + pressure_drop_needle;

g = 9.81;
switch orientation
    case 'horizontal'
        orientation_factor = 0;
    case 'upward'
        orientation_factor = +1;  % flow against gravity
    case 'downward'
        orientation_factor = -1;  % gravity assists flow
end

if include_hydrostatic
    hydrostatic_pressure_drop = orientation_factor * rho * g * (Ls + Ln);
else
    hydrostatic_pressure_drop = 0;
end

total_system_pressure_drop = total_viscous_pressure_drop + hydrostatic_pressure_drop;

% --- NEW: Absolute Pressure Calculations ---
P_atm = 101325; % 1 atm in Pascals
abs_pressure_exit = P_atm; 
abs_pressure_needle_inlet = abs_pressure_exit + pressure_drop_needle;
abs_pressure_syringe_inlet = abs_pressure_needle_inlet + pressure_drop_syringe;

fprintf('STEP 2 - Pressure Analysis:\n');
fprintf('  |dp/dz| in needle:  %.3e Pa/m\n', abs_dpdz_needle);
fprintf('  |dp/dz| in syringe: %.3e Pa/m\n', abs_dpdz_syringe);
fprintf('  Pressure drop across needle:  %.6f Pa\n', pressure_drop_needle);
fprintf('  Pressure drop across syringe: %.6f Pa\n', pressure_drop_syringe);
fprintf('  Total viscous pressure drop:  %.6f Pa\n', total_viscous_pressure_drop);
fprintf('  Hydrostatic contribution:     %.6f Pa\n', hydrostatic_pressure_drop);
fprintf('  Total estimated system pressure drop: %.6f Pa\n\n', total_system_pressure_drop);
fprintf('  --- Absolute Pressures (Viscous) ---\n');
fprintf('  Exit pressure (Air):          %.2f Pa (1 atm)\n', abs_pressure_exit);
fprintf('  Needle inlet pressure:        %.2f Pa\n', abs_pressure_needle_inlet);
fprintf('  Syringe inlet (piston) pres.: %.2f Pa\n\n', abs_pressure_syringe_inlet);

%% -------------------- STEP 3: SHEAR RATES AND SHEAR STRESSES --------------------
gamma_dot_app_needle = 4 * Q / (pi * Rn^3);
gamma_dot_app_syringe = 4 * Q / (pi * Rs^3);

% Rabinowitsch correction for power-law fluid
gamma_dot_true_needle = gamma_dot_app_needle * (3*n + 1) / (4*n);
gamma_dot_true_syringe = gamma_dot_app_syringe * (3*n + 1) / (4*n);

tau_wall_needle = abs_dpdz_needle * Rn / 2;
tau_wall_syringe = abs_dpdz_syringe * Rs / 2;

fprintf('STEP 3 - Shear Rate and Stress Analysis:\n');
fprintf('  Needle:\n');
fprintf('    Apparent wall shear rate: %.6f s^-1\n', gamma_dot_app_needle);
fprintf('    True wall shear rate:     %.6f s^-1\n', gamma_dot_true_needle);
fprintf('    Wall shear stress:        %.6f Pa\n', tau_wall_needle);
fprintf('  Syringe:\n');
fprintf('    Apparent wall shear rate: %.6f s^-1\n', gamma_dot_app_syringe);
fprintf('    True wall shear rate:     %.6f s^-1\n', gamma_dot_true_syringe);
fprintf('    Wall shear stress:        %.6f Pa\n\n', tau_wall_syringe);

%% -------------------- STEP 4: VELOCITY PROFILES --------------------
r_needle = linspace(0, Rn, num_points);
r_syringe = linspace(0, Rs, num_points);

% Fully developed velocity profile for power-law fluid
u_needle = (n/(n+1)) * (abs_dpdz_needle/(2*K))^(1/n) .* ...
           (Rn^((n+1)/n) - r_needle.^((n+1)/n));

u_syringe = (n/(n+1)) * (abs_dpdz_syringe/(2*K))^(1/n) .* ...
            (Rs^((n+1)/n) - r_syringe.^((n+1)/n));

u_max_needle = u_needle(1);
u_max_syringe = u_syringe(1);

u_avg_needle = Q / (pi * Rn^2);
u_avg_syringe = Q / (pi * Rs^2);

fprintf('STEP 4 - Velocity Profile Analysis:\n');
fprintf('  Needle:\n');
fprintf('    Max velocity: %.6f m/s\n', u_max_needle);
fprintf('    Avg velocity: %.6f m/s\n', u_avg_needle);
fprintf('    Max/Avg ratio: %.6f\n', u_max_needle/u_avg_needle);
fprintf('  Syringe:\n');
fprintf('    Max velocity: %.6f m/s\n', u_max_syringe);
fprintf('    Avg velocity: %.6f m/s\n', u_avg_syringe);
fprintf('    Max/Avg ratio: %.6f\n\n', u_max_syringe/u_avg_syringe);

%% -------------------- STEP 5: REYNOLDS NUMBER --------------------
D_needle = 2*Rn;
D_syringe = 2*Rs;

Re_needle = (rho * u_avg_needle^(2-n) * D_needle^n) / ...
            (K * 8^(n-1) * ((3*n+1)/(4*n))^n);

Re_syringe = (rho * u_avg_syringe^(2-n) * D_syringe^n) / ...
             (K * 8^(n-1) * ((3*n+1)/(4*n))^n);

Re_crit = 2100 * (n^0.75);  % empirical estimate only

if Re_needle < Re_crit && Re_syringe < Re_crit
    flow_regime_message = 'Likely laminar based on generalized Reynolds-number estimate';
else
    flow_regime_message = 'Flow may be transitional; check assumptions and correlations carefully';
end

fprintf('STEP 5 - Flow Regime Analysis:\n');
fprintf('  Generalized Reynolds number (needle):  %.6f\n', Re_needle);
fprintf('  Generalized Reynolds number (syringe): %.6f\n', Re_syringe);
fprintf('  Critical Reynolds estimate: %.6f\n', Re_crit);
fprintf('  Interpretation: %s\n\n', flow_regime_message);

%% -------------------- STEP 6: RADIAL SHEAR-RATE DISTRIBUTION --------------------
gamma_dot_needle = (abs_dpdz_needle/(2*K))^(1/n) .* r_needle.^(1/n);
gamma_dot_syringe = (abs_dpdz_syringe/(2*K))^(1/n) .* r_syringe.^(1/n);

gamma_dot_needle(1) = 0;
gamma_dot_syringe(1) = 0;

%% -------------------- STEP 7: RADIAL SHEAR-STRESS DISTRIBUTION --------------------
tau_needle = K .* gamma_dot_needle.^n;
tau_syringe = K .* gamma_dot_syringe.^n;

fprintf('STEP 6-7 - Radial Distribution Summary:\n');
fprintf('  Maximum shear rate in needle:  %.6f s^-1\n', max(gamma_dot_needle));
fprintf('  Maximum shear stress in needle: %.6f Pa\n', max(tau_needle));
fprintf('  Maximum shear rate in syringe:  %.6f s^-1\n', max(gamma_dot_syringe));
fprintf('  Maximum shear stress in syringe: %.6f Pa\n\n', max(tau_syringe));

%% -------------------- STORE RESULTS --------------------
results.sample_name = sample_name;
results.safe_name = safe_name;
results.output_folder = output_folder;

results.model.type = 'Power Law';
results.model.assumptions = { ...
    'steady flow', ...
    'laminar flow', ...
    'incompressible fluid', ...
    'fully developed flow in syringe and needle', ...
    'no contraction/entrance losses', ...
    'no wall slip', ...
    'no yield stress', ...
    'no viscoelasticity', ...
    'no thixotropy', ...
    'exit pressure is 1 atm'};

results.notes.contraction_loss_included = false;
results.notes.wall_slip_included = false;
results.notes.yield_stress_included = false;
results.notes.hydrostatic_included = include_hydrostatic;
results.notes.orientation = orientation;
results.notes.reynolds_comment = flow_regime_message;

results.system_params.Rs = Rs;
results.system_params.Rn = Rn;
results.system_params.Ln = Ln;
results.system_params.Ls = Ls;
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
results.pressure.drop_syringe = pressure_drop_syringe;
results.pressure.drop_total_viscous = total_viscous_pressure_drop;
results.pressure.drop_hydrostatic = hydrostatic_pressure_drop;
results.pressure.drop_total_system = total_system_pressure_drop;

% --- NEW: Storing Absolute Pressures ---
results.pressure.abs_exit = abs_pressure_exit;
results.pressure.abs_needle_inlet = abs_pressure_needle_inlet;
results.pressure.abs_syringe_inlet = abs_pressure_syringe_inlet;

results.shear.gamma_dot_app_needle = gamma_dot_app_needle;
results.shear.gamma_dot_true_needle = gamma_dot_true_needle;
results.shear.gamma_dot_app_syringe = gamma_dot_app_syringe;
results.shear.gamma_dot_true_syringe = gamma_dot_true_syringe;
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

%% -------------------- PLOTS --------------------
if plot_results
    create_bioprinting_plots_improved(results, save_figures);
end

%% -------------------- SAVE DATA --------------------
if save_data
    save_bioprinting_data_improved(results);
end

fprintf('=== ANALYSIS COMPLETE ===\n');
end

%% ========================================================================
function create_bioprinting_plots_improved(results, save_figs)
if nargin < 2
    save_figs = false;
end

Rs = results.system_params.Rs;
Rn = results.system_params.Rn;
r_needle = results.profiles.r_needle;
r_syringe = results.profiles.r_syringe;
u_needle = results.profiles.u_needle;
u_syringe = results.profiles.u_syringe;
gamma_dot_needle = results.profiles.gamma_dot_needle;
gamma_dot_syringe = results.profiles.gamma_dot_syringe;
tau_needle = results.profiles.tau_needle;
tau_syringe = results.profiles.tau_syringe;

sample_name = results.sample_name;
safe_name = results.safe_name;
output_folder = results.output_folder;

%% Figure 1: Velocity profiles
f1 = figure('Position', [100, 100, 1200, 400]);

subplot(1,2,1)
plot(r_needle*1000, u_needle*1000, 'r-', 'LineWidth', 2);
xlabel('Radial Position (mm)');
ylabel('Velocity (mm/s)');
title(['Velocity Profile - Needle (', sample_name, ')']);
grid on;
xlim([0, Rn*1000]);

subplot(1,2,2)
plot(r_syringe*1000, u_syringe*1000, 'b-', 'LineWidth', 2);
xlabel('Radial Position (mm)');
ylabel('Velocity (mm/s)');
title(['Velocity Profile - Syringe (', sample_name, ')']);
grid on;
xlim([0, Rs*1000]);

sgtitle(['Power Law Velocity Profiles - ', sample_name]);

if save_figs
    exportgraphics(f1, fullfile(output_folder, sprintf('%s - Velocity Profile.png', safe_name)), 'Resolution', 1200);
end

%% Figure 2: Shear-rate profiles
f2 = figure('Position', [150, 150, 1200, 400]);

subplot(1,2,1)
plot(r_needle*1000, gamma_dot_needle, 'r-', 'LineWidth', 2);
xlabel('Radial Position (mm)');
ylabel('Shear Rate (s^-1)');
title(['Shear Rate Profile - Needle (', sample_name, ')']);
grid on;
xlim([0, Rn*1000]);

subplot(1,2,2)
plot(r_syringe*1000, gamma_dot_syringe, 'b-', 'LineWidth', 2);
xlabel('Radial Position (mm)');
ylabel('Shear Rate (s^-1)');
title(['Shear Rate Profile - Syringe (', sample_name, ')']);
grid on;
xlim([0, Rs*1000]);

sgtitle(['Shear Rate Distribution - ', sample_name]);

if save_figs
    exportgraphics(f2, fullfile(output_folder, sprintf('%s - Shear Rate Profile.png', safe_name)), 'Resolution', 1200);
end

%% Figure 3: Shear-stress profiles
f3 = figure('Position', [200, 200, 1200, 400]);

subplot(1,2,1)
plot(r_needle*1000, tau_needle, 'r-', 'LineWidth', 2);
xlabel('Radial Position (mm)');
ylabel('Shear Stress (Pa)');
title(['Shear Stress Profile - Needle (', sample_name, ')']);
grid on;
xlim([0, Rn*1000]);

subplot(1,2,2)
plot(r_syringe*1000, tau_syringe, 'b-', 'LineWidth', 2);
xlabel('Radial Position (mm)');
ylabel('Shear Stress (Pa)');
title(['Shear Stress Profile - Syringe (', sample_name, ')']);
grid on;
xlim([0, Rs*1000]);

sgtitle(['Shear Stress Distribution - ', sample_name]);

if save_figs
    exportgraphics(f3, fullfile(output_folder, sprintf('%s - Shear Stress Profiles.png', safe_name)), 'Resolution', 1200);
end

%% Figure 4: Summary
f4 = figure('Position', [250, 250, 1100, 650]);

subplot(2,2,1)
z1 = linspace(0, results.system_params.Ls, 101);
z2 = linspace(results.system_params.Ls, results.system_params.Ls + results.system_params.Ln, 101);

% --- NEW: Corrected continuous absolute pressure plotting ---
P1 = results.pressure.abs_syringe_inlet - results.pressure.abs_dpdz_syringe * z1;
P2 = results.pressure.abs_needle_inlet - results.pressure.abs_dpdz_needle * (z2 - results.system_params.Ls);

plot(z1*1000, P1, 'b-', 'LineWidth', 2); hold on;
plot(z2*1000, P2, 'r-', 'LineWidth', 2);
xlabel('Axial Position (mm)');
ylabel('Absolute Pressure (Pa)');
title(['Absolute Pressure Distribution - ', sample_name]);
legend('Syringe', 'Needle', 'Location', 'best');
grid on;

subplot(2,2,2)
bar([results.reynolds.Re_syringe, results.reynolds.Re_needle]); hold on;
yline(results.reynolds.Re_crit_estimated, 'k--', 'LineWidth', 2);
set(gca, 'XTickLabel', {'Syringe','Needle'});
ylabel('Generalized Reynolds Number');
title(['Flow Regime Estimate - ', sample_name]);
grid on;

subplot(2,2,3)
velocities = [results.flow.u_avg_syringe*1000, results.flow.u_max_syringe*1000, ...
              results.flow.u_avg_needle*1000, results.flow.u_max_needle*1000];
bar(velocities);
set(gca, 'XTickLabel', {'Syr Avg','Syr Max','Need Avg','Need Max'});
ylabel('Velocity (mm/s)');
title(['Velocity Comparison - ', sample_name]);
grid on;

subplot(2,2,4)
bar([results.shear.tau_wall_syringe, results.shear.tau_wall_needle]);
set(gca, 'XTickLabel', {'Syringe Wall','Needle Wall'});
ylabel('Wall Shear Stress (Pa)');
title(['Wall Shear Stress - ', sample_name]);
grid on;

sgtitle(['Bioprinting System Summary - ', sample_name]);

if save_figs
    exportgraphics(f4, fullfile(output_folder, sprintf('%s - System Overview.png', safe_name)), 'Resolution', 1200);
end

end

%% ========================================================================
function save_bioprinting_data_improved(results)
filename = fullfile(results.output_folder, [results.safe_name, '_data.txt']);
fid = fopen(filename, 'w');

if fid == -1
    warning('Could not create output file: %s', filename);
    return;
end

fprintf(fid, '=== BIOPRINTING ANALYSIS RESULTS ===\n');
fprintf(fid, 'Sample: %s\n\n', results.sample_name);

fprintf(fid, 'Model:\n');
fprintf(fid, '  Type: %s\n', results.model.type);
fprintf(fid, '  Assumptions:\n');
for i = 1:numel(results.model.assumptions)
    fprintf(fid, '    - %s\n', results.model.assumptions{i});
end
fprintf(fid, '\n');

fprintf(fid, 'System Parameters:\n');
fprintf(fid, '  Syringe radius (mm): %.6f\n', results.system_params.Rs*1000);
fprintf(fid, '  Needle radius (mm): %.6f\n', results.system_params.Rn*1000);
fprintf(fid, '  Needle length (mm): %.6f\n', results.system_params.Ln*1000);
fprintf(fid, '  Syringe length (mm): %.6f\n', results.system_params.Ls*1000);
fprintf(fid, '  Piston velocity (mm/s): %.6f\n', results.system_params.Vp*1000);
fprintf(fid, '  Density (kg/m^3): %.6f\n', results.system_params.rho);
fprintf(fid, '  Power-law index n: %.6f\n', results.system_params.n);
fprintf(fid, '  Consistency index K (Pa·s^n): %.6f\n\n', results.system_params.K);

fprintf(fid, 'Flow Analysis:\n');
fprintf(fid, '  Volumetric flow rate Q (m^3/s): %.9e\n', results.flow.Q);
fprintf(fid, '  Syringe avg velocity (m/s): %.9e\n', results.flow.u_avg_syringe);
fprintf(fid, '  Syringe max velocity (m/s): %.9e\n', results.flow.u_max_syringe);
fprintf(fid, '  Needle avg velocity (m/s): %.9e\n', results.flow.u_avg_needle);
fprintf(fid, '  Needle max velocity (m/s): %.9e\n\n', results.flow.u_max_needle);

fprintf(fid, 'Pressure Analysis:\n');
fprintf(fid, '  Syringe |dp/dz| (Pa/m): %.9e\n', results.pressure.abs_dpdz_syringe);
fprintf(fid, '  Needle  |dp/dz| (Pa/m): %.9e\n', results.pressure.abs_dpdz_needle);
fprintf(fid, '  Syringe pressure drop (Pa): %.9e\n', results.pressure.drop_syringe);
fprintf(fid, '  Needle pressure drop (Pa): %.9e\n', results.pressure.drop_needle);
fprintf(fid, '  Total viscous pressure drop (Pa): %.9e\n', results.pressure.drop_total_viscous);
fprintf(fid, '  Hydrostatic contribution (Pa): %.9e\n', results.pressure.drop_hydrostatic);
fprintf(fid, '  Total estimated system pressure drop (Pa): %.9e\n', results.pressure.drop_total_system);
% --- NEW: Absolute pressure logs ---
fprintf(fid, '  Absolute exit pressure (Pa): %.9e\n', results.pressure.abs_exit);
fprintf(fid, '  Absolute needle inlet pressure (Pa): %.9e\n', results.pressure.abs_needle_inlet);
fprintf(fid, '  Absolute syringe inlet pressure (Pa): %.9e\n\n', results.pressure.abs_syringe_inlet);

fprintf(fid, 'Shear Analysis:\n');
fprintf(fid, '  Syringe apparent wall shear rate (s^-1): %.9e\n', results.shear.gamma_dot_app_syringe);
fprintf(fid, '  Syringe true wall shear rate (s^-1): %.9e\n', results.shear.gamma_dot_true_syringe);
fprintf(fid, '  Syringe wall shear stress (Pa): %.9e\n', results.shear.tau_wall_syringe);
fprintf(fid, '  Needle apparent wall shear rate (s^-1): %.9e\n', results.shear.gamma_dot_app_needle);
fprintf(fid, '  Needle true wall shear rate (s^-1): %.9e\n', results.shear.gamma_dot_true_needle);
fprintf(fid, '  Needle wall shear stress (Pa): %.9e\n\n', results.shear.tau_wall_needle);

fprintf(fid, 'Reynolds Analysis:\n');
fprintf(fid, '  Syringe generalized Reynolds number: %.9e\n', results.reynolds.Re_syringe);
fprintf(fid, '  Needle generalized Reynolds number: %.9e\n', results.reynolds.Re_needle);
fprintf(fid, '  Critical Reynolds estimate: %.9e\n', results.reynolds.Re_crit_estimated);
fprintf(fid, '  Comment: %s\n\n', results.notes.reynolds_comment);

fprintf(fid, 'Radial Profiles - Needle:\n');
fprintf(fid, 'r_needle(m),u_needle(m/s),gamma_dot_needle(s^-1),tau_needle(Pa)\n');
for i = 1:length(results.profiles.r_needle)
    fprintf(fid, '%.9e,%.9e,%.9e,%.9e\n', ...
        results.profiles.r_needle(i), ...
        results.profiles.u_needle(i), ...
        results.profiles.gamma_dot_needle(i), ...
        results.profiles.tau_needle(i));
end

fprintf(fid, '\nRadial Profiles - Syringe:\n');
fprintf(fid, 'r_syringe(m),u_syringe(m/s),gamma_dot_syringe(s^-1),tau_syringe(Pa)\n');
for i = 1:length(results.profiles.r_syringe)
    fprintf(fid, '%.9e,%.9e,%.9e,%.9e\n', ...
        results.profiles.r_syringe(i), ...
        results.profiles.u_syringe(i), ...
        results.profiles.gamma_dot_syringe(i), ...
        results.profiles.tau_syringe(i));
end

fclose(fid);
fprintf('Data saved to file: %s\n', filename);

end