function results = bioprinting_algorithm_cross_v2(Rs, Rn_in, Rn_out, Ln, Ls, Vp, eta0, etainf, lambda, m, rho, varargin)
%BIOPRINTING_ALGORITHM_CROSS
% Steady, laminar syringe-needle flow analysis for a Cross model fluid
% in extrusion-based bioprinting, accounting for tapered nozzle geometry.
%
% INPUTS
%   Rs      - Syringe internal radius (m)
%   Rn_in   - Needle/nozzle inlet radius (m)
%   Rn_out  - Needle/nozzle outlet radius (m)
%   Ln      - Needle length (m)
%   Ls      - Syringe length (m)
%   Vp      - Piston velocity (m/s)
%   eta0    - Zero-shear viscosity (Pa·s)
%   etainf  - Infinite-shear viscosity (Pa·s)
%   lambda  - Cross time constant (s)
%   m       - Cross rate constant / flow behavior index (dimensionless)
%   rho     - Fluid density (kg/m^3)
%
% OPTIONAL NAME-VALUE PAIRS
%   'PlotResults'       - true/false (default: true)
%   'NumPoints'         - number of spatial discretization points (default: 100)
%   'Name'              - sample name (default: "Cross_Fluid")
%   'SaveData'          - true/false for TXT export (default: true)
%   'SaveFigures'       - true/false for figure export (default: true)
%   'Orientation'       - 'horizontal', 'upward', or 'downward'
%   'OutputFolder'      - folder for exported files (default: pwd)
%   'IncludeHydrostatic'- true/false (default: true)

%% -------------------- INPUT PARSING --------------------
p = inputParser;
addRequired(p, 'Rs',  @(x) isnumeric(x) && x > 0);
addRequired(p, 'Rn_in', @(x) isnumeric(x) && x > 0);
addRequired(p, 'Rn_out', @(x) isnumeric(x) && x > 0);
addRequired(p, 'Ln',  @(x) isnumeric(x) && x > 0);
addRequired(p, 'Ls',  @(x) isnumeric(x) && x > 0);
addRequired(p, 'Vp',  @(x) isnumeric(x) && x > 0);
addRequired(p, 'eta0', @(x) isnumeric(x) && x > 0);
addRequired(p, 'etainf',@(x) isnumeric(x) && x >= 0);
addRequired(p, 'lambda',@(x) isnumeric(x) && x > 0);
addRequired(p, 'm',   @(x) isnumeric(x) && x > 0);
addRequired(p, 'rho', @(x) isnumeric(x) && x > 0);
addParameter(p, 'PlotResults', true, @islogical);
addParameter(p, 'NumPoints', 100, @isnumeric);
addParameter(p, 'Name', 'Cross_Fluid', @(x) ischar(x) || isstring(x));
addParameter(p, 'SaveData', true, @islogical);
addParameter(p, 'SaveFigures', true, @islogical);
addParameter(p, 'Orientation', 'horizontal', @(x) any(validatestring(x, {'horizontal','upward','downward'})));
addParameter(p, 'OutputFolder', pwd, @(x) ischar(x) || isstring(x));
addParameter(p, 'IncludeHydrostatic', true, @islogical);
parse(p, Rs, Rn_in, Rn_out, Ln, Ls, Vp, eta0, etainf, lambda, m, rho, varargin{:});

plot_results       = p.Results.PlotResults;
num_points         = p.Results.NumPoints;
sample_name        = char(p.Results.Name);
save_data          = p.Results.SaveData;
save_figures       = p.Results.SaveFigures;
orientation        = validatestring(p.Results.Orientation, {'horizontal','upward','downward'});
output_folder      = char(p.Results.OutputFolder);
include_hydrostatic = p.Results.IncludeHydrostatic;

if ~exist(output_folder, 'dir'), mkdir(output_folder); end
safe_name = regexprep(sample_name, '[<>:"/\\|?*]', '_');

%% -------------------- HEADER --------------------
fprintf('=== BIOPRINTING FLOW ANALYSIS - CROSS MODEL (TAPERED NOZZLE) ===\n\n');
fprintf('Sample: %s\n\n', sample_name);

%% -------------------- STEP 1: FLOW RATE --------------------
Q = pi * Rs^2 * Vp;
fprintf('STEP 1 - Flow Rate Calculation:\n');
fprintf('  Q = %.6f mm^3/s\n\n', Q*1e9);

%% -------------------- STEP 2: PRESSURE ANALYSIS (NUMERICAL) --------------------
fprintf('STEP 2 - Pressure Analysis (Numerical Integration):\n');
fprintf('  Solving Cross model equations... this may take a few seconds.\n');

% Syringe pressure drop (constant radius)
abs_dpdz_syringe = solve_dpdz_cross(Q, Rs, eta0, etainf, lambda, m, num_points);
pressure_drop_syringe = abs_dpdz_syringe * Ls;

% Tapered Nozzle pressure drop (variable radius)
z_needle = linspace(0, Ln, num_points);
R_z = Rn_in - ((Rn_in - Rn_out) / Ln) .* z_needle; % Linear taper profile
dpdz_needle_z = zeros(1, num_points);
for i = 1:num_points
    dpdz_needle_z(i) = solve_dpdz_cross(Q, R_z(i), eta0, etainf, lambda, m, 50);
end
pressure_drop_needle = trapz(z_needle, dpdz_needle_z);
abs_dpdz_needle_exit = dpdz_needle_z(end); % Max gradient is at the exit

% Hydrostatics
g = 9.81;
switch orientation
    case 'horizontal', orientation_factor = 0;
    case 'upward',     orientation_factor = +1;
    case 'downward',   orientation_factor = -1;
end
if include_hydrostatic
    hydrostatic_pressure_drop = orientation_factor * rho * g * (Ls + Ln);
else
    hydrostatic_pressure_drop = 0;
end

total_viscous_pressure_drop = pressure_drop_syringe + pressure_drop_needle;
total_system_pressure_drop = total_viscous_pressure_drop + hydrostatic_pressure_drop;

P_atm = 101325;
abs_pressure_exit = P_atm; 
abs_pressure_needle_inlet = abs_pressure_exit + pressure_drop_needle;
abs_pressure_syringe_inlet = abs_pressure_needle_inlet + pressure_drop_syringe;

fprintf('  Pressure drop across tapered needle:  %.2f Pa\n', pressure_drop_needle);
fprintf('  Pressure drop across syringe:         %.2f Pa\n', pressure_drop_syringe);
fprintf('  Total estimated system pressure drop: %.2f Pa\n', total_system_pressure_drop);
beta_contraction = (Rs / Rn_out)^2;
fprintf('\n  *** LOWER BOUND WARNING ***\n');
fprintf('  Fully-developed laminar flow assumed (no entrance losses, no wall slip).\n');
fprintf('  Contraction area ratio at syringe-needle junction: beta = %.0f\n', beta_contraction);
fprintf('  For beta >> 1, Bagley entrance correction may add 20-50%% to needle pressure.\n');
fprintf('  Report these pressures as lower bounds in the manuscript.\n\n');

%% -------------------- STEP 3 & 4: PROFILES & SHEAR ANALYSIS --------------------
[r_syringe, u_syringe, gamma_dot_syringe, tau_syringe] = calc_radial_profiles(abs_dpdz_syringe, Rs, eta0, etainf, lambda, m, num_points);
[r_needle, u_needle, gamma_dot_needle, tau_needle]     = calc_radial_profiles(abs_dpdz_needle_exit, Rn_out, eta0, etainf, lambda, m, num_points);

gamma_dot_wall_needle = max(gamma_dot_needle);
gamma_dot_wall_syringe = max(gamma_dot_syringe);
tau_wall_needle = max(tau_needle);
tau_wall_syringe = max(tau_syringe);

u_avg_needle_exit = Q / (pi * Rn_out^2);
u_avg_syringe = Q / (pi * Rs^2);

fprintf('STEP 3 & 4 - Shear and Velocity at Needle Exit (R = %.3f mm):\n', Rn_out*1000);
fprintf('  Wall shear rate:   %.2f s^-1\n', gamma_dot_wall_needle);
fprintf('  Wall shear stress: %.2f Pa\n', tau_wall_needle);
fprintf('  Max velocity:      %.2f mm/s\n\n', max(u_needle)*1000);

%% -------------------- STEP 5: REYNOLDS NUMBER --------------------
eta_wall_needle = tau_wall_needle / gamma_dot_wall_needle;
eta_wall_syringe = tau_wall_syringe / gamma_dot_wall_syringe;
Re_needle = (rho * u_avg_needle_exit * (2*Rn_out)) / eta_wall_needle;
Re_syringe = (rho * u_avg_syringe * (2*Rs)) / eta_wall_syringe;

%% -------------------- STORE RESULTS --------------------
results.sample_name = sample_name;
results.safe_name = safe_name;
results.output_folder = output_folder;
results.system_params = struct('Rs',Rs, 'Rn_in',Rn_in, 'Rn_out',Rn_out, 'Ln',Ln, 'Ls',Ls, 'Vp',Vp, 'rho',rho);
results.rheology = struct('eta0',eta0, 'etainf',etainf, 'lambda',lambda, 'm',m);
results.flow.Q = Q;
results.flow.u_avg_needle_exit = u_avg_needle_exit;
results.flow.u_avg_syringe = u_avg_syringe;
results.pressure = struct('drop_needle',pressure_drop_needle, 'drop_syringe',pressure_drop_syringe, ...
    'abs_exit',abs_pressure_exit, 'abs_needle_inlet',abs_pressure_needle_inlet, 'abs_syringe_inlet',abs_pressure_syringe_inlet, ...
    'drop_total_system',total_system_pressure_drop, 'drop_hydrostatic',hydrostatic_pressure_drop);
results.shear = struct('gamma_dot_wall_needle',gamma_dot_wall_needle, 'gamma_dot_wall_syringe',gamma_dot_wall_syringe, ...
    'tau_wall_needle',tau_wall_needle, 'tau_wall_syringe',tau_wall_syringe);
results.reynolds = struct('Re_needle',Re_needle, 'Re_syringe',Re_syringe);
results.profiles = struct('r_needle',r_needle, 'u_needle',u_needle, 'gamma_dot_needle',gamma_dot_needle, 'tau_needle',tau_needle, ...
    'r_syringe',r_syringe, 'u_syringe',u_syringe, 'gamma_dot_syringe',gamma_dot_syringe, 'tau_syringe',tau_syringe, ...
    'z_needle',z_needle, 'dpdz_needle_z',dpdz_needle_z);

%% -------------------- PLOTS --------------------
if plot_results
    f = figure('Position', [100, 100, 1200, 800]);
    
    % Velocity
    subplot(2,2,1);
    plot(r_needle*1000, u_needle*1000, 'r-', 'LineWidth', 2); hold on;
    plot(r_syringe*1000, u_syringe*1000, 'b--', 'LineWidth', 2);
    xlabel('Radial Position (mm)'); ylabel('Velocity (mm/s)');
    title('Velocity Profiles'); legend('Needle Exit', 'Syringe'); grid on;
    
    % Shear Rate
    subplot(2,2,2);
    plot(r_needle*1000, gamma_dot_needle, 'r-', 'LineWidth', 2); hold on;
    plot(r_syringe*1000, gamma_dot_syringe, 'b--', 'LineWidth', 2);
    xlabel('Radial Position (mm)'); ylabel('Shear Rate (s^{-1})');
    title('Shear Rate Distribution'); grid on;
    
    % Absolute Pressure (Tapered Integration)
    subplot(2,2,3);
    z1 = linspace(0, Ls, 100);
    P1 = abs_pressure_syringe_inlet - abs_dpdz_syringe * z1;
    
    P2 = zeros(1, num_points);
    P2(1) = abs_pressure_needle_inlet;
    for i = 2:num_points
        P2(i) = abs_pressure_needle_inlet - trapz(z_needle(1:i), dpdz_needle_z(1:i));
    end
    z2 = z_needle + Ls;
    
    plot(z1*1000, P1/1000, 'b-', 'LineWidth', 2); hold on;
    plot(z2*1000, P2/1000, 'r-', 'LineWidth', 2);
    xlabel('Axial Position (mm)'); ylabel('Absolute Pressure (kPa)');
    title('System Pressure Profile'); legend('Syringe', 'Tapered Needle'); grid on;
    
    % Apparent Viscosity Profile
    subplot(2,2,4);
    eta_app_needle = tau_needle ./ max(gamma_dot_needle, 1e-10);
    plot(r_needle*1000, eta_app_needle, 'g-', 'LineWidth', 2);
    xlabel('Radial Position (mm)'); ylabel('\eta_{app} (Pa\cdots)');
    title('Viscosity across Needle Exit'); grid on;
    
    sgtitle(['Cross Model Analysis: ', sample_name]);
    if save_figures
        exportgraphics(f, fullfile(output_folder, sprintf('%s_Overview.png', safe_name)), 'Resolution', 600);
    end
end

%% -------------------- SAVE DATA --------------------
if save_data
    save_bioprinting_data_cross(results);
end

fprintf('=== ANALYSIS COMPLETE ===\n');
end

%% ================= HELPER FUNCTIONS =================
function dpdz = solve_dpdz_cross(Q_target, R, eta0, etainf, lambda, m, N)
    % Guaranteed bracket: at dpdz_lo the Cross fluid flows LESS than Q_target
    % (eta~eta0 so high resistance); at dpdz_hi it flows MORE (eta~etainf so low resistance).
    dpdz_lo = (8 * etainf * Q_target) / (pi * R^4);
    dpdz_hi = (8 * eta0   * Q_target) / (pi * R^4);
    options = optimset('Display','off', 'TolX', 1e-4);
    dpdz = fzero(@(dp) calc_Q_discrete(dp, R, eta0, etainf, lambda, m, N) - Q_target, ...
                 [dpdz_lo, dpdz_hi], options);
end

function Q_calc = calc_Q_discrete(dpdz, R, eta0, etainf, lambda, m, N)
    r = linspace(0, R, N);
    tau_r = (dpdz * r) / 2;
    g_r = zeros(size(r));
    for i = 2:N
        % g_lo~0: eta(g)*g ~ 0 < tau  ->  objective negative
        % g_hi=tau/etainf: eta(g_hi)*g_hi > etainf*(tau/etainf) = tau  ->  objective positive
        g_lo = 1e-10;
        g_hi = tau_r(i) / etainf;
        obj = @(g) (etainf + (eta0 - etainf)./(1 + (lambda*abs(g)).^m)).*g - tau_r(i);
        g_r(i) = fzero(obj, [g_lo, g_hi]);
    end
    Q_calc = pi * trapz(r, r.^2 .* g_r);
end

function [r, u, g_r, tau_r] = calc_radial_profiles(dpdz, R, eta0, etainf, lambda, m, N)
    r = linspace(0, R, N);
    tau_r = (dpdz * r) / 2;
    g_r = zeros(size(r));
    for i = 2:N
        g_lo = 1e-10;
        g_hi = tau_r(i) / etainf;
        obj = @(g) (etainf + (eta0 - etainf)./(1 + (lambda*abs(g)).^m)).*g - tau_r(i);
        g_r(i) = fzero(obj, [g_lo, g_hi]);
    end
    u = zeros(size(r));
    for i = 1:N-1
        u(i) = trapz(r(i:end), g_r(i:end));
    end
end

function save_bioprinting_data_cross(results)
    filename = fullfile(results.output_folder, [results.safe_name, '_data.txt']);
    fid = fopen(filename, 'w');

    if fid == -1
        warning('Could not create output file: %s', filename);
        return;
    end

    fprintf(fid, '=== BIOPRINTING ANALYSIS RESULTS ===\n');
    fprintf(fid, 'Sample: %s\n\n', results.sample_name);

    fprintf(fid, 'Model:\n');
    fprintf(fid, '  Type: Cross Model (Tapered Nozzle)\n\n');

    fprintf(fid, 'System Parameters:\n');
    fprintf(fid, '  Syringe internal radius (mm): %.6f\n', results.system_params.Rs*1000);
    fprintf(fid, '  Needle inlet radius (mm): %.6f\n', results.system_params.Rn_in*1000);
    fprintf(fid, '  Needle outlet radius (mm): %.6f\n', results.system_params.Rn_out*1000);
    fprintf(fid, '  Needle length (mm): %.6f\n', results.system_params.Ln*1000);
    fprintf(fid, '  Syringe length (mm): %.6f\n', results.system_params.Ls*1000);
    fprintf(fid, '  Piston velocity (mm/s): %.6f\n', results.system_params.Vp*1000);
    fprintf(fid, '  Density (kg/m^3): %.6f\n\n', results.system_params.rho);

    fprintf(fid, 'Rheology Parameters (Cross Model):\n');
    fprintf(fid, '  Zero-shear viscosity eta0 (Pa.s): %.6f\n', results.rheology.eta0);
    fprintf(fid, '  Infinite-shear viscosity etainf (Pa.s): %.6f\n', results.rheology.etainf);
    fprintf(fid, '  Cross time constant lambda (s): %.6f\n', results.rheology.lambda);
    fprintf(fid, '  Cross rate constant m: %.6f\n\n', results.rheology.m);

    fprintf(fid, 'Flow Analysis:\n');
    fprintf(fid, '  Volumetric flow rate Q (m^3/s): %.9e\n', results.flow.Q);
    fprintf(fid, '  Syringe avg velocity (m/s): %.9e\n', results.flow.u_avg_syringe);
    fprintf(fid, '  Needle exit avg velocity (m/s): %.9e\n\n', results.flow.u_avg_needle_exit);

    fprintf(fid, 'Pressure Analysis (LOWER BOUNDS -- fully-developed flow assumed):\n');
    fprintf(fid, '  *** Entrance losses not modelled. Contraction beta = (Rs/Rn_out)^2 = %.0f ***\n', ...
            (results.system_params.Rs / results.system_params.Rn_out)^2);
    fprintf(fid, '  *** For beta >> 1, Bagley correction may add 20-50%% to needle pressure. ***\n');
    fprintf(fid, '  Syringe pressure drop (Pa): %.9e\n', results.pressure.drop_syringe);
    fprintf(fid, '  Needle pressure drop (Pa): %.9e\n', results.pressure.drop_needle);
    fprintf(fid, '  Hydrostatic contribution (Pa): %.9e\n', results.pressure.drop_hydrostatic);
    fprintf(fid, '  Total estimated system pressure drop (Pa): %.9e\n', results.pressure.drop_total_system);
    fprintf(fid, '  Absolute exit pressure (Pa): %.9e\n', results.pressure.abs_exit);
    fprintf(fid, '  Absolute needle inlet pressure (Pa): %.9e\n', results.pressure.abs_needle_inlet);
    fprintf(fid, '  Absolute syringe inlet pressure (Pa): %.9e\n\n', results.pressure.abs_syringe_inlet);

    fprintf(fid, 'Shear Analysis (Max Wall Values):\n');
    fprintf(fid, '  Syringe wall shear rate (s^-1): %.9e\n', results.shear.gamma_dot_wall_syringe);
    fprintf(fid, '  Syringe wall shear stress (Pa): %.9e\n', results.shear.tau_wall_syringe);
    fprintf(fid, '  Needle exit wall shear rate (s^-1): %.9e\n', results.shear.gamma_dot_wall_needle);
    fprintf(fid, '  Needle exit wall shear stress (Pa): %.9e\n\n', results.shear.tau_wall_needle);

    fprintf(fid, 'Reynolds Analysis:\n');
    fprintf(fid, '  Syringe Reynolds number: %.9e\n', results.reynolds.Re_syringe);
    fprintf(fid, '  Needle exit Reynolds number: %.9e\n\n', results.reynolds.Re_needle);

    fprintf(fid, 'Radial Profiles - Needle Exit:\n');
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