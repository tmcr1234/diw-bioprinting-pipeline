clear; clc; close all;

% --- Define samples using Power Law constants and individual densities ---

% Unsheared samples
samples = struct(...
    'Name', {'C15', 'C15 Gira 5.5', 'Bozzano''s Hair Gel'}, ...
    'K',    {92.5224, 267.102, 144.322}, ...
    'n',    {0.397475, 0.24029, 0.234029}, ...
    'Vp',   0.025, ... % <--- This single scalar automatically applies to ALL elements!
    'rho',  {898, 898, 885} ...
);

% Sheared samples
samples_sheared = struct(...
    'Name', {'C15', 'C15 Gira 5.5', 'Bozzano''s Hair Gel'}, ...
    'K',    {39.4123, 92.5405, 141.494}, ...
    'n',    {0.387113, 0.30149, 0.225653}, ...
    'Vp',   0.025, ... % Applies to all
    'rho',  {898, 898, 885} ...
);


% --- Syringe & needle geometry ---
Rs = 14.3e-3/2;        % Syringe radius (m)
Rn = 0.000515/2;       % Needle radius (m)
Ln = 1.25 * 0.0254;    % Needle length (m)
Ls = 90e-3;            % Syringe length (m)



% --- User-defined options ---
orientation  = 'downward';   % 'horizontal', 'upward', or 'downward'
output_folder = fullfile(pwd, 'Vp25-Sheared');

% Create folder if it does not exist
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
results_all = cell(numel(samples_sheared),1);

for i = 1:numel(samples_sheared)
    s = samples_sheared(i); % <--- Changed {} to () here!
    
    results = bioprinting_algorithm_3( ...
        Rs, Rn, Ln, Ls, s.Vp, s.K, s.n, s.rho, ...
        'Name', s.Name, ...
        'PlotResults', true, ...
        'SaveData', true, ...
        'SaveFigures', true, ...
        'Orientation', orientation, ...
        'OutputFolder', output_folder, ...
        'IncludeHydrostatic', true);
        
    results_all{i} = results; % Curly braces are still correct here since results_all is a cell array
end

% Plot all samples together
create_multi_sample_plots(results_all, output_folder, orientation);

function create_multi_sample_plots(results_all, output_folder, orientation)
    colors = lines(numel(results_all));

    f = figure('Position',[100,100,1200,400]);

    % --- Velocity Profiles in Needle ---
    subplot(1,2,1);
    hold on;
    for i = 1:numel(results_all)
        r = results_all{i}.profiles.r_needle * 1000;
        u = results_all{i}.profiles.u_needle * 1000;
        plot(r, u, 'LineWidth', 2, 'Color', colors(i,:), ...
             'DisplayName', results_all{i}.sample_name);
    end
    xlabel('Radial Position (mm)');
    ylabel('Velocity (mm/s)');
    title('Needle Velocity Profiles');
    legend('show');
    grid on;

    % --- Velocity Profiles in Syringe ---
    subplot(1,2,2);
    hold on;
    for i = 1:numel(results_all)
        r = results_all{i}.profiles.r_syringe * 1000;
        u = results_all{i}.profiles.u_syringe * 1000;
        plot(r, u, 'LineWidth', 2, 'Color', colors(i,:), ...
             'DisplayName', results_all{i}.sample_name);
    end
    xlabel('Radial Position (mm)');
    ylabel('Velocity (mm/s)');
    title('Syringe Velocity Profiles');
    legend('show');
    grid on;

    sgtitle(sprintf('Velocity Profiles Comparison - Multiple Samples (%s)', orientation));

    % --- Automatically save figure into output folder ---
    filename = fullfile(output_folder, ...
        sprintf('all_samples_velocity_%s_%s.png', orientation, datestr(now,'yyyymmdd_HHMMSS')));

    exportgraphics(f, filename, 'Resolution', 1200);
    disp(['Figure saved as: ', filename]);
end