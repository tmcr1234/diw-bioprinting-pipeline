%RUN_SOLVER_V4  One-shot master driver for the full DIW extrusion sweep.
%
% Joins the legacy per-sample workflow (radial velocity / shear-rate /
% shear-stress profiles + System Overview PNGs + _data.txt files, produced
% by bioprinting_algorithm_3.m and bioprinting_algorithm_cross_v2.m) with
% the v3 slicer-lookup workflow (slicer_lookup_*.csv + 4-panel summary PNGs,
% produced by bioprinting_algorithm_v3.m) into a single execution.
%
% One run of this script reproduces, in one folder tree, every output that
% previously required ~14 manual reruns of run_solver_improved.m and
% run_solver_Cross_v2.m, plus the v3 outputs, plus a long-format master
% summary CSV.
%
% v4.1 (2026-05-21) iteration:
%   - Per-ink Ramp{1,2} fields are now OPTIONAL. To declare an ink that has
%     only no-pre-shear data, simply omit its Ramp2 field (or set it to []).
%     The driver skips the missing slice silently for every layer (legacy
%     per-sample PNGs, slicer CSVs, master summary).
%   - Each subset {ramp, Vp, needle, model} now also emits
%     all_samples_velocity_<orientation>.png -- a cross-ink overlay of the
%     needle and syringe radial velocity profiles, matching the legacy
%     create_multi_sample_plots output but at cfg.fig_dpi resolution and
%     with a deterministic (non-timestamped) filename.
%
% OUTPUT LAYOUT (under output_v4/)
%   <Ramp>/Vp_<v>_mmps/PL_<needle>/      legacy per-sample 4 PNGs + _data.txt
%                                        + all_samples_velocity_<orientation>.png
%   <Ramp>/Vp_<v>_mmps/Cross_<needle>/   same, for the Cross solver
%   slicer_lookup/<Ramp>/                v3 slicer CSV + 4-panel summary PNG
%   master_summary_v4.csv                long format, one row per
%                                        {ramp,ink,needle,Vp,model}
%
% Author: T.M.C. Rodrigues (PEMM/COPPE/UFRJ) -- 2026-05-21
%

clear; clc; close all;
addpath(fileparts(mfilename('fullpath')));

%% ================== CONFIGURATION ==================
cfg.output_root       = fullfile('output_v4');
cfg.orientation       = 'downward';      % 'horizontal' | 'upward' | 'downward'
cfg.include_hydrostat = true;
cfg.save_legacy_plots = true;            % per-sample 4 PNGs + _data.txt + all_samples PNG
cfg.save_slicer       = true;            % v3-style slicer CSV + summary PNG per (ramp,ink,needle)
cfg.fig_dpi           = 300;             % used for the all-samples overlay (legacy = 1200)
cfg.ramps_to_run      = {'Ramp1','Ramp2'};   % drop one to skip that ramp entirely
cfg.show_figures      = false;           % true = pop figures on screen during the run; false = headless

% Headless figure mode: every figure created during this script run is built
% off-screen and saved straight to PNG. Restored on exit (including errors
% and Ctrl-C) via onCleanup.
prev_fig_visible = get(groot, 'DefaultFigureVisible');
if ~cfg.show_figures
    set(groot, 'DefaultFigureVisible', 'off');
end
restore_fig_visible = onCleanup(@() set(groot, 'DefaultFigureVisible', prev_fig_visible)); %#ok<NASGU>

% Bench-realistic piston-velocity sweep (mm/s).
% Maps to head speeds of ~2-30 mm/s via Vp * (Rs/Rn)^2 amplification.
% Legacy values (3..25 mm/s) live on a mechanically unreachable operating
% point per CLAUDE.md memory v3_vs_old_simulation_audit.md -- they are NOT
% used here. To audit-reproduce the old runs, set legacy_Vp_mode = true.
cfg.legacy_Vp_mode = false;
if cfg.legacy_Vp_mode
    cfg.Vp_mm_s = [3 5 7 10 15 20 25];                     % mm/s, unreachable
else
    cfg.Vp_mm_s = [0.003 0.005 0.007 0.01 0.015 0.02 0.03 0.04];  % mm/s, bench
end

%% ================== SAMPLE DATABASE ==================
% Unified rheology + recovery + density per ink.
% Ramp 1 = no pre-shear (Fit_Muitos_Modelos_v4 / FitAll-Ramp1-v4.txt)
% Ramp 2 = post pre-shear at 200 1/s for 300 s (FitAll-Ramp2-v4.txt)
% Names follow CLAUDE.md Section 6.2 (canonical, not the v3 underscored aliases).
%
% RAMP OPTIONALITY:
%   To declare an ink that has NO Ramp 2 (post-pre-shear) data,
%   simply omit the .Ramp2 field, or set it to []:
%       inks(N).Ramp2 = [];
%   The driver will skip that ink in every Ramp 2 output without errors.

% --- C15 (CMC baseline) ---
inks(1).name     = 'C15';
inks(1).rho      = 898;
inks(1).Rrec_pct = 55;
inks(1).Ramp1    = struct('K_PL',92.5224, 'n_PL',0.397475, ...
                          'eta0',578.26,  'etaInf',0.0162966, ...
                          'lambda',3.71132,'m_Cross',0.810869);
inks(1).Ramp2    = struct('K_PL',39.4123, 'n_PL',0.387113, ...
                          'eta0',254.259, 'etaInf',0.001, ...
                          'lambda',3.91943,'m_Cross',0.81443);

% --- C15 Gira 5.5 (CMC + sunflower-oil nanoemulsion) ---
inks(2).name     = 'C15 Gira 5.5';
inks(2).rho      = 898;
inks(2).Rrec_pct = 40;
inks(2).Ramp1    = struct('K_PL',267.102, 'n_PL',0.24029, ...
                          'eta0',2439.5,  'etaInf',0.001, ...
                          'lambda',4.09905,'m_Cross',0.987578);
inks(2).Ramp2    = struct('K_PL',92.5405, 'n_PL',0.30149, ...
                          'eta0',951.474, 'etaInf',0.001, ...
                          'lambda',7.96028,'m_Cross',0.854576);

% --- Bozzano's Hair Gel (commercial reference) ---
inks(3).name     = 'Bozzano''s Hair Gel';
inks(3).rho      = 885;
inks(3).Rrec_pct = 96;
inks(3).Ramp1    = struct('K_PL',144.322, 'n_PL',0.234029, ...
                          'eta0',6354.55, 'etaInf',0.932091, ...
                          'lambda',73.0454,'m_Cross',0.865734);
inks(3).Ramp2    = struct('K_PL',141.494, 'n_PL',0.225653, ...
                          'eta0',57503.1, 'etaInf',0.816576, ...
                          'lambda',1544.3, 'm_Cross',0.825573);

%% ================== GEOMETRIES ==================
Rs  = 14.3e-3 / 2;
Ls  = 90e-3;

geom_21G = struct('Rs',Rs, 'R_n',0.515e-3/2, 'L_n',31.75e-3, ...
                  'Ls',Ls, 'label','21G', 'h_factor',0.7);
geom_22G = struct('Rs',Rs, 'R_n',0.413e-3/2, 'L_n',25.4e-3, ...
                  'Ls',Ls, 'label','22G', 'h_factor',0.7);
geoms = [geom_21G, geom_22G];

%% ================== PREP OUTPUT ROOT ==================
if ~exist(cfg.output_root, 'dir'), mkdir(cfg.output_root); end
fprintf('[v4] output_root = %s\n', cfg.output_root);
fprintf('[v4] orientation = %s | hydrostatic = %d\n', ...
        cfg.orientation, cfg.include_hydrostat);
fprintf('[v4] Vp sweep (mm/s) = [%s]\n', num2str(cfg.Vp_mm_s));
fprintf('[v4] ramps           = %s\n', strjoin(cfg.ramps_to_run, ', '));

t0 = tic;

%% ================== LEGACY-STYLE PER-SAMPLE LOOP ==================
% Loop order: ramp -> Vp -> needle -> ink, with PL and Cross collected in
% parallel cell arrays so the all-samples overlay can be emitted per
% {ramp, Vp, needle, model} subset after the ink loop closes.
master_rows = {};
for rr = 1:numel(cfg.ramps_to_run)
    ramp = cfg.ramps_to_run{rr};
    for vv = 1:numel(cfg.Vp_mm_s)
        Vp_mmps = cfg.Vp_mm_s(vv);
        Vp_mps  = Vp_mmps * 1e-3;
        % Folder tag uses the numeric Vp value in mm/s with '.' -> 'p':
        %   Vp = 10    -> 'Vp10mmps'
        %   Vp = 0.01  -> 'Vp0p01mmps'
        %   Vp = 0.003 -> 'Vp0p003mmps'
        vp_tag = sprintf('Vp%smmps', strrep(sprintf('%g',Vp_mmps),'.','p'));

        for jj = 1:numel(geoms)
            g = geoms(jj);

            pl_results    = {};   % collected across inks for this {ramp, Vp, needle}
            cross_results = {};

            folder_PL = fullfile(cfg.output_root, ramp, vp_tag, sprintf('PL_%s',    g.label));
            folder_C  = fullfile(cfg.output_root, ramp, vp_tag, sprintf('Cross_%s', g.label));

            for ii = 1:numel(inks)
                ink = inks(ii);
                if ~has_ramp(ink, ramp)
                    fprintf('[skip] %s has no %s data\n', ink.name, ramp);
                    continue;
                end
                rh = ink.(ramp);
                tag = sprintf('%s | %s | %s | %s', ramp, vp_tag, ink.name, g.label);

                % ----- Power Law (straight needle, legacy algorithm) -----
                if cfg.save_legacy_plots
                    fprintf('\n>>> [PL] %s  ->  %s\n', tag, folder_PL);
                    res_PL = bioprinting_algorithm_3( ...
                        g.Rs, g.R_n, g.L_n, g.Ls, Vp_mps, ...
                        rh.K_PL, rh.n_PL, ink.rho, ...
                        'Name',              ink.name, ...
                        'PlotResults',       true, ...
                        'SaveData',          true, ...
                        'SaveFigures',       true, ...
                        'Orientation',       cfg.orientation, ...
                        'OutputFolder',      folder_PL, ...
                        'IncludeHydrostatic',cfg.include_hydrostat);
                    pl_results{end+1} = res_PL;   %#ok<AGROW>
                    close all;   % prevent GUI flooding + RAM leak across iterations
                else
                    res_PL = [];
                end

                % ----- Cross (degenerate to straight cylinder, Rn_in = Rn_out) -----
                if cfg.save_legacy_plots
                    fprintf('\n>>> [Cross] %s  ->  %s\n', tag, folder_C);
                    res_C = bioprinting_algorithm_cross_v2( ...
                        g.Rs, g.R_n, g.R_n, g.L_n, g.Ls, Vp_mps, ...
                        rh.eta0, rh.etaInf, rh.lambda, rh.m_Cross, ink.rho, ...
                        'Name',              ink.name, ...
                        'PlotResults',       true, ...
                        'SaveData',          true, ...
                        'SaveFigures',       true, ...
                        'Orientation',       cfg.orientation, ...
                        'OutputFolder',      folder_C, ...
                        'IncludeHydrostatic',cfg.include_hydrostat);
                    cross_results{end+1} = res_C;   %#ok<AGROW>
                    close all;
                else
                    res_C = [];
                end

                % ----- Long-format master summary rows -----
                if ~isempty(res_PL)
                    master_rows(end+1,:) = { ramp, ink.name, g.label, Vp_mmps, ...
                        'PowerLaw', ...
                        res_PL.flow.Q*1e9, ...
                        res_PL.flow.u_avg_needle*1e3, res_PL.flow.u_max_needle*1e3, ...
                        res_PL.shear.tau_wall_needle, ...
                        res_PL.shear.gamma_dot_true_needle, ...
                        res_PL.pressure.drop_total_system/1e3, ...
                        res_PL.reynolds.Re_needle }; %#ok<AGROW>
                end
                if ~isempty(res_C)
                    % Field names differ from PL: Cross uses .u_avg_needle_exit
                    % (no u_max needle), .shear.gamma_dot_wall_needle (no
                    % Rabinowitsch-corrected "true" rate -- Cross integrates the
                    % discrete radial profile directly), and DOES populate
                    % .reynolds.Re_needle. See bioprinting_algorithm_cross_v2.m
                    % lines 145-159.
                    master_rows(end+1,:) = { ramp, ink.name, g.label, Vp_mmps, ...
                        'Cross', ...
                        res_C.flow.Q*1e9, ...
                        res_C.flow.u_avg_needle_exit*1e3, NaN, ...   % no u_max in Cross
                        res_C.shear.tau_wall_needle, ...
                        res_C.shear.gamma_dot_wall_needle, ...
                        res_C.pressure.drop_total_system/1e3, ...
                        res_C.reynolds.Re_needle }; %#ok<AGROW>
                end
            end

            % ----- Cross-ink overlay per subset (one per model) -----
            if cfg.save_legacy_plots
                overlay_tag_PL = sprintf('%s | %s | PL %s',    ramp, vp_tag, g.label);
                overlay_tag_C  = sprintf('%s | %s | Cross %s', ramp, vp_tag, g.label);
                create_all_samples_plot(pl_results,    folder_PL, cfg.orientation, ...
                                        overlay_tag_PL, cfg.fig_dpi);
                create_all_samples_plot(cross_results, folder_C,  cfg.orientation, ...
                                        overlay_tag_C,  cfg.fig_dpi);
            end
        end
    end
end

%% ================== V3 SLICER-LOOKUP LAYER ==================
% One call per (ramp, ink, needle) over the FULL Vp vector.
% Writes slicer_lookup_<ink>_<needle>.csv and 4-panel summary PNG.
% Inks without that ramp's parameters are skipped silently.
if cfg.save_slicer
    slicer_root = fullfile(cfg.output_root, 'slicer_lookup');
    if ~exist(slicer_root, 'dir'), mkdir(slicer_root); end

    for rr = 1:numel(cfg.ramps_to_run)
        ramp = cfg.ramps_to_run{rr};

        % Pre-filter inks that actually have this ramp
        ramp_inks_idx = [];
        for ii = 1:numel(inks)
            if has_ramp(inks(ii), ramp), ramp_inks_idx(end+1) = ii; end %#ok<AGROW>
        end
        if isempty(ramp_inks_idx)
            fprintf('[slicer] no inks have %s data; skipping that ramp\n', ramp);
            continue;
        end

        slicer_ramp = fullfile(slicer_root, ramp);
        if ~exist(slicer_ramp, 'dir'), mkdir(slicer_ramp); end

        for ii = ramp_inks_idx
            ink = inks(ii);
            rh  = ink.(ramp);
            v3sample = struct( ...
                'name',     ink.name, ...
                'K_PL',     rh.K_PL,    'n_PL',    rh.n_PL, ...
                'eta0',     rh.eta0,    'etaInf',  rh.etaInf, ...
                'lambda',   rh.lambda,  'm_Cross', rh.m_Cross, ...
                'Rrec_pct', ink.Rrec_pct);

            for jj = 1:numel(geoms)
                g = geoms(jj);
                fprintf('\n>>> [Slicer] %s | %s | %s\n', ramp, ink.name, g.label);
                bioprinting_algorithm_v3(v3sample, g, cfg.Vp_mm_s*1e-3, ...
                    'OutputFolder', slicer_ramp, ...
                    'SaveCSV',      true, ...
                    'PlotResults',  true);
            end
        end
    end
end

%% ================== MASTER SUMMARY CSV (LONG FORMAT) ==================
if ~isempty(master_rows)
    T = cell2table(master_rows, 'VariableNames', { ...
        'ramp','ink','needle','Vp_mm_s','model', ...
        'Q_mm3_s','u_avg_needle_mm_s','u_max_needle_mm_s', ...
        'tau_wall_needle_Pa','gamma_w_true_needle_invs', ...
        'dP_total_kPa','Re_needle'});
    master_csv = fullfile(cfg.output_root, 'master_summary_v4.csv');
    writetable(T, master_csv);
    fprintf('\n[v4] master summary -> %s  (%d rows)\n', master_csv, height(T));
end

fprintf('\n=== v4 DONE in %.1f s ===\n', toc(t0));
fprintf('Outputs under: %s/\n', cfg.output_root);

%% ================== LOCAL FUNCTIONS ==================

function tf = has_ramp(ink, ramp)
%HAS_RAMP  True if ink declares non-empty parameters for the given ramp.
    tf = isfield(ink, ramp) && ~isempty(ink.(ramp)) && isstruct(ink.(ramp));
end

function create_all_samples_plot(results_cell, output_folder, orientation, title_tag, fig_dpi)
%CREATE_ALL_SAMPLES_PLOT  Cross-ink overlay of needle + syringe velocity profiles
% for one {ramp, Vp, needle, model} subset.
%
% Mirrors the legacy create_multi_sample_plots nested function found in
% run_solver_improved.m and run_solver_Cross_v2.m, with three differences:
%   - drops samples whose profiles contain NaN/Inf (e.g., Cross collapse at
%     Q > Q_max), so a broken sample never aborts the whole overlay,
%   - uses a deterministic filename (overwrites on rerun) instead of a
%     timestamp suffix,
%   - figure stays Visible=off + close-after-export to avoid GUI flooding.

    % Filter out invalid profiles (Cross with NaN tau_w, etc.)
    valid = false(1, numel(results_cell));
    for k = 1:numel(results_cell)
        r = results_cell{k};
        if isstruct(r) && isfield(r, 'profiles') ...
                && isfield(r.profiles, 'u_needle') ...
                && isfield(r.profiles, 'u_syringe') ...
                && all(isfinite(r.profiles.u_needle)) ...
                && all(isfinite(r.profiles.u_syringe))
            valid(k) = true;
        end
    end
    results_cell = results_cell(valid);
    if isempty(results_cell)
        fprintf('[skip overlay] %s: no valid profiles\n', title_tag);
        return;
    end
    if ~exist(output_folder, 'dir'), mkdir(output_folder); end

    colors = lines(numel(results_cell));
    fig = figure('Visible','off','Position',[100 100 1200 400]);

    subplot(1,2,1); hold on;
    for k = 1:numel(results_cell)
        r = results_cell{k};
        plot(r.profiles.r_needle*1000, r.profiles.u_needle*1000, ...
             'LineWidth', 2, 'Color', colors(k,:), 'DisplayName', r.sample_name);
    end
    xlabel('Radial Position (mm)'); ylabel('Velocity (mm/s)');
    title('Needle Velocity Profiles'); grid on;
    legend('show', 'Interpreter','none', 'Location','best');

    subplot(1,2,2); hold on;
    for k = 1:numel(results_cell)
        r = results_cell{k};
        plot(r.profiles.r_syringe*1000, r.profiles.u_syringe*1000, ...
             'LineWidth', 2, 'Color', colors(k,:), 'DisplayName', r.sample_name);
    end
    xlabel('Radial Position (mm)'); ylabel('Velocity (mm/s)');
    title('Syringe Velocity Profiles'); grid on;
    legend('show', 'Interpreter','none', 'Location','best');

    sgtitle(sprintf('Velocity Profiles -- %s', title_tag), 'Interpreter','none');

    filename = fullfile(output_folder, sprintf('all_samples_velocity_%s.png', orientation));
    exportgraphics(fig, filename, 'Resolution', fig_dpi);
    close(fig);
    fprintf('[Saved overlay] %s\n', filename);
end
