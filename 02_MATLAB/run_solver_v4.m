%RUN_SOLVER_V4  One-shot master driver for the full DIW extrusion sweep.
%
% Calls the unified superset solver bioprinting_algorithm_v4 once per
% (ramp x ink x needle) over the full piston-velocity vector. A single run
% produces, in one folder tree, every output that previously required the
% three separate solvers (bioprinting_algorithm_3, _cross_v2, _v3):
%
%   output_v4/<Ramp>/Vp_<v>_mmps/<ink>_<needle>_data.txt    combined PL+Cross
%   output_v4/<Ramp>/Vp_<v>_mmps/<ink>_<needle>_profiles.png radial profiles
%   output_v4/<Ramp>/Vp_<v>_mmps/<ink>_<needle>_system.png   pressure / Re / eta
%   output_v4/<Ramp>/slicer_lookup_<ink>_<needle>.csv        across-Vp slicer table
%   output_v4/<Ramp>/<ink>_<needle>_summary.png              across-Vp trends
%   output_v4/master_summary_v4.csv                          long format
%
% v4.2 (2026-06-01): replaced the legacy 3-algorithm orchestration with a
% single bioprinting_algorithm_v4 call. Geometry is straight needle only.
%
% Author: T.M.C. Rodrigues (PEMM/COPPE/UFRJ) -- 2026-06-01

clear; clc; close all;
addpath(fileparts(mfilename('fullpath')));

%% ================== CONFIGURATION ==================
cfg.output_root       = fullfile('output_v4');
cfg.orientation       = 'downward';      % 'horizontal' | 'upward' | 'downward'
cfg.include_hydrostat = true;
cfg.save_data         = true;            % per-Vp combined _data.txt
cfg.save_figures      = true;            % per-Vp + summary PNGs
cfg.save_csv          = true;            % across-Vp slicer CSV
cfg.num_points        = 200;             % radial discretisation
cfg.ramps_to_run      = {'Ramp1','Ramp2'};   % drop one to skip that ramp
cfg.show_figures      = false;           % false = headless

% Headless figure mode (restored on exit, incl. errors / Ctrl-C).
prev_fig_visible = get(groot, 'DefaultFigureVisible');
if ~cfg.show_figures
    set(groot, 'DefaultFigureVisible', 'off');
end
restore_fig_visible = onCleanup(@() set(groot, 'DefaultFigureVisible', prev_fig_visible)); %#ok<NASGU>

% Bench-realistic piston-velocity sweep (mm/s). Maps to head speeds of
% ~2-30 mm/s via Vp*(Rs/Rn)^2 amplification. Legacy 3..25 mm/s values live
% on a mechanically unreachable operating point (CLAUDE.md memory
% v3_vs_old_simulation_audit.md). Set legacy_Vp_mode = true to reproduce them.
cfg.legacy_Vp_mode = false;
if cfg.legacy_Vp_mode
    cfg.Vp_mm_s = [3 5 7 10 15 20 25];
else
    cfg.Vp_mm_s = [0.003 0.005 0.007 0.01 0.015 0.02 0.03 0.04];
end

%% ================== SAMPLE DATABASE (project-local) ==================
% The fitted ink parameters live OUTSIDE this shared file so Export/ can be a
% single shared clone symlinked across projects (see Export/CLAUDE.md ->
% "Shared Export / inks_local.m"). Each project supplies its own inks_local.m
% in its PROJECT ROOT (the cwd when you run this). Start from the template:
% Export/02_MATLAB/inks_local.template.m
%
% Contract for inks_local():
%   inks(i).name     char
%   inks(i).rho      kg/m^3
%   inks(i).Rrec_pct structural recovery at the deposition wall shear rate
%                    (~200 s^-1 at the nominal Vp=0.01 mm/s, 21G point;
%                     read from the project's Relatorio_Recovery.txt)
%   inks(i).Ramp1    struct(K_PL,n_PL,eta0,etaInf,lambda,m_Cross)  no pre-shear
%   inks(i).Ramp2    struct(...)  post pre-shear 200 1/s x 300 s; [] if absent
if exist(fullfile(pwd,'inks_local.m'),'file') ~= 2
    error('run_solver_v4:noInks', ...
        ['inks_local.m not found in the current folder:\n  %s\n' ...
         'cd to your PROJECT ROOT (the folder containing Export/) and create\n' ...
         'inks_local.m there. Copy Export/02_MATLAB/inks_local.template.m as a start.'], pwd);
end
inks = inks_local();

%% ================== GEOMETRIES ==================
Rs  = 14.3e-3 / 2;
Ls  = 90e-3;
geom_21G = struct('Rs',Rs,'R_n',0.515e-3/2,'L_n',31.75e-3,'Ls',Ls,'label','21G','h_factor',0.7);
geom_22G = struct('Rs',Rs,'R_n',0.413e-3/2,'L_n',25.4e-3,'Ls',Ls,'label','22G','h_factor',0.7);
geoms = [geom_21G, geom_22G];

%% ================== PREP OUTPUT ROOT ==================
if ~exist(cfg.output_root,'dir'), mkdir(cfg.output_root); end
fprintf('[v4] output_root = %s\n', cfg.output_root);
fprintf('[v4] orientation = %s | hydrostatic = %d\n', cfg.orientation, cfg.include_hydrostat);
fprintf('[v4] Vp sweep (mm/s) = [%s]\n', num2str(cfg.Vp_mm_s));
fprintf('[v4] ramps = %s\n', strjoin(cfg.ramps_to_run, ', '));
t0 = tic;

%% ================== MAIN LOOP ==================
master_rows = {};
for rr = 1:numel(cfg.ramps_to_run)
    ramp = cfg.ramps_to_run{rr};
    ramp_folder = fullfile(cfg.output_root, ramp);
    if ~exist(ramp_folder,'dir'), mkdir(ramp_folder); end

    for ii = 1:numel(inks)
        ink = inks(ii);
        if ~has_ramp(ink, ramp)
            fprintf('[skip] %s has no %s data\n', ink.name, ramp);
            continue;
        end
        rh = ink.(ramp);
        sample = struct('name',ink.name, ...
            'K_PL',rh.K_PL,'n_PL',rh.n_PL, ...
            'eta0',rh.eta0,'etaInf',rh.etaInf,'lambda',rh.lambda,'m_Cross',rh.m_Cross, ...
            'rho',ink.rho,'Rrec_pct',ink.Rrec_pct);

        for jj = 1:numel(geoms)
            g = geoms(jj);
            fprintf('\n>>> [v4] %s | %s | %s\n', ramp, ink.name, g.label);
            res = bioprinting_algorithm_v4(sample, g, cfg.Vp_mm_s*1e-3, ...
                'OutputFolder',       ramp_folder, ...
                'Orientation',        cfg.orientation, ...
                'IncludeHydrostatic', cfg.include_hydrostat, ...
                'SaveData',           cfg.save_data, ...
                'SaveFigures',        cfg.save_figures, ...
                'SaveCSV',            cfg.save_csv, ...
                'NumPoints',          cfg.num_points, ...
                'PlotResults',        true);

            % ----- long-format master summary rows (one per Vp x model) -----
            for vv = 1:numel(cfg.Vp_mm_s)
                % Machine geometry travels WITH every row. A reviewer flagged
                % the absence of these numbers from a manuscript as the single
                % biggest reproducibility blocker, even though every pressure
                % and velocity figure depends on them. They were always written
                % per-run to the _data.txt files; the summary CSV that tables
                % are actually built from dropped them. So did beta, which is
                % the heuristic die-swell closure - see the block above
                % beta_PL in bioprinting_algorithm_v4.m.
                geomcols = { res.geom.Rs*1e3, res.geom.R_n*1e3, res.geom.L_n*1e3, ...
                             res.geom.Ls*1e3, res.geom.h_factor, ink.rho };
                master_rows(end+1,:) = [ { ramp, ink.name, g.label, cfg.Vp_mm_s(vv), 'PowerLaw' }, ...
                    geomcols, ...
                    { res.PL.Q(vv)*1e9, res.PL.u_avg_needle(vv)*1e3, res.PL.u_max_needle(vv)*1e3, ...
                    res.PL.tau_w_needle(vv), res.PL.gamma_w_needle(vv), ...
                    res.PL.dP_total(vv)/1e3, res.PL.Re_needle(vv), ...
                    res.slicer.PL.beta, res.slicer.PL.w_line_mm, res.slicer.PL.h_layer_mm, ...
                    res.slicer.PL.v_print_mm_s(vv), res.slicer.PL.k_flow } ]; %#ok<AGROW>
                master_rows(end+1,:) = [ { ramp, ink.name, g.label, cfg.Vp_mm_s(vv), 'Cross' }, ...
                    geomcols, ...
                    { res.Cross.Q(vv)*1e9, res.Cross.u_avg_needle(vv)*1e3, res.Cross.u_max_needle(vv)*1e3, ...
                    res.Cross.tau_w_needle(vv), res.Cross.gamma_w_needle(vv), ...
                    res.Cross.dP_total(vv)/1e3, res.Cross.Re_needle(vv), ...
                    res.slicer.CR.beta, res.slicer.CR.w_line_mm, res.slicer.CR.h_layer_mm, ...
                    res.slicer.CR.v_print_mm_s(vv), res.slicer.CR.k_flow } ]; %#ok<AGROW>
            end
            close all;
        end
    end
end

%% ================== MASTER SUMMARY CSV ==================
if ~isempty(master_rows)
    T = cell2table(master_rows, 'VariableNames', { ...
        'ramp','ink','needle','Vp_mm_s','model', ...
        'Rs_mm','Rn_mm','Ln_mm','Ls_mm','h_factor','rho_kg_m3', ...
        'Q_mm3_s','u_avg_needle_mm_s','u_max_needle_mm_s', ...
        'tau_wall_needle_Pa','gamma_w_needle_invs','dP_total_kPa','Re_needle', ...
        'beta_swell','w_line_mm','h_layer_mm','v_print_mm_s','k_flow'});
    master_csv = fullfile(cfg.output_root, 'master_summary_v4.csv');
    writetable(T, master_csv);
    fprintf('\n[v4] master summary -> %s  (%d rows)\n', master_csv, height(T));
end

fprintf('\n=== v4 DONE in %.1f s ===\nOutputs under: %s/\n', toc(t0), cfg.output_root);

%% ================== LOCAL FUNCTIONS ==================
function tf = has_ramp(ink, ramp)
    tf = isfield(ink, ramp) && ~isempty(ink.(ramp)) && isstruct(ink.(ramp));
end
