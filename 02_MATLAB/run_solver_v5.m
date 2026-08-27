%RUN_SOLVER_V5  Driver for the anchored-closure solver (Tanner die swell).
%
% Same sweep as run_solver_v4 — one call per (ramp x ink x needle) over the
% piston-velocity vector — but through bioprinting_algorithm_v5, which
% replaces v4's two unsourced heuristics with:
%
%   die swell        Tanner (1970) closure keyed on the MEASURED N1
%   critical Reynolds Ryan & Johnson (1959)
%
% WHAT IS AND IS NOT DUPLICATED FROM run_solver_v4
% ------------------------------------------------
% Only the ORCHESTRATION loop is repeated here. No physics, no pressure
% decomposition and — importantly — no slicer formula lives in this file.
% w_line, v_print and k_flow are computed once, inside the solver, exactly
% as they are for v4. If you find yourself about to paste a formula into
% this file, you are in the wrong file.
%
% WHERE N1 COMES FROM
% -------------------
% Each ink supplies inks(i).N1_wall_Pa — the first normal-stress difference
% at the NEEDLE wall shear rate, in Pa. Produce it with:
%
%   python3 Export/01_Python/extract_N1_tanner.py \
%       --force-folder "./Reologia/Flow with force" \
%       --master-summary output_v4/master_summary_v4.csv --needle 21G --vp 0.01
%
% which reads a force-augmented export, tares the transducer, excludes any
% edge-fracture window, interpolates N1 to the wall shear rate the v4 solver
% actually predicts, and prints a paste-ready inks_local.m line.
%
% An ink with no N1_wall_Pa is NOT an error. It runs on the v4 heuristic,
% is tagged 'v5-fallback-v4-heuristic' in every output, and is listed in the
% closure report at the end of the run. Mixing anchored and heuristic inks
% in one sweep is allowed precisely so a partial measurement campaign is
% usable — but a manuscript must not quote the two as if they were the same.
%
%   output_v5/<Ramp>/...                        as v4
%   output_v5/master_summary_v5.csv             + closure, N1, beta columns
%   output_v5/closure_report_v5.txt             per-ink disclosure block
%
% Author: T.M.C. Rodrigues (PEMM/COPPE/UFRJ) -- 2026-08-26

clear; clc; close all;
addpath(fileparts(mfilename('fullpath')));

%% ================== CONFIGURATION ==================
cfg.output_root       = fullfile('output_v5');
cfg.orientation       = 'downward';
cfg.include_hydrostat = true;
cfg.save_data         = true;
cfg.save_figures      = true;
cfg.save_csv          = true;
cfg.num_points        = 200;
cfg.ramps_to_run      = {'Ramp1','Ramp2'};
cfg.show_figures      = false;
% Operating point whose wall shear stress anchors the Tanner ratio. beta is
% one number per run, so one point has to be chosen; make it the point you
% actually print at, and make it the SAME one you gave extract_N1_tanner.py.
cfg.beta_ref_Vp_mm_s  = 0.01;

prev_fig_visible = get(groot, 'DefaultFigureVisible');
if ~cfg.show_figures
    set(groot, 'DefaultFigureVisible', 'off');
end
restore_fig_visible = onCleanup(@() set(groot, 'DefaultFigureVisible', prev_fig_visible)); %#ok<NASGU>

cfg.legacy_Vp_mode = false;
if cfg.legacy_Vp_mode
    cfg.Vp_mm_s = [3 5 7 10 15 20 25];
else
    cfg.Vp_mm_s = [0.003 0.005 0.007 0.01 0.015 0.02 0.03 0.04];
end

%% ================== SAMPLE DATABASE (project-local) ==================
% Contract as run_solver_v4, PLUS the optional field:
%   inks(i).N1_wall_Pa   first normal-stress difference at the needle wall
%                        shear rate (Pa). Omit or leave [] to fall back.
if exist(fullfile(pwd,'inks_local.m'),'file') ~= 2
    error('run_solver_v5:noInks', ...
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

%% ================== PREP ==================
if ~exist(cfg.output_root,'dir'), mkdir(cfg.output_root); end
fprintf('[v5] output_root = %s\n', cfg.output_root);
fprintf('[v5] beta anchored at Vp = %g mm/s\n', cfg.beta_ref_Vp_mm_s);

n_anchored = 0; n_fallback = 0;
for ii = 1:numel(inks)
    if isfield(inks(ii),'N1_wall_Pa') && ~isempty(inks(ii).N1_wall_Pa) ...
            && all(isfinite(inks(ii).N1_wall_Pa))
        fprintf('[v5]   %-18s N1_wall = %.4g Pa  -> Tanner\n', inks(ii).name, inks(ii).N1_wall_Pa);
        n_anchored = n_anchored + 1;
    else
        fprintf('[v5]   %-18s no N1     -> v4 HEURISTIC fallback\n', inks(ii).name);
        n_fallback = n_fallback + 1;
    end
end
if n_fallback > 0
    fprintf(['[v5] %d of %d ink(s) have no measured N1 and will run on the\n' ...
             '     unsourced v4 die-swell closure. Their slicer numbers are not\n' ...
             '     comparable with the anchored ones - see closure_report_v5.txt.\n'], ...
             n_fallback, numel(inks));
end
t0 = tic;

%% ================== MAIN LOOP ==================
master_rows = {};
closure_rows = {};
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
        if isfield(ink,'N1_wall_Pa') && ~isempty(ink.N1_wall_Pa)
            sample.N1_wall_Pa = ink.N1_wall_Pa;
        end

        for jj = 1:numel(geoms)
            g = geoms(jj);
            fprintf('\n>>> [v5] %s | %s | %s\n', ramp, ink.name, g.label);
            res = bioprinting_algorithm_v5(sample, g, cfg.Vp_mm_s*1e-3, ...
                'OutputFolder',       ramp_folder, ...
                'Orientation',        cfg.orientation, ...
                'IncludeHydrostatic', cfg.include_hydrostat, ...
                'SaveData',           cfg.save_data, ...
                'SaveFigures',        cfg.save_figures, ...
                'SaveCSV',            cfg.save_csv, ...
                'NumPoints',          cfg.num_points, ...
                'PlotResults',        true, ...
                'BetaRefVp_mm_s',     cfg.beta_ref_Vp_mm_s);

            br = res.beta_report;
            closure_rows(end+1,:) = { ramp, ink.name, g.label, res.closure, ...
                br.N1_wall_Pa, br.tau_wall_PL_Pa, ...
                br.beta_heuristic_PL, br.beta_used_PL, ...
                br.k_flow_heuristic_PL, br.k_flow_used_PL, ...
                br.extrusion_multiplier_heuristic_PL, ...
                br.extrusion_multiplier_used_PL }; %#ok<AGROW>

            for vv = 1:numel(cfg.Vp_mm_s)
                geomcols = { res.geom.Rs*1e3, res.geom.R_n*1e3, res.geom.L_n*1e3, ...
                             res.geom.Ls*1e3, res.geom.h_factor, ink.rho };
                master_rows(end+1,:) = [ { ramp, ink.name, g.label, cfg.Vp_mm_s(vv), 'PowerLaw', res.closure }, ...
                    geomcols, ...
                    { res.PL.Q(vv)*1e9, res.PL.u_avg_needle(vv)*1e3, res.PL.u_max_needle(vv)*1e3, ...
                    res.PL.tau_w_needle(vv), res.PL.gamma_w_needle(vv), ...
                    res.PL.dP_total(vv)/1e3, res.PL.Re_needle(vv), res.Re_crit, ...
                    br.N1_wall_Pa, res.slicer.PL.beta, res.slicer.PL.w_line_mm, ...
                    res.slicer.PL.h_layer_mm, res.slicer.PL.v_print_mm_s(vv), ...
                    res.slicer.PL.k_flow } ]; %#ok<AGROW>
                master_rows(end+1,:) = [ { ramp, ink.name, g.label, cfg.Vp_mm_s(vv), 'Cross', res.closure }, ...
                    geomcols, ...
                    { res.Cross.Q(vv)*1e9, res.Cross.u_avg_needle(vv)*1e3, res.Cross.u_max_needle(vv)*1e3, ...
                    res.Cross.tau_w_needle(vv), res.Cross.gamma_w_needle(vv), ...
                    res.Cross.dP_total(vv)/1e3, res.Cross.Re_needle(vv), res.Re_crit, ...
                    br.N1_wall_Pa, res.slicer.CR.beta, res.slicer.CR.w_line_mm, ...
                    res.slicer.CR.h_layer_mm, res.slicer.CR.v_print_mm_s(vv), ...
                    res.slicer.CR.k_flow } ]; %#ok<AGROW>
            end
            close all;
        end
    end
end

%% ================== MASTER SUMMARY CSV ==================
if ~isempty(master_rows)
    T = cell2table(master_rows, 'VariableNames', { ...
        'ramp','ink','needle','Vp_mm_s','model','closure', ...
        'Rs_mm','Rn_mm','Ln_mm','Ls_mm','h_factor','rho_kg_m3', ...
        'Q_mm3_s','u_avg_needle_mm_s','u_max_needle_mm_s', ...
        'tau_wall_needle_Pa','gamma_w_needle_invs','dP_total_kPa', ...
        'Re_needle','Re_crit','N1_wall_Pa', ...
        'beta_swell','w_line_mm','h_layer_mm','v_print_mm_s','k_flow'});
    master_csv = fullfile(cfg.output_root, 'master_summary_v5.csv');
    writetable(T, master_csv);
    fprintf('\n[v5] master summary -> %s  (%d rows)\n', master_csv, height(T));
end

%% ================== CLOSURE DISCLOSURE REPORT ==================
% Generated, not remembered. This is the block a Methods section needs in
% order to state which die-swell closure produced which slicer number.
if ~isempty(closure_rows)
    rep = fullfile(cfg.output_root, 'closure_report_v5.txt');
    fid = fopen(rep, 'w');
    fprintf(fid, '=== DIE-SWELL CLOSURE REPORT (bioprinting_algorithm_v5) ===\n\n');
    fprintf(fid, 'beta anchored at the operating point Vp = %g mm/s.\n\n', cfg.beta_ref_Vp_mm_s);
    fprintf(fid, 'Closures:\n');
    fprintf(fid, '  v5-tanner                  Tanner 1970 (doi:10.1002/pol.1970.160081203)\n');
    fprintf(fid, '                             keyed on the MEASURED N1 at the needle wall.\n');
    fprintf(fid, '                             N1 -> 0 gives the 0.13 inelastic floor of\n');
    fprintf(fid, '                             Nickell, Tanner & Yamada 1974\n');
    fprintf(fid, '                             (doi:10.1017/S0022112074001339).\n');
    fprintf(fid, '  v5-fallback-v4-heuristic   beta = 0.30*(1-n). NO PUBLISHED SOURCE.\n');
    fprintf(fid, '                             Reported only because no N1 was supplied.\n\n');
    fprintf(fid, ['NOTE ON THE ESTIMATE. Tanner describes a FREE JET leaving a long\n' ...
                  'capillary. A DIW road is deposited against a substrate and is still\n' ...
                  'swelling when it lands, so a measured road width is expected to\n' ...
                  'exceed this beta. Report both; the gap is a definition mismatch,\n' ...
                  'not a fitting error.\n\n']);
    fprintf(fid, ['NOTE ON N1. The min-baseline tare treats the quietest point of the\n' ...
                  'ramp as zero normal stress, so N1 is short by N1(gamma_dot_min) -\n' ...
                  'a downward, conservative bias (1.6%% on a synthetic check).\n\n']);
    fprintf(fid, '%-8s %-18s %-7s %-26s %10s %10s %9s %9s %10s %10s\n', ...
        'ramp','ink','needle','closure','N1_w(Pa)','tau_w(Pa)', ...
        'beta_v4','beta_used','EM_v4','EM_used');
    fprintf(fid, '%s\n', repmat('-', 1, 132));
    for k = 1:size(closure_rows,1)
        c = closure_rows(k,:);
        fprintf(fid, '%-8s %-18s %-7s %-26s %10.4g %10.4g %9.4f %9.4f %10.4f %10.4f\n', ...
            c{1}, c{2}, c{3}, c{4}, c{5}, c{6}, c{7}, c{8}, c{11}, c{12});
    end
    fclose(fid);
    fprintf('[v5] closure report -> %s\n', rep);
    type(rep);
end

fprintf('\n=== v5 DONE in %.1f s ===\nOutputs under: %s/\n', toc(t0), cfg.output_root);

%% ================== LOCAL FUNCTIONS ==================
function tf = has_ramp(ink, ramp)
    tf = isfield(ink, ramp) && ~isempty(ink.(ramp)) && isstruct(ink.(ramp));
end
