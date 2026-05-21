%RUN_SOLVER_V3  Unified PL + Cross driver for all inks and needle sizes.
%
% Produces, per ink x needle:
%   - slicer_lookup_<sample>_<needle>.csv  (PL + Cross side-by-side)
%   - plots_<sample>_<needle>.png          (4-panel summary)
%
% After running, hand the CSVs to an intern using the SOP at
% Latex/Printing_Parameters_SOP_v1.tex
%
% Author: T.M.C. Rodrigues -- 2026-05-14

clear; clc; close all;

addpath(fileparts(mfilename('fullpath')));

%% -------------------- SAMPLE DATABASE (Ramp 1, no pre-shear) --------------------
% Values from Fit_Muitos_Modelos_v4 (FitAll-Ramp1-v4.txt) + Recovery_v1
samples(1).name     = 'C15';
samples(1).K_PL     = 92.52;
samples(1).n_PL     = 0.397;
samples(1).eta0     = 578.3;
samples(1).etaInf   = 0.0163;
samples(1).lambda   = 3.711;
samples(1).m_Cross  = 0.811;
samples(1).Rrec_pct = 55;    % at gamma_dot ~ 150 1/s (between 100 and 150 Hz recovery)

samples(2).name     = 'C15_SF_5_5';
samples(2).K_PL     = 267.1;
samples(2).n_PL     = 0.240;
samples(2).eta0     = 2439.5;
samples(2).etaInf   = 0.001;
samples(2).lambda   = 4.099;
samples(2).m_Cross  = 0.988;
samples(2).Rrec_pct = 40;    % at ~150 Hz (noisy band; conservative midpoint)

samples(3).name     = 'Hair Gel';
samples(3).K_PL     = 144.3;
samples(3).n_PL     = 0.234;
samples(3).eta0     = 6354.6;
samples(3).etaInf   = 0.932;
samples(3).lambda   = 73.05;
samples(3).m_Cross  = 0.866;
samples(3).Rrec_pct = 96;    % near-complete recovery, robust at all shear

%% -------------------- GEOMETRIES --------------------
% Common syringe: BD 10 mL (Rs = 14.3/2 mm)
Rs = 14.3e-3 / 2;

geom_21G.Rs       = Rs;
geom_21G.R_n      = 0.515e-3 / 2;
geom_21G.L_n      = 31.75e-3;
geom_21G.label    = '21G';
geom_21G.h_factor = 0.7;

geom_22G.Rs       = Rs;
geom_22G.R_n      = 0.413e-3 / 2;
geom_22G.L_n      = 25.4e-3;
geom_22G.label    = '22G';
geom_22G.h_factor = 0.7;

geoms = [geom_21G, geom_22G];

%% -------------------- PISTON VELOCITIES --------------------
% Bench-realistic sweep: v_print = (pi/2.8) * (Rs/Rn)^2 / (1+beta) * Vp
% For 21G the amplification is ~700-865x (ink-dependent via beta_swell),
% so Vp = 0.003-0.04 mm/s maps to head speeds of ~2-30 mm/s,
% which is the lab DIW printer's working envelope. Values above 0.05 mm/s
% imply head speeds >35 mm/s and are mechanically unreachable on this rig.
Vp_mm_s = [0.003 0.005 0.007 0.01 0.015 0.02 0.03 0.04];   % piston velocities mm/s
Vp_vec  = Vp_mm_s * 1e-3;                              % to m/s

output_folder = fullfile('output_v3');
if ~exist(output_folder, 'dir'), mkdir(output_folder); end

%% -------------------- RUN ALL COMBINATIONS --------------------
all_results = cell(length(samples), length(geoms));
for ii = 1:length(samples)
    for jj = 1:length(geoms)
        fprintf('\n=== %s on %s ===\n', samples(ii).name, geoms(jj).label);
        res = bioprinting_algorithm_v3(samples(ii), geoms(jj), Vp_vec, ...
            'OutputFolder', output_folder, 'SaveCSV', true, 'PlotResults', true);
        all_results{ii, jj} = res;
    end
end

%% -------------------- MASTER SUMMARY CSV --------------------
% Single CSV that the SOP refers to: predicted slicer params at
% a single reference piston velocity (Vp = 0.01 mm/s -> moderate flow)
ref_Vp_idx = find(abs(Vp_mm_s - 0.01) < 1e-9, 1);
if isempty(ref_Vp_idx), ref_Vp_idx = 3; end

% Column order chosen so the operator sees Vp (input) immediately next to
% v_print (the resulting head speed) -- the only two numbers they really
% need to map a slicer command to the lab DIW behaviour.
rows = {};
for ii = 1:length(samples)
    for jj = 1:length(geoms)
        r = all_results{ii, jj};
        rows(end+1, :) = { ...
            r.sample.name, r.geom.label, ...
            r.Vp_mm_s(ref_Vp_idx), ...
            r.PL.v_print_mm_s(ref_Vp_idx), r.Cross.v_print_mm_s(ref_Vp_idx), ...
            r.Q_mm3_s(ref_Vp_idx), r.PL.v_mean_mm_s(ref_Vp_idx), ...
            r.PL.dP_kPa(ref_Vp_idx),    r.PL.w_line_mm,    r.PL.k_flow, ...
            r.Cross.dP_kPa(ref_Vp_idx), r.Cross.w_line_mm, r.Cross.k_flow, ...
            r.PL.h_layer_mm};  %#ok<AGROW>
    end
end
T_summary = cell2table(rows, 'VariableNames', { ...
    'sample','needle','Vp_mm_s', ...
    'PL_v_print_mm_s','Cross_v_print_mm_s', ...
    'Q_mm3_s','v_mean_mm_s', ...
    'PL_dP_kPa','PL_w_line_mm','PL_k_flow', ...
    'Cross_dP_kPa','Cross_w_line_mm','Cross_k_flow', ...
    'h_layer_mm'});
summary_csv = fullfile(output_folder, 'master_summary_v3.csv');
writetable(T_summary, summary_csv);
fprintf('\n[Saved] %s\n', summary_csv);

fprintf('\n--- DONE ---\nCSVs and plots in %s/\n', output_folder);
