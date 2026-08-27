function validate_v5()
%VALIDATE_V5  Regression check for bioprinting_algorithm_v5 against v4.
%
% Mirrors validate_v4.m in intent. Three things are asserted, and they are
% the three that would otherwise be discovered by a reviewer:
%
%   (1) FALLBACK IS INERT. v5 with no sample.N1_wall_Pa must reproduce v4
%       bit-for-bit on every slicer and flow quantity. If this fails, the
%       override plumbing added to v4 changed v4's own numbers, and every
%       previously reported run is suspect.
%   (2) CLOSURE MATHS. tanner_beta and ryan_johnson_recrit hit their known
%       analytical limits (N1 -> 0 gives exactly 0.13; n = 1 gives 2099.4).
%   (3) THE SWAP IS REPORTED. With N1 supplied, v5 must produce a different
%       beta from v4 AND populate beta_report with both, so the disclosure
%       is generated rather than recalled.
%
% Synthetic sample throughout - this runs anywhere, with no rheology data
% present. It does NOT substitute for the audited-dataset regression, which
% needs the C10-C25 force exports that are not on this machine.
%
% Author: T.M.C. Rodrigues (PEMM/COPPE/UFRJ) -- 2026-08-26

fprintf('\n=== validate_v5 ===\n');
tol = 1e-12;
nfail = 0;
global VV5_NCHK; VV5_NCHK = 0;

% Synthetic ink, roughly C15-SF5.5 shaped. Values are invented; only the
% agreement between solvers is under test.
sample = struct('name','SYNTH', ...
    'K_PL',267.102,'n_PL',0.24029, ...
    'eta0',2439.5,'etaInf',0.001,'lambda',4.09905,'m_Cross',0.987578, ...
    'rho',898,'Rrec_pct',36.20);
geom = struct('Rs',14.3e-3/2,'R_n',0.515e-3/2,'L_n',31.75e-3,'Ls',90e-3, ...
    'label','21G','h_factor',0.7);
Vp = [0.005 0.010 0.020] * 1e-3;      % m/s

quiet = {'SaveData',false,'SaveFigures',false,'SaveCSV',false,'PlotResults',false};

%% (1) FALLBACK IS INERT ------------------------------------------------
r4 = bioprinting_algorithm_v4(sample, geom, Vp, quiet{:});
ws = warning('off','bioprinting:v5:noNormalStress');
r5 = bioprinting_algorithm_v5(sample, geom, Vp, quiet{:});
warning(ws);

fields_vec = {'Q','dP_total','tau_w_needle','gamma_w_needle','Re_needle', ...
              'u_avg_needle','u_max_needle'};
for k = 1:numel(fields_vec)
    f = fields_vec{k};
    nfail = nfail + chk(sprintf('PL.%s', f), r4.PL.(f), r5.PL.(f), tol);
    nfail = nfail + chk(sprintf('Cross.%s', f), r4.Cross.(f), r5.Cross.(f), tol);
end
nfail = nfail + chk('slicer.PL.beta',    r4.slicer.PL.beta,    r5.slicer.PL.beta,    tol);
nfail = nfail + chk('slicer.PL.k_flow',  r4.slicer.PL.k_flow,  r5.slicer.PL.k_flow,  tol);
nfail = nfail + chk('slicer.PL.v_print', r4.slicer.PL.v_print_mm_s, r5.slicer.PL.v_print_mm_s, tol);
nfail = nfail + chk('slicer.CR.k_flow',  r4.slicer.CR.k_flow,  r5.slicer.CR.k_flow,  tol);

if ~strcmp(r5.closure,'v5-fallback-v4-heuristic')
    fprintf('  [FAIL] fallback tag = %s\n', r5.closure); nfail = nfail + 1;
else
    fprintf('  [ok]   fallback tagged: %s\n', r5.closure);
end

%% (2) CLOSURE MATHS ----------------------------------------------------
% N1 -> 0 must give exactly the 0.13 inelastic floor (Nickell 1974).
b0 = local_tanner(0, 500);
nfail = nfail + chk('tanner_beta(N1=0)', 0.13, b0, 1e-12);

% n = 1 must recover the Newtonian pipe value (2099.4, not 2100 exactly).
rc1 = local_ryan(1.0);
nfail = nfail + chk('ryan_johnson(n=1)', 6464*3^1.5/16, rc1, 1e-9);

% SHAPE check. An earlier internal note claimed the heuristic was merely
% "wrong-signed" against Ryan-Johnson. It is not: Ryan-Johnson is
% NON-MONOTONIC, peaking near n = 0.42, and stays inside a narrow band. The
% heuristic is what has a spurious monotonic collapse. Assert the real
% shape so nobody "corrects" the correlation back into a monotone rise.
nsweep = linspace(0.05, 1.0, 20000);
rsweep = local_ryan(nsweep);
[rmax, imax] = max(rsweep);
if abs(nsweep(imax) - 0.4165) > 0.01
    fprintf('  [FAIL] Ryan-Johnson maximum at n=%.4f, expected ~0.4165\n', nsweep(imax));
    nfail = nfail + 1;
else
    fprintf('  [ok]   Ryan-Johnson is non-monotonic: max %.1f at n=%.4f\n', rmax, nsweep(imax));
end
if (max(rsweep)/min(rsweep)) > 3
    fprintf('  [FAIL] Ryan-Johnson band wider than expected (%.2fx)\n', max(rsweep)/min(rsweep));
    nfail = nfail + 1;
else
    fprintf('  [ok]   Ryan-Johnson band %.0f-%.0f (%.2fx) vs heuristic %.0f-%.0f (%.2fx)\n', ...
        min(rsweep), max(rsweep), max(rsweep)/min(rsweep), ...
        2100*0.05^0.75, 2100, 1/0.05^0.75);
end
% The two working inks, against the numbers quoted in the source comments.
for nn = [0.24029 0.0879]
    ratio = local_ryan(nn) / (2100*nn^0.75);
    fprintf('         n=%.4f: heuristic %.0f vs Ryan-Johnson %.0f -> %.2fx under\n', ...
        nn, 2100*nn^0.75, local_ryan(nn), ratio);
end

%% (3) THE SWAP IS REPORTED ---------------------------------------------
s2 = sample; s2.N1_wall_Pa = 2527;    % Pa, order of a real CMC gel at the wall
r5t = bioprinting_algorithm_v5(s2, geom, Vp, quiet{:});
if ~strcmp(r5t.closure,'v5-tanner')
    fprintf('  [FAIL] tanner tag = %s\n', r5t.closure); nfail = nfail + 1;
else
    fprintf('  [ok]   tanner tagged: %s\n', r5t.closure);
end
br = r5t.beta_report;
if abs(br.beta_used_PL - br.beta_heuristic_PL) < 1e-6
    fprintf('  [FAIL] tanner beta did not differ from heuristic\n'); nfail = nfail + 1;
else
    fprintf('  [ok]   beta  heuristic %.4f -> tanner %.4f\n', ...
        br.beta_heuristic_PL, br.beta_used_PL);
    fprintf('         k_flow        %.4f -> %.4f\n', ...
        br.k_flow_heuristic_PL, br.k_flow_used_PL);
    fprintf('         extr. mult.   %.4f -> %.4f  (%+.1f%%)\n', ...
        br.extrusion_multiplier_heuristic_PL, br.extrusion_multiplier_used_PL, ...
        100*(br.extrusion_multiplier_used_PL/br.extrusion_multiplier_heuristic_PL - 1));
end

%% ----------------------------------------------------------------------
if nfail == 0
    fprintf('\n=== validate_v5: ALL CHECKS PASSED (%d numeric identities + 5 behavioural) ===\n\n', VV5_NCHK);
else
    fprintf('\n=== validate_v5: %d CHECK(S) FAILED ===\n\n', nfail);
end
end

function nf = chk(label, a, b, tol)
    global VV5_NCHK; VV5_NCHK = VV5_NCHK + 1;
    d = max(abs(a(:) - b(:)));
    scale = max(1, max(abs(a(:))));
    if d/scale > tol
        fprintf('  [FAIL] %-24s max abs diff %.3e (rel %.3e)\n', label, d, d/scale);
        nf = 1;
    else
        nf = 0;
    end
end

% Mirrors of the v5 local functions, so this file can test them without
% v5 having to export them. If these two ever diverge from v5's copies the
% test is worthless - keep them literally identical.
function beta = local_tanner(N1_w, tau_w)
    beta = 0.13 + (1 + 0.5*(N1_w/(2*tau_w))^2)^(1/6) - 1;
end
function Re_c = local_ryan(n)
    Re_c = 6464 * n .* (2+n).^((2+n)./(1+n)) ./ (3*n + 1).^2;
end
