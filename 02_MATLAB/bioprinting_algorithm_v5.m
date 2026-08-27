function results = bioprinting_algorithm_v5(sample, geom, Vp_vec, varargin)
%BIOPRINTING_ALGORITHM_V5  v4 physics with physically anchored closures.
%
% v5 changes exactly TWO things relative to v4, both of which v4 flagged in
% its own source as in-house heuristics with no published support:
%
%   1. DIE SWELL. v4 used beta = 0.30*(1-n), which has no source and is
%      wrong-signed (it makes swell grow as the ink becomes more shear
%      thinning, whereas a lower n means flow closer to plug flow, less
%      stored recoverable strain, and therefore LESS swell). v5 uses the
%      Tanner closure keyed on the MEASURED first normal-stress difference:
%
%          d/D  = 0.13 + [ 1 + 0.5*( N1_w / (2*tau_w) )^2 ]^(1/6)
%          beta = d/D - 1
%
%      (Tanner 1970, J Polym Sci A-2 8(12):2067-2078,
%       doi:10.1002/pol.1970.160081203; the 0.13 inelastic floor is
%       Nickell, Tanner & Yamada 1974, J Fluid Mech 65(1):189-206,
%       doi:10.1017/S0022112074001339.)
%
%      With N1_w = 0 the closure collapses to beta = 0.13 exactly, so the
%      inelastic floor is not a separate branch - it is the N1 -> 0 limit.
%
%      HONEST LIMITATION, DO NOT TUNE AROUND IT: Tanner describes a FREE JET
%      leaving a long capillary. A DIW road is deposited against a substrate
%      at a finite stand-off and is still swelling when it lands. On the one
%      dataset where this group has both N1 and a measured road width, the
%      Tanner estimate remains below the measured swell. That residual is a
%      definition mismatch between "free-jet diameter" and "deposited road
%      width", not a bad fit. Do not adjust the 0.13 or the exponent to
%      close it; report both numbers instead.
%
%   2. CRITICAL REYNOLDS. v4 used Re_c = 2100*n^0.75, also unsourced. v5
%      uses Ryan & Johnson (1959, AIChE J 5(4):433-435,
%      doi:10.1002/aic.690050407):
%
%          Re_c = 6464*n*(2+n)^((2+n)/(1+n)) / (3n+1)^2
%
%      which reduces to 2099.4 at n = 1. Note it is NON-MONOTONIC - it
%      peaks at 2397 near n = 0.42 and stays within 992-2397 over the whole
%      practical range - whereas the heuristic collapses monotonically to
%      222 by n = 0.05. So the heuristic under-predicts by 3.13x at
%      n = 0.24 and 4.31x at n = 0.088. This changes only the laminar/
%      transitional verdict string and the Re_crit line on the flow-regime
%      figure; the generalised Reynolds numbers here are ~1e-4, so the
%      verdict stays "laminar" under either correlation. The point of the
%      swap is that the reported number is now defensible, not that any
%      conclusion moves.
%
% EVERYTHING ELSE IS v4. v5 does not reimplement the physics - it resolves
% the two closures and then delegates to bioprinting_algorithm_v4 through
% its 'BetaPL' / 'BetaCross' / 'RecritFcn' override parameters. That is
% deliberate: forking 600 lines of identical solver code is how the two
% copies of this pipeline drift apart. validate_v5.m asserts that v5 with
% no N1 supplied reproduces v4 bit-for-bit.
%
% --------------------------------------------------------------------
% INPUTS
%   sample, geom, Vp_vec - as bioprinting_algorithm_v4, PLUS optionally:
%
%       sample.N1_wall_Pa - measured first normal-stress difference at the
%           needle wall shear rate (Pa). Scalar. Obtain it from a
%           force-augmented flow curve; see antpar_io.read_flow_curve_with_force
%           and normal_stress_from_force in the Python layer, which apply
%           N1 = 2F/(pi*R^2) and the min-baseline tare.
%
%   If sample.N1_wall_Pa is absent, empty or NaN, v5 falls back to the v4
%   heuristic AND SAYS SO in results.closure and in every report it writes.
%   The fallback exists so v5 is runnable on inks that have no normal-force
%   measurement; it is not an endorsement of the heuristic.
%
% NAME-VALUE (in addition to every v4 option, which is passed through)
%   'BetaRefVp_mm_s' (default: median of Vp_vec) - the operating point whose
%       wall shear stress anchors the Tanner ratio. beta is a single number
%       per run in this solver, so one operating point has to be chosen;
%       make it the deposition point you actually print at.
%
% OUTPUT
%   results - the v4 results struct, plus:
%       .closure        'v5-tanner' or 'v5-fallback-v4-heuristic'
%       .beta_report    struct comparing the heuristic and Tanner betas and
%                       their slicer consequences, so the disclosure a
%                       manuscript needs is generated rather than recalled.
%
% Author: T.M.C. Rodrigues (PEMM/COPPE/UFRJ) -- 2026-08-26
% --------------------------------------------------------------------

%% -------------------- SPLIT OFF v5-ONLY OPTIONS --------------------
[betaRefVp, v4args] = pop_option(varargin, 'BetaRefVp_mm_s', median(Vp_vec)*1e3);

n  = sample.n_PL;
mC = sample.m_Cross;

beta_heur_PL = max(0, 0.30*(1 - n));
beta_heur_CR = max(0, 0.30*(1 - mC));

%% -------------------- PASS 1: resolve tau_w --------------------------
% The Tanner ratio needs the wall shear stress, which only the solver knows.
% Pass 1 runs silently with the v4 closure purely to obtain tau_w; nothing
% it writes is kept.
has_N1 = isfield(sample,'N1_wall_Pa') && ~isempty(sample.N1_wall_Pa) ...
         && all(isfinite(sample.N1_wall_Pa));

if has_N1
    N1_report = sample.N1_wall_Pa;
else
    N1_report = NaN;   % MATLAB evaluates both branches of a call, so this
end                    % cannot be folded into the report line below.
if has_N1
    probe = bioprinting_algorithm_v4(sample, geom, Vp_vec, ...
        'SaveData', false, 'SaveFigures', false, 'SaveCSV', false, ...
        'PlotResults', false);
    [~, iref] = min(abs(probe.Vp_mm_s - betaRefVp));
    tau_w_PL = probe.PL.tau_w_needle(iref);
    tau_w_CR = probe.Cross.tau_w_needle(iref);

    beta_PL = tanner_beta(sample.N1_wall_Pa, tau_w_PL);
    beta_CR = tanner_beta(sample.N1_wall_Pa, tau_w_CR);
    tag = 'v5-tanner';
else
    iref = NaN; tau_w_PL = NaN; tau_w_CR = NaN;
    beta_PL = beta_heur_PL;
    beta_CR = beta_heur_CR;
    tag = 'v5-fallback-v4-heuristic';
end

%% -------------------- PASS 2: the real run ---------------------------
results = bioprinting_algorithm_v4(sample, geom, Vp_vec, v4args{:}, ...
    'BetaPL', beta_PL, 'BetaCross', beta_CR, ...
    'RecritFcn', @ryan_johnson_recrit, 'ClosureTag', tag);

%% -------------------- DISCLOSURE BLOCK -------------------------------
% Generated, not remembered. Paste straight into a Methods/limitations
% paragraph: it carries both closures and what the swap costs at the bench.
results.closure = tag;
results.beta_report = struct( ...
    'closure',            tag, ...
    'beta_heuristic_PL',  beta_heur_PL, ...
    'beta_heuristic_CR',  beta_heur_CR, ...
    'beta_used_PL',       beta_PL, ...
    'beta_used_CR',       beta_CR, ...
    'N1_wall_Pa',         N1_report, ...
    'tau_wall_PL_Pa',     tau_w_PL, ...
    'tau_wall_CR_Pa',     tau_w_CR, ...
    'beta_ref_Vp_mm_s',   betaRefVp, ...
    'k_flow_heuristic_PL', kflow(beta_heur_PL, sample.Rrec_pct), ...
    'k_flow_used_PL',      kflow(beta_PL,      sample.Rrec_pct), ...
    'extrusion_multiplier_heuristic_PL', 1/kflow(beta_heur_PL, sample.Rrec_pct), ...
    'extrusion_multiplier_used_PL',      1/kflow(beta_PL,      sample.Rrec_pct));

if strcmp(tag, 'v5-fallback-v4-heuristic')
    warning('bioprinting:v5:noNormalStress', ...
        ['[%s] No sample.N1_wall_Pa supplied - v5 fell back to the v4 ' ...
         'die-swell HEURISTIC beta = 0.30*(1-n). The slicer numbers in ' ...
         'this run inherit a closure with no published source and the ' ...
         'wrong sign. Measure N1 before quoting them.'], sample.name);
end
end

%% ==================================================================
%% LOCAL FUNCTIONS
%% ==================================================================
function beta = tanner_beta(N1_w, tau_w)
%TANNER_BETA  Die-swell ratio minus one, from the measured N1 (Tanner 1970).
%   d/D  = 0.13 + [ 1 + 0.5*( N1_w/(2*tau_w) )^2 ]^(1/6)
%   beta = d/D - 1
%   N1_w = 0 gives exactly 0.13 (the Nickell 1974 inelastic floor).
%   Free-jet definition - see the limitation note in the header.
    if ~isfinite(tau_w) || tau_w <= 0
        error('tanner_beta:badStress', 'tau_w must be finite and positive.');
    end
    if N1_w < 0
        error('tanner_beta:negativeN1', ...
            ['N1_w = %g Pa is negative. A negative first normal-stress ' ...
             'difference in steady shear on a CMC gel almost always means ' ...
             'an un-tared transducer, not a real measurement - apply the ' ...
             'min-baseline correction first (see baseline_correct).'], N1_w);
    end
    ratio = N1_w / (2*tau_w);
    beta  = 0.13 + (1 + 0.5*ratio^2)^(1/6) - 1;
end

function Re_c = ryan_johnson_recrit(n)
%RYAN_JOHNSON_RECRIT  Critical generalised Reynolds number for a power-law
%   fluid in a pipe. Ryan & Johnson 1959, AIChE J 5(4):433-435,
%   doi:10.1002/aic.690050407.
%       Re_c = 6464*n*(2+n)^((2+n)/(1+n)) / (3n+1)^2
%   Reduces to 2099.4 at n = 1. NON-MONOTONIC: peaks at 2397 near n = 0.42
%   and stays in 992-2397 over n in [0.05,1], where the v4 heuristic it
%   replaces collapses monotonically to 222.
    Re_c = 6464 * n .* (2+n).^((2+n)./(1+n)) ./ (3*n + 1).^2;
end

function k = kflow(beta, Rrec_pct)
    k = (1+beta)^2 * sqrt(Rrec_pct/100);
end

function [val, rest] = pop_option(args, name, default)
%POP_OPTION  Remove one name-value pair from a varargin cell, leaving the
%   remainder for pass-through. Case-insensitive on the name.
    val = default; rest = args; i = 1; keep = true(1, numel(args));
    while i < numel(args)
        if (ischar(args{i}) || isstring(args{i})) && strcmpi(char(args{i}), name)
            val = args{i+1};
            keep(i) = false; keep(i+1) = false;
            i = i + 2;
        else
            i = i + 1;
        end
    end
    rest = args(keep);
end
