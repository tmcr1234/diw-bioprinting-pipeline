function inks = inks_local()
%INKS_LOCAL  Project-local fitted ink parameters for run_solver_v4.
%   Copy this file to your PROJECT ROOT (the folder containing Export/) and
%   rename it inks_local.m, then fill in YOUR ink parameters. run_solver_v4
%   loads it from the current folder (pwd) so the Export/ tree stays a single
%   shared clone across projects.
%
%   Paste the K_PL/n_PL (Power-Law) and eta0/etaInf/lambda/m_Cross (Cross)
%   numbers straight from your latest fit output
%   (Analises/Python/Results/FitAll-*.txt). Update them the SAME DAY you re-fit.
%
%   Ramp1 = no pre-shear; Ramp2 = post pre-shear (200 1/s x 300 s); set a ramp
%   to [] if that ink lacks it. Rrec_pct = structural recovery at the
%   deposition wall shear rate (~200 s^-1), from your Relatorio_Recovery.txt.

inks(1).name     = 'MyInk';
inks(1).rho      = 1012;          % kg/m^3
inks(1).Rrec_pct = 40.0;
inks(1).Ramp1    = struct('K_PL',0,'n_PL',1,'eta0',0,'etaInf',0,'lambda',0,'m_Cross',1);
inks(1).Ramp2    = [];            % e.g. struct('K_PL',...) or [] if absent

% inks(2).name = '...'; inks(2).rho = ...; ...   % add as many inks as needed
end
