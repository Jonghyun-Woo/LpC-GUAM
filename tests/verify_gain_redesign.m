% verify_gain_redesign - does regenerating the gains reproduce the stored table?
%
% Run this FIRST, with RSLQRConfig.update_gains = true and Qlon0/Rlon0/Wlon0
% left at their original values. If the regenerated gains match the ones
% shipped in trim_table_Poly_ConcatVer4p0.mat, the design path (set_gains ->
% ctrl_lon -> lqr) is wired correctly and any later change to Q can be
% trusted. If they do not match, tuning Q is meaningless because the design
% being solved is not the one that produced the table.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

if ~RSLQRConfig.update_gains
    fprintf(['\nRSLQRConfig.update_gains is FALSE - the stored gains are being\n' ...
             'loaded, so there is nothing to compare. Set it to true first.\n']);
    return;
end

fprintf('Qlon0 = [%s]\n', strjoin(compose('%g', RSLQRConfig.Qlon0(:)'), ' '));
fprintf('Rlon0 = [%s]\n', strjoin(compose('%g', RSLQRConfig.Rlon0(:)'), ' '));
fprintf('\nregenerating 28 x 3 = 84 design points ...\n');

cfg = Config('trim_schedule', struct('filter_mode','off','T_seg',2.0));
tic;  ctl = Controller(cfg.controller, cfg.sim.dt);  t = toc;
new = ctl.baseline_controller;
fprintf('done in %.0f s\n', t);

old = load('trim_table_Poly_ConcatVer4p0.mat');

% Two classes of matrix, two tolerances.
%   gains      : outputs of the Riccati solve and the hardcoded servo
%                structure. Deterministic given Q/R, so they must agree to
%                solver precision.
%   linearised : Ap, Bp and the allocation effectiveness B come from a
%                central-difference of the aero model. The perturbation step
%                was changed at some point (see the commented-out "Ben's
%                original hardcoded steps" in poly_aero_wrapper_Mod_du), so a
%                few 1e-4 of disagreement is finite-difference noise, not a
%                wiring fault. Only a gross mismatch matters here.
gains = {'Ki','Kx','Kv','F','G','C','Cv','W'};
lins  = {'B','Ap','Bp'};
TOL_G = 1e-3;      % relative, gains
TOL_L = 1e-2;      % relative, numerically differentiated matrices

fprintf('\n%-6s %-11s %14s %14s %10s\n', ...
        'matrix', 'class', 'max |new-old|', 'rel. error', 'verdict');
fprintf('%s\n', repmat('-', 1, 60));
worst_g = 0;  worst_l = 0;
for a = 1:numel(gains) + numel(lins)
    if a <= numel(gains)
        n = gains{a};   cls = 'gain';         tol = TOL_G;
    else
        n = lins{a-numel(gains)};  cls = 'linearised';  tol = TOL_L;
    end
    A = new.LON.(n);
    Bm = old.([n '_lon_interp']);
    if ~isequal(size(A), size(Bm))
        fprintf('%-6s %-11s %14s %14s %10s\n', n, cls, 'SIZE MISMATCH', ...
                mat2str(size(A)), mat2str(size(Bm)));
        continue;
    end
    d   = max(abs(A(:) - Bm(:)));
    rel = d / max(max(abs(Bm(:))), eps);
    if a <= numel(gains), worst_g = max(worst_g, rel);
    else,                 worst_l = max(worst_l, rel);  end
    if rel < 1e-12,   v = 'exact';
    elseif rel < tol, v = 'ok';
    else,             v = 'DIFFERS';  end
    fprintf('%-6s %-11s %14.3e %14.3e %10s\n', n, cls, d, rel, v);
end

fprintf('\nworst relative error   gains %.3e (tol %.0e) | linearised %.3e (tol %.0e)\n', ...
        worst_g, TOL_G, worst_l, TOL_L);
if worst_g < TOL_G && worst_l < TOL_L
    fprintf(['\nPASS - the design path reproduces the shipped gains.\n' ...
             'Qlon0 can be changed and the result trusted.\n']);
elseif worst_g >= TOL_G
    fprintf(['\nFAIL - the GAINS differ. The design being solved is not the\n' ...
             'one that produced the table; do not tune Q until this is fixed.\n']);
else
    fprintf(['\nCHECK - gains agree but the linearised matrices are far off.\n' ...
             'Likely a different aero model or perturbation step, not a\n' ...
             'wiring fault, but worth confirming before relying on Ap/Bp.\n']);
end

%% Where the gains sit, for reference when tuning
fprintf('\nKi diagonal vs sqrt(Q/R) at a few trims (should agree closely)\n');
fprintf('%5s %8s | %9s %9s %9s | %9s %9s %9s\n', ...
        'trim','u [f/s]','Ki(1,1)','Ki(2,2)','Ki(3,3)', ...
        'sqrt Q1','sqrt Q2','sqrt Q3');
q = RSLQRConfig.Qlon0;  r = RSLQRConfig.Rlon0;
for k = [1 8 16 28]
    K = new.LON.Ki(:,:,k,2);            % WH index 2 (w = 0), the row actually used
    fprintf('%5d %8.1f | %9.4f %9.4f %9.4f | %9.4f %9.4f %9.4f\n', ...
            k, new.UH(k), K(1,1), K(2,2), K(3,3), ...
            sqrt(q(1)/r(1)), sqrt(q(2)/r(2)), sqrt(q(3)/r(3)));
end
