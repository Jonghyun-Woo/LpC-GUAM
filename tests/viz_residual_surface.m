% Observation view of the model-mismatch disturbance fit, for ALL UH at WH3.
% Each channel is shown over the plane of the TWO states that drive it most
% (ranked from the fitted response-surface coefficients):
%   colored surface = fitted e_max_c(z) = phi(dev/half)*beta   (the fit)
%   black scatter   = actual max_u|e_c(dev,u)| at sampled deviations on the plane
% lon channels: dX=e_u, dZ=e_w, dM=e_q ;  lat: dY=e_v, dL=e_p, dN=e_r.
% Figures saved per UH/axis under reachability_data/mc_verify/resid_surfaces/.

clear; clc;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

% ------------------------------- settings --------------------------------
stackDir = fullfile(root, 'reachability_data', 'guam_timestack', 'BRT');
outDir   = fullfile(root, 'reachability_data', 'mc_verify');
figDir   = fullfile(outDir, 'resid_surfaces');
WH_IDX   = 2;
NG       = 25;      % fit-surface grid nodes per axis
N_S      = 250;     % scatter (actual residual) samples per channel
N_NU     = 12;      % input samples per deviation
DT       = 0.01;
FT2M     = 0.3048;  % ft -> m (SI) for display
% -------------------------------------------------------------------------

if ~exist(figDir, 'dir'), mkdir(figDir); end
gridFile  = @(axis_name) fullfile(stackDir, sprintf('BRT_%s_MAT', upper(axis_name)), ...
                sprintf('GUAM_%s_BRT_UH1_WH%d_stack.mat', upper(axis_name), WH_IDX));
brt       = struct('lon', brt_axis_config('lon', gridFile('lon')), ...
                   'lat', brt_axis_config('lat', gridFile('lat')));
trimTable = load(fullfile(root, 'controller', 'trim_table_Poly_ConcatVer4p0.mat'));
GUAM      = LpC_GUAM(Config('althold', struct('dt', DT)));
QF        = load(fullfile(outDir, 'guam_disturbance_quadfit_WH2.mat'));
UH_ALL    = QF.UH_LIST(:)';

axinfo = struct( ...
    'name', {'lon', 'lat'}, ...
    'slab', {{'\Deltau [m/s]','\Deltaw [m/s]','\Deltaq [rad/s]','\Delta\theta [rad]'}, ...
             {'\Deltav [m/s]','\Deltap [rad/s]','\Deltar [rad/s]','\Delta\phi [rad]'}}, ...
    'clab', {{'dX','dZ','dM'}, {'dY','dL','dN'}}, ...
    'cunit',{{'m/s^2','m/s^2','rad/s^2'}, {'m/s^2','rad/s^2','rad/s^2'}}, ...
    'sscale',{[FT2M FT2M 1 1], [FT2M 1 1 1]}, ...
    'cscale',{[FT2M FT2M 1], [FT2M 1 1]});

for iu = 1:numel(UH_ALL)
    UH = UH_ALL(iu);
    for ai = 1:2
        ax = axinfo(ai).name;  gi = brt.(ax);
        XU0 = trimTable.XU0_interp(:, UH, WH_IDX);
        sched = struct();
        sched.X0 = XU0(1:12);  sched.U0 = XU0(13:25);
        sched.Ap = trimTable.(sprintf('Ap_%s_interp', ax))(:, :, UH, WH_IDX);
        sched.Bp = trimTable.(sprintf('Bp_%s_interp', ax))(:, :, UH, WH_IDX);
        sched.axis = ax;
        [sched.lb, sched.ub] = LivenessFilter.input_bounds(sched.U0, gi);
        prow = gi.state_rows;  half = (gi.grid_max - gi.grid_min) / 2;
        half_fit = half;   % OLS fitted with z=dev/half (full grid half)
        state0 = [0; 0; -100; sched.X0(1:3); sched.X0(10:12); sched.X0(4:6)];
        trimRef = state0(prow);
        beta = QF.(sprintf('beta_%s', ax))(:, :, iu);
        us = sched.lb' + (sched.ub - sched.lb)' .* rand(N_NU, numel(sched.lb));

        fig = figure('Color', 'w', 'Position', [50 50 1550 480], 'Visible', 'off');
        tl = tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
        for c = 1:3
            imp = state_importance(beta(:, c));
            [~, ord] = sort(imp, 'descend');
            d1 = ord(1);  d2 = ord(2);
            sc1 = axinfo(ai).sscale(d1);  sc2 = axinfo(ai).sscale(d2);  scz = axinfo(ai).cscale(c);

            % fitted surface over the dominant plane
            g1 = linspace(gi.grid_min(d1), gi.grid_max(d1), NG);
            g2 = linspace(gi.grid_min(d2), gi.grid_max(d2), NG);
            [X1, X2] = ndgrid(g1, g2);
            Zf = zeros(NG*NG, 4);  Zf(:, d1) = X1(:) / half_fit(d1);  Zf(:, d2) = X2(:) / half_fit(d2);
            Bfit = max(0, reshape(design_quad(Zf) * beta(:, c), NG, NG));

            % actual residual scatter on the same plane
            s1 = gi.grid_min(d1) + (gi.grid_max(d1) - gi.grid_min(d1)) * rand(N_S, 1);
            s2 = gi.grid_min(d2) + (gi.grid_max(d2) - gi.grid_min(d2)) * rand(N_S, 1);
            yS = zeros(N_S, 1);
            parfor s = 1:N_S
                dev = zeros(4, 1);  dev(d1) = s1(s);  dev(d2) = s2(s);
                st = state0;  st(prow) = trimRef + dev;
                em = 0;
                for j = 1:N_NU
                    [eng, srf] = actuator_cmd(sched, us(j, :)');
                    fnl = GUAM.state_derivative(st, eng, srf);
                    e   = fnl(prow) - (sched.Ap * dev + sched.Bp * us(j, :)');
                    em  = max(em, abs(e(c)));
                end
                yS(s) = em;
            end

            nexttile;  hold on;  grid on;  box on;
            surf(X1 * sc1, X2 * sc2, Bfit * scz, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
            % colour scatter by whether they fall inside the fitting region
            in_fit = abs(s1) <= half_fit(d1) & abs(s2) <= half_fit(d2);
            if any(in_fit)
                scatter3(s1(in_fit) * sc1,  s2(in_fit) * sc2,  yS(in_fit) * scz,  9, [0 0 0],        'filled', 'MarkerFaceAlpha', 0.55);
            end
            if any(~in_fit)
                scatter3(s1(~in_fit) * sc1, s2(~in_fit) * sc2, yS(~in_fit) * scz, 9, [0.75 0.1 0.1], 'filled', 'MarkerFaceAlpha', 0.35);
            end
            % fitting region boundary (dashed rectangle at z=0)
            bx = sc1 * half_fit(d1) * [-1  1  1 -1 -1];
            by = sc2 * half_fit(d2) * [-1 -1  1  1 -1];
            plot3(bx, by, zeros(1,5), 'k--', 'LineWidth', 1.2);
            xlabel(axinfo(ai).slab{d1});  ylabel(axinfo(ai).slab{d2});
            zlabel(sprintf('%s [%s]', axinfo(ai).clab{c}, axinfo(ai).cunit{c}));
            view(38, 26);  colormap(parula);
        end
        exportgraphics(fig, fullfile(figDir, sprintf('resid_%s_dom_UH%02d.png', ax, UH)), ...
                       'Resolution', 150);
        close(fig);
    end
    fprintf('UH%02d (%.0f ft/s) done\n', UH, trimTable.UH(UH));
end
fprintf('saved %d figures -> %s\n', 2*numel(UH_ALL), figDir);


% ===================== helpers ============================================
function imp = state_importance(b)
    lin = [2 3 4 5];  sq = [6 7 8 9];
    cross = {[10 11 12], [10 13 14], [11 13 15], [12 14 15]};
    imp = zeros(1, 4);
    for j = 1:4
        imp(j) = abs(b(lin(j))) + abs(b(sq(j))) + sum(abs(b(cross{j})));
    end
end

function Phi = design_quad(Z)
    z1 = Z(:,1);  z2 = Z(:,2);  z3 = Z(:,3);  z4 = Z(:,4);  o = ones(size(z1));
    Phi = [o, z1, z2, z3, z4, z1.^2, z2.^2, z3.^2, z4.^2, ...
           z1.*z2, z1.*z3, z1.*z4, z2.*z3, z2.*z4, z3.*z4];
end

function [eng, srf] = actuator_cmd(sched, u)
    U0 = sched.U0;
    if strcmp(sched.axis, 'lon')
        eng  = U0(5:13) + [u(1:8); u(9)];
        flap = U0(1) + u(11);  ail = U0(2);
        ele  = U0(3) + u(10);  rud = U0(4);
    else
        eng  = U0(5:13) + [u(1:8); 0];
        flap = U0(1);          ail = U0(2) + u(9);
        ele  = U0(3);          rud = U0(4) + u(10);
    end
    srf = [flap - ail; flap + ail; ele; ele; rud];
end
