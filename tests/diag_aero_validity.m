% Aero-validity sweep: perturb each effector from trim over a wide range and
% evaluate the polynomial aero model. The effector's PRIMARY response (the
% force/moment component it drives most near trim) should stay monotone over the
% valid envelope; where its slope reverses (control reversal) or blows up, the
% aero database has left its valid range -- the practical cap on the control
% perturbation Delta. Marks the current yml Delta and the physical limit.
% Figures per UH under reachable_data/mc_verify/aero_validity/.

clear; clc;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

% ------------------------------- settings --------------------------------
outDir   = fullfile(root, 'reachable_data', 'mc_verify');
figDir   = fullfile(outDir, 'aero_validity');
UH_SWEEP = 1:20;        % all trim airspeeds 1-20
DETAIL   = [1 10 20];   % which UH get the full 6-panel sweep figure
WH_IDX   = 3;
NP       = 241;
ALT      = 100;     % ft
% -------------------------------------------------------------------------

if ~exist(figDir, 'dir'), mkdir(figDir); end
trimTable = load(fullfile(root, 'controller', 'trim_table_Poly_ConcatVer4p0.mat'));
GUAM      = LpC_GUAM(Config('althold', struct('dt', 0.01)));
[rho, a]  = GUAM.environment.atmosphere(ALT);
units = GUAM.units;  veh = GUAM.vehicleConfig;

rpm = 60 / (2*pi);  r2d = 180/pi;
comp = {'Fx','Fy','Fz','Mx','My','Mz'};
% effector: perturbation directions on engine(9)/surface(5), display scale, current Delta, sweep half-range, phys limit
eff = struct( ...
  'name', {'lift','pusher','flap','aileron','elevator','rudder'}, ...
  'deng', {[ones(8,1);0], [zeros(8,1);1], zeros(9,1), zeros(9,1), zeros(9,1), zeros(9,1)}, ...
  'dsrf', {zeros(5,1), zeros(5,1), [1;1;0;0;0], [-1;1;0;0;0], [0;0;1;1;0], [0;0;0;0;1]}, ...
  'unit', {'RPM','RPM','deg','deg','deg','deg'}, ...
  'prim', {3, 1, 3, 4, 5, 6}, ...                                  % design primary: Fz,Fx,Fz,Mx,My,Mz
  'disp', {rpm, rpm, r2d, r2d, r2d, r2d}, ...
  'dcur', {12.5664, 31.4159, 0.1745, 0.1745, 0.1745, 0.1745}, ...   % current yml Delta (rad/s | rad)
  'half', {600/rpm, 800/rpm, 40/r2d, 40/r2d, 40/r2d, 40/r2d}, ...   % sweep half-range
  'plim', {1600/rpm, 2000/rpm, 30/r2d, 30/r2d, 30/r2d, 30/r2d});    % physical limit (abs)

UHv = trimTable.UH(UH_SWEEP);
nE = numel(eff);  nU = numel(UH_SWEEP);
Lp = nan(nE, nU);  Ln = nan(nE, nU);         % linear-valid Delta (+ / -) per effector per UH

fprintf('%-4s %-9s %-6s %10s %10s %10s   (units: RPM|deg)\n', 'UH','effector','prim','curDelta','aVal+','aVal-');
for iu = 1:nU
    uh = UH_SWEEP(iu);
    XU0 = trimTable.XU0_interp(:, uh, WH_IDX);
    X0 = XU0(1:12);  U0 = XU0(13:25);
    x_aero = [X0(1:3); X0(10:12)];
    eng0 = U0(5:13);
    srf0 = [U0(1)-U0(2); U0(1)+U0(2); U0(3); U0(3); U0(4)];
    makefig = ismember(uh, DETAIL);
    if makefig
        fig = figure('Color','w','Position',[40 40 1500 780],'Visible','off');
        tl = tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
    end
    for ie = 1:nE
        del = linspace(-eff(ie).half, eff(ie).half, NP)';
        FM = zeros(NP, 6);  vmask = false(NP, 1);
        for k = 1:NP
            eng = eng0 + del(k)*eff(ie).deng;
            srf = srf0 + del(k)*eff(ie).dsrf;
            eng = max(eng, 0);                             % rotors cannot spin backwards
            [Fb, Mb, Val] = run_LpC_aero(x_aero, eng, srf, rho, a, units, veh);
            FM(k,:) = [Fb(:)', Mb(:)'];
            vmask(k) = ~any(logical(Val(:)));              % valid = both Validity flags false
        end
        pc = eff(ie).prim;  P = FM(:, pc);
        tp = firstinvalid(del, vmask, +1);                 % aero-DB validity edge (authoritative)
        tn = firstinvalid(del, vmask, -1);
        Lp(ie, iu) = tp * eff(ie).disp;  Ln(ie, iu) = tn * eff(ie).disp;

        if makefig
            nexttile;  hold on; grid on; box on;
            plot(del*eff(ie).disp, P, 'b-', 'LineWidth', 1.3);  yl = ylim;
            plot( eff(ie).dcur*eff(ie).disp*[1 1], yl, 'k--');
            plot(-eff(ie).dcur*eff(ie).disp*[1 1], yl, 'k--');
            if ~isnan(tp), plot(tp*eff(ie).disp*[1 1], yl, 'm-', 'LineWidth', 1.3); end
            if ~isnan(tn), plot(tn*eff(ie).disp*[1 1], yl, 'm-', 'LineWidth', 1.3); end
            xlabel(sprintf('\\Delta %s [%s]', eff(ie).name, eff(ie).unit));  ylabel(comp{pc});
            title(sprintf('%s -> %s   (curr \\Delta = %.0f %s)', eff(ie).name, comp{pc}, ...
                          eff(ie).dcur*eff(ie).disp, eff(ie).unit));
            fprintf('%-4d %-9s %-6s %10.0f %10s %10s\n', uh, eff(ie).name, comp{pc}, ...
                    eff(ie).dcur*eff(ie).disp, fmt(tp, eff(ie).disp), fmt(tn, eff(ie).disp));
        end
    end
    if makefig
        title(tl, sprintf('Aero-validity sweep  (UH%d = %.0f ft/s, WH3):  blue = design-primary,  dashed = current \\Delta,  magenta = aero-DB validity edge (frozen forces past here)', ...
                          uh, trimTable.UH(uh)));
        exportgraphics(fig, fullfile(figDir, sprintf('aero_validity_UH%02d.png', uh)), 'Resolution', 150);
        close(fig);  fprintf('\n');
    end
end

% ---- summary: linear-valid Delta vs airspeed for all UH 1-20 --------------
figS = figure('Color','w','Position',[40 40 1500 780],'Visible','off');
tlS = tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
for ie = 1:nE
    nexttile;  hold on; grid on; box on;
    plot(UHv, Lp(ie,:),  'b.-', 'LineWidth', 1.1);
    plot(UHv, -Ln(ie,:), 'r.-', 'LineWidth', 1.1);
    yline(eff(ie).dcur*eff(ie).disp, 'k--', 'current');
    xlabel('UH [ft/s]');  ylabel(sprintf('valid \\Delta %s [%s]', eff(ie).name, eff(ie).unit));
    title(eff(ie).name);
    if ie == 1, legend({'lin+','|lin-|','current'}, 'Location','best'); end
end
title(tlS, 'Aero-DB validity-valid \Delta vs airspeed (UH1-20, WH3);  gaps = no speed/prop trip (surfaces, or beyond sweep)');
exportgraphics(figS, fullfile(figDir, 'aero_validity_summary.png'), 'Resolution', 150);
close(figS);
fprintf('saved sweeps + summary -> %s\n', figDir);


% ===================== helpers ============================================
function d = firstinvalid(del, vmask, dir)
    % first |del| in direction dir where the aero Validity flags trip (invalid)
    d = NaN;
    if dir > 0, idx = find(del > 0); else, idx = flip(find(del < 0)); end
    for k = idx(:)'
        if ~vmask(k), d = del(k); return; end
    end
end

function s = fmt(d, disp)
    if isnan(d), s = '   --'; else, s = sprintf('%.0f', d*disp); end
end
