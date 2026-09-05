% Quantify the linear-vs-nonlinear model residual as a DISTURBANCE ENVELOPE.
%
%   e(dev,u) = f_nl(x0+dev, u_trim+u) - (Ap*dev + Bp*u)   [state-rate units]
%   Y(dev)   = max_u |e(dev,u)|                            (bang-bang input set)
%
% For each channel/UH: evaluate Y over a uniform grid of state deviations
% (full grid), fit a full 2nd-order polynomial in normalised z = dev./half,
% then INFLATE each channel by s = max_k Y/fit so the fit upper-bounds Y:
%   e_c(z) = phi(z) * (s_c * beta_c)
%   phi(z) = [1, z1..z4, z1^2..z4^2, z1z2,z1z3,z1z4,z2z3,z2z4,z3z4]  (15 terms)

clear; clc;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

% ------------------------------- settings --------------------------------
stackDir    = fullfile(root, 'reachable_data');
outDir      = fullfile(root, 'reachable_data', 'mc_verify');
UH_LIST     = 1:20;
WH_IDX      = 3;
NG_PER_DIM  = 7;       % uniform grid per dim: 7^4 = 2401 state points
N_NU        = 24;      % random input samples per deviation (bang-bang envelope)
DT          = 0.01;
REGION_FRAC = 1.0;     % fit over the FULL grid so r-growth is captured (prev trim-only: 0.5)
FIT_FLOOR   = 1e-3;    % ignore near-zero fit when computing inflation scalar
% -------------------------------------------------------------------------

if ~exist(outDir, 'dir'), mkdir(outDir); end
brt       = brt_setup(read_yml(fullfile(stackDir, 'guam_analysis_config.yml')));
trimTable = load(fullfile(root, 'controller', 'trim_table_Poly_ConcatVer4p0.mat'));
GUAM      = LpC_GUAM(Config('althold', struct('dt', DT)));

axes_ = struct('name', {'lon', 'lat'}, 'chan', {{'u','w','q','th'}, {'v','p','r','ph'}});
D = struct('lon', struct([]), 'lat', struct([]));

fprintf('Envelope Y = max_u |f_nl(x0+dev,u_trim+u) - (Ap*dev + Bp*u)|\n');
fprintf('OLS quad fit + per-channel inflation to upper-bound Y.\n');
fprintf('cover%% = %% grid where envelope >= Y (target ~100);  r-chan highlighted:\n\n');

for ai = 1:2
    axisName = axes_(ai).name;
    gridInfo = brt.(axisName);
    half     = (gridInfo.gmax - gridInfo.gmin) / 2;
    cen      = (gridInfo.gmax + gridInfo.gmin) / 2;
    prow     = gridInfo.prow;

    fprintf('== %s ==  channels [%s]\n', upper(axisName), strjoin(axes_(ai).chan, ' '));
    for ui = 1:numel(UH_LIST)
        uh  = UH_LIST(ui);
        XU0 = trimTable.XU0_interp(:, uh, WH_IDX);
        U0  = XU0(13:25);
        Ap  = trimTable.(sprintf('Ap_%s_interp', axisName))(:, :, uh, WH_IDX);
        Bp  = trimTable.(sprintf('Bp_%s_interp', axisName))(:, :, uh, WH_IDX);

        state0  = [0; 0; -100; XU0(1:3); XU0(10:12); XU0(4:6)];
        trimRef = state0(prow);

        % perturbation bounds about trim (bang-bang box), then random samples
        sched      = struct('axis', axisName, 'U0', U0);
        uTrim      = U0(gridInfo.U0_idx);
        sched.lb   = max(gridInfo.pl, uTrim + gridInfo.Dlo) - uTrim;
        sched.ub   = min(gridInfo.pu, uTrim + gridInfo.Dhi) - uTrim;
        us         = sched.lb' + (sched.ub - sched.lb)' .* rand(N_NU, numel(sched.lb));

        % uniform grid over inner REGION_FRAC of each dim
        edges = arrayfun(@(lo,hi) linspace(lo,hi,NG_PER_DIM)', ...
                         cen - REGION_FRAC*half, cen + REGION_FRAC*half, ...
                         'UniformOutput', false);
        [G1,G2,G3,G4] = ndgrid(edges{1}, edges{2}, edges{3}, edges{4});
        devs    = [G1(:), G2(:), G3(:), G4(:)];
        N_DELTA = size(devs, 1);

        % ===== 이전(trim-only) 방식 — 보존: 입력 섭동 없이 트림에서만 잔차 =====
        % eng = U0(5:13);
        % srf = [U0(1)-U0(2); U0(1)+U0(2); U0(3); U0(3); U0(4)];
        % Y = zeros(N_DELTA, 4);
        % parfor k = 1:N_DELTA
        %     st  = state0;  st(prow) = trimRef + devs(k,:)';
        %     fnl = GUAM.state_derivative(st, eng, srf);
        %     e   = fnl(prow) - Ap * devs(k,:)';   % Bp*u 없음, u=trim 고정
        %     Y(k,:) = abs(e)';
        % end
        % ===== 현재(envelope) 방식 — Y = max_u |f_nl - (Ap*dev + Bp*u)| (bang-bang) =====
        Y = zeros(N_DELTA, 4);
        parfor k = 1:N_DELTA
            st = state0;  st(prow) = trimRef + devs(k,:)';
            em = zeros(1, 4);
            for j = 1:N_NU
                [eng, srf] = actuator_cmd(sched, us(j,:)');
                fnl = GUAM.state_derivative(st, eng, srf);
                e   = fnl(prow) - (Ap * devs(k,:)' + Bp * us(j,:)');
                em  = max(em, abs(e)');
            end
            Y(k,:) = em;
        end

        zs   = (devs - cen') ./ half';
        Phi  = design_quad(zs);
        % ===== 이전 방식 — OLS만 (inflation 없음, 상계 아님) =====
        % beta = zeros(15, 4);  R2 = zeros(1, 4);
        % for c = 1:4
        %     beta(:, c) = Phi \ Y(:, c);
        %     res = Y(:, c) - Phi * beta(:, c);
        %     R2(c) = 1 - sum(res.^2) / max(sum((Y(:,c) - mean(Y(:,c))).^2), eps);
        % end
        % ===== 현재 — OLS + per-channel inflation scalar 로 상계화 =====
        beta = zeros(15, 4);  R2 = zeros(1, 4);  infl = ones(1, 4);  cover = zeros(1, 4);
        for c = 1:4
            b = Phi \ Y(:, c);
            res  = Y(:, c) - Phi * b;
            R2(c) = 1 - sum(res.^2) / max(sum((Y(:,c) - mean(Y(:,c))).^2), eps);
            F   = max(Phi * b, 0);                       % OLS shape, clipped >=0
            sel = F > FIT_FLOOR * max(F);
            infl(c) = max([1, max(Y(sel, c) ./ F(sel))]);% scalar so s*F >= Y
            beta(:, c) = infl(c) * b;
            cover(c) = 100 * mean(max(Phi * beta(:,c), 0) >= Y(:, c));
        end
        rec = struct('uh', uh, 'beta', beta, 'R2', R2, 'infl', infl, 'cover', cover);
        if isempty(D.(axisName)), D.(axisName) = rec; else, D.(axisName)(ui) = rec; end

        rc = 3;  % r / q channel index of interest per axis
        fprintf('  UH%-2d  %s: R2=%.2f infl=%.2f cover=%3.0f%% (Ymax=%.3f)\n', ...
                uh, axes_(ai).chan{rc}, R2(rc), infl(rc), cover(rc), max(Y(:,rc)));
    end
    fprintf('\n');
end

% ---- save ----------------------------------------------------------------
nM = numel(UH_LIST);
beta_lon = zeros(15, 4, nM);  beta_lat = zeros(15, 4, nM);
infl_lon = zeros(4, nM);      infl_lat = zeros(4, nM);
for ui = 1:nM
    beta_lon(:, :, ui) = D.lon(ui).beta;  infl_lon(:, ui) = D.lon(ui).infl';
    beta_lat(:, :, ui) = D.lat(ui).beta;  infl_lat(:, ui) = D.lat(ui).infl';
end
half_lon = ((brt.lon.gmax - brt.lon.gmin) / 2)';
half_lat = ((brt.lat.gmax - brt.lat.gmin) / 2)';
feature  = 'phi(z)=[1, z1..z4, z1^2..z4^2, z1z2,z1z3,z1z4,z2z3,z2z4,z3z4], z=dev/half; beta inflated to upper-bound max_u|e|';
save(fullfile(outDir, 'guam_disturbance_quadfit.mat'), ...
     'beta_lon', 'beta_lat', 'half_lon', 'half_lat', 'infl_lon', 'infl_lat', 'UH_LIST', 'feature');
fprintf('saved -> %s\n', fullfile(outDir, 'guam_disturbance_quadfit.mat'));
fprintf('  channel order lon [u w q theta], lat [v p r phi];  beta(:,c,uh) = 15x1 (inflated)\n');


% ===================== helpers ============================================
function Phi = design_quad(Z)
    z1=Z(:,1); z2=Z(:,2); z3=Z(:,3); z4=Z(:,4); o=ones(size(z1));
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

function brt = brt_setup(yml)
    dyn = yml.dynamics;
    lon = yml.longitudinal;
    brt.lon = axis_setup(lon, {'u', 'w', 'q', 'theta'});
    brt.lon.prow   = [4; 6; 11; 8];
    brt.lon.U0_idx = [5:12, 13, 3, 1];
    brt.lon.Dlo = [repmat(lon.input_min_Pi, 8, 1); lon.input_min_Pi_p; lon.input_min_delta_e; lon.input_min_delta_f];
    brt.lon.Dhi = [repmat(lon.input_max_Pi, 8, 1); lon.input_max_Pi_p; lon.input_max_delta_e; lon.input_max_delta_f];
    brt.lon.pl  = [repmat(dyn.input_min_Pi, 8, 1); dyn.input_min_Pi_p; dyn.input_min_delta_e; dyn.input_min_delta_f];
    brt.lon.pu  = [repmat(dyn.input_max_Pi, 8, 1); dyn.input_max_Pi_p; dyn.input_max_delta_e; dyn.input_max_delta_f];
    lat = yml.lateral;
    brt.lat = axis_setup(lat, {'v', 'p', 'r', 'phi'});
    brt.lat.prow   = [5; 10; 12; 7];
    brt.lat.U0_idx = [5:12, 2, 4];
    brt.lat.Dlo = [repmat(lat.input_min_Pi, 8, 1); lat.input_min_delta_a; lat.input_min_delta_r];
    brt.lat.Dhi = [repmat(lat.input_max_Pi, 8, 1); lat.input_max_delta_a; lat.input_max_delta_r];
    brt.lat.pl  = [repmat(dyn.input_min_Pi, 8, 1); dyn.input_min_delta_a; dyn.input_min_delta_r];
    brt.lat.pu  = [repmat(dyn.input_max_Pi, 8, 1); dyn.input_max_delta_a; dyn.input_max_delta_r];
end

function gridInfo = axis_setup(sec, names)
    gridInfo.gmin = cellfun(@(s) sec.(['grid_min_' s]), names)';
    gridInfo.gmax = cellfun(@(s) sec.(['grid_max_' s]), names)';
end

function yml = read_yml(path)
    lines = regexp(fileread(path), '\r?\n', 'split');
    yml = struct();  sec = '';
    for i = 1:numel(lines)
        ln = regexprep(lines{i}, '#.*$', '');
        if isempty(strtrim(ln)), continue; end
        tok = regexp(ln, '^(\s*)([^:\s][^:]*):\s*(.*)$', 'tokens', 'once');
        if isempty(tok), continue; end
        indent = numel(tok{1});  key = strtrim(tok{2});  val = strtrim(tok{3});
        if indent == 0
            if isempty(val), sec = key;  yml.(sec) = struct(); end
            continue;
        end
        if isempty(sec), continue; end
        num = str2double(val);
        if isnan(num), yml.(sec).(key) = val; else, yml.(sec).(key) = num; end
    end
end
