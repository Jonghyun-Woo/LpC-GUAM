% Model-mismatch DISTURBANCE ENVELOPE, per-axis method (see channel order below).
%
%   lon (trim-only):  e(dev) = f_nl(x0+dev, u_trim) - Ap*dev
%                     fit over inner REGION_FRAC=0.5 grid, plain OLS (no inflation).
%   lat (envelope):   Y(dev) = max_u |f_nl(x0+dev, u_trim+u) - (Ap*dev + Bp*u)|  (bang-bang)
%                     fit over the FULL grid, OLS + per-channel inflation to upper-bound Y.
%
% Why they differ: the lat r-channel residual grows toward the grid edge, so it needs a
% full-grid upper-bounding envelope; the lon inflation blows up 20-27x at extreme UH (q-chan)
% and lon was already a valid certificate, so lon is left trim-only.
%
%   e_c(z) = phi(z) * beta_c,  phi(z) = [1, z1..z4, z1^2..z4^2, pairwise]  (15 terms), z=dev/half

clear; clc;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

% ------------------------------- settings --------------------------------
stackDir    = fullfile(root, 'reachability_data', 'guam_timestack', 'BRT');
outDir      = fullfile(root, 'reachability_data', 'mc_verify');
UH_LIST     = 1:20;
WH_IDX      = 2;
NG_PER_DIM  = 7;       % uniform grid per dim: 7^4 = 2401 state points
N_NU        = 24;      % random input samples per deviation (lat bang-bang envelope)
DT          = 0.01;
FIT_FLOOR   = 1e-3;    % ignore near-zero fit when computing inflation scalar
REGION_FRAC = struct('lon', 0.5, 'lat', 1.0);  % lon: inner tube (trim-only); lat: full grid
% -------------------------------------------------------------------------

if ~exist(outDir, 'dir'), mkdir(outDir); end
gridFile  = @(axis_name) fullfile(stackDir, sprintf('BRT_%s_MAT', upper(axis_name)), ...
                sprintf('GUAM_%s_BRT_UH1_WH%d_stack.mat', upper(axis_name), WH_IDX));
brt       = struct('lon', brt_axis_config('lon', gridFile('lon')), ...
                   'lat', brt_axis_config('lat', gridFile('lat')));
trimTable = load(fullfile(root, 'controller', 'trim_table_Poly_ConcatVer4p0.mat'));
GUAM      = LpC_GUAM(Config('althold', struct('dt', DT)));

axes_ = struct('name', {'lon', 'lat'}, 'chan', {{'u','w','q','th'}, {'v','p','r','ph'}});
D = struct('lon', struct([]), 'lat', struct([]));

fprintf('lon = trim-only OLS,  lat = max-over-u envelope + inflation\n');
fprintf('r/q-chan per UH (infl=1 & low cover on lon is expected -- OLS is not an upper bound):\n\n');

for ai = 1:2
    axisName = axes_(ai).name;
    is_lon   = strcmp(axisName, 'lon');
    gridInfo = brt.(axisName);
    half     = (gridInfo.grid_max - gridInfo.grid_min) / 2;
    cen      = (gridInfo.grid_max + gridInfo.grid_min) / 2;
    prow     = gridInfo.state_rows;
    frac     = REGION_FRAC.(axisName);

    fprintf('== %s ==  channels [%s]\n', upper(axisName), strjoin(axes_(ai).chan, ' '));
    for ui = 1:numel(UH_LIST)
        uh  = UH_LIST(ui);
        XU0 = trimTable.XU0_interp(:, uh, WH_IDX);
        U0  = XU0(13:25);
        Ap  = trimTable.(sprintf('Ap_%s_interp', axisName))(:, :, uh, WH_IDX);
        Bp  = trimTable.(sprintf('Bp_%s_interp', axisName))(:, :, uh, WH_IDX);

        state0  = [0; 0; -100; XU0(1:3); XU0(10:12); XU0(4:6)];
        trimRef = state0(prow);

        % uniform grid over inner frac of each dim
        edges = arrayfun(@(lo,hi) linspace(lo,hi,NG_PER_DIM)', ...
                         cen - frac*half, cen + frac*half, 'UniformOutput', false);
        [G1,G2,G3,G4] = ndgrid(edges{1}, edges{2}, edges{3}, edges{4});
        devs    = [G1(:), G2(:), G3(:), G4(:)];
        N_DELTA = size(devs, 1);

        if is_lon
            % ----- trim-only: no input perturbation, residual at u_trim -----
            eng = U0(5:13);
            srf = [U0(1)-U0(2); U0(1)+U0(2); U0(3); U0(3); U0(4)];
            Y = zeros(N_DELTA, 4);
            parfor k = 1:N_DELTA
                st  = state0;  st(prow) = trimRef + devs(k,:)';
                fnl = GUAM.state_derivative(st, eng, srf);
                e   = fnl(prow) - Ap * devs(k,:)';
                Y(k,:) = abs(e)';
            end
        else
            % ----- envelope: Y = max_u |f_nl - (Ap*dev + Bp*u)| over bang-bang box -----
            sched      = struct('axis', axisName, 'U0', U0);
            [sched.lb, sched.ub] = LivenessFilter.input_bounds(U0, gridInfo);
            us         = sched.lb' + (sched.ub - sched.lb)' .* rand(N_NU, numel(sched.lb));
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
        end

        zs   = (devs - cen') ./ half';
        Phi  = design_quad(zs);
        beta = zeros(15, 4);  R2 = zeros(1, 4);  infl = ones(1, 4);  cover = zeros(1, 4);
        for c = 1:4
            b = Phi \ Y(:, c);
            res  = Y(:, c) - Phi * b;
            R2(c) = 1 - sum(res.^2) / max(sum((Y(:,c) - mean(Y(:,c))).^2), eps);
            if is_lon
                beta(:, c) = b;                              % OLS only (no upper bound)
            else
                F   = max(Phi * b, 0);                       % OLS shape, clipped >=0
                sel = F > FIT_FLOOR * max(F);
                infl(c) = max([1, max(Y(sel, c) ./ F(sel))]);% scalar so s*F >= Y
                beta(:, c) = infl(c) * b;
            end
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
half_lon = ((brt.lon.grid_max - brt.lon.grid_min) / 2)';
half_lat = ((brt.lat.grid_max - brt.lat.grid_min) / 2)';
feature  = 'phi(z)=[1, z1..z4, z1^2..z4^2, pairwise], z=dev/half; lon trim-only OLS, lat max_u envelope + inflation';
save(fullfile(outDir, 'guam_disturbance_quadfit.mat'), ...
     'beta_lon', 'beta_lat', 'half_lon', 'half_lat', 'infl_lon', 'infl_lat', 'UH_LIST', 'feature');
fprintf('saved -> %s\n', fullfile(outDir, 'guam_disturbance_quadfit.mat'));
fprintf('  channel order lon [u w q theta], lat [v p r phi];  lon OLS, lat inflated\n');


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
