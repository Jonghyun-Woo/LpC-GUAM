% diag_e3_full_guam_caseA_currentU_linear.m

% Purpose:
%   Compare filter-OFF baseline and filter-ON anchored linear CBF-like filter
%   under off-trim initial perturbations.
close all; clear all;

addpath(genpath('C:\Users\grape\OneDrive\CAU\AISL\LpC-GUAM\Refactoring'));
cd('C:\Users\grape\OneDrive\CAU\AISL\LpC-GUAM');

S = load('tables/trim/trim_table_Poly_ConcatVer4p0.mat');

cfg.M        = 1500;
cfg.whAnchor = S.WH(3);
cfg.gamma    = 5.0;
cfg.margin   = 0.1;
cfg.rateMode = 'linear';
cfg.tol      = 1e-6;
cfg.fdTol    = 1e-3;
cfg.divW     = 200;
cfg.makePlots = true;
cfg.showPlots = true; 
cfg.plotDir   = fullfile(pwd, 'diag_plots', 'full_guam_caseA_currentU_linear');

if cfg.makePlots && ~exist(cfg.plotDir, 'dir')
    mkdir(cfg.plotDir);
end

% This must match the RSLQR implementation.
% If RSLQR passes raw state(4) to the filter schedule, set false.
% If RSLQR clips current-u to trim-table UH range, set true.
cfg.clipUH   = true;

% perts = { 11, 0.0,        'q=0.0'; ...
%           11, 0.1,        'q=0.1'; ...
%           11, 0.2,        'q=0.2'; ...
%            8, deg2rad(5), 'th=5deg' };

perts = {11, 0.08,        'q=0.08'};

fprintf('\nFull GUAM Case A: current-u + WH3-anchor, linear rate\n');
fprintf('WH3 anchor = %.6f ft/s, gamma=%.2f, margin=%.3g, clipUH=%d\n\n', ...
        cfg.whAnchor, cfg.gamma, cfg.margin, cfg.clipUH);

fprintf('%-9s %-4s %-8s %-8s %-6s %-7s %-8s %-8s %-6s %-7s %-7s %-8s %-8s %-9s %-8s %-8s %s\n', ...
        'pert','mode','V0','Vmax','exitK','gridEx','nAct','actBef','1stAct', ...
        'cbf%','fd%','grid%','dVerrMx','VmisMax','w_fin','th_fin','verdict');

for p = 1:size(perts, 1)
    idx  = perts{p, 1};
    val  = perts{p, 2};
    name = perts{p, 3};

    for onoff = [false true]
        R = run_case(idx, val, onoff, cfg);

        if onoff
            md = 'ON';
        else
            md = 'OFF';
        end

        fprintf('%-9s %-4s %-8.3f %-8.3f %-6d %-7d %-8d %-8d %-6d %-7.1f %-7.1f %-8.1f %-8.2e %-9.2e %-8.1f %-8.1f %s\n', ...
                name, md, R.V0, R.Vmax, R.exitK, R.gridExitK, ...
                R.nAct, R.actBeforeExit, R.firstActK, ...
                100*R.cbfRate, 100*R.fdRate, 100*R.gridRate, ...
                R.dVerrMax, R.VmisMax, R.wfin, R.thfin, R.verdict);
        if cfg.makePlots
            plot_full_guam_trace(R.trace, cfg, name, md);
        end
    end
    fprintf('\n');
end

% -------------------------------------------------------------------------
function R = run_case(idx, val, onoff, cfg)
    g = LpC_GUAM('lon_brt_verify');
    f = g.controller.liveness_lon;

    % Set anchor for both OFF and ON so diagnostics and controller share config.
    g.controller.liveness_wh_anchor = cfg.whAnchor;

    if onoff
        f.mode        = 'blend';
        f.rate_mode   = cfg.rateMode;
        f.gamma       = cfg.gamma;
        f.live_margin = cfg.margin;
    else
        f.mode = 'off';
    end

    g.reset();
    f.reset_counters();

    % Apply off-trim initial perturbation.
    g.state(idx) = g.state(idx) + val;

    rt = g.refTraj;
    dt = g.simConfig.dt;

    % Initial diagnostic query.
    [Vinit, okInit, frameInit] = query_filter_frame(g, cfg);
    if okInit && frameInit.insideGrid
        V0 = Vinit;
    else
        V0 = NaN;
    end

    nu = g.controller.liveness_lon.spec.nu;
    trace = init_full_trace(cfg.M, nu);

    Vmax      = -inf;
    exitK     = 0;
    gridExitK = 0;
    div       = false;

    nAct          = 0;
    firstActK     = 0;
    actBeforeExit = 0;

    nQuery = 0;
    nOK    = 0;
    nGrid  = 0;
    nNoCov = 0;

    cbfNum = 0;
    cbfDen = 0;
    firstCBFfailK = 0;

    fdNum = 0;
    fdDen = 0;
    firstFDfailK = 0;

    dVerrMax = 0;
    dVerrSum = 0;
    dVerrDen = 0;

    VmisMax = 0;

    for k = 1:cfg.M
        ref.pos     = rt.pos(:, k);
        ref.vel     = rt.vel(:, k);
        ref.chi     = rt.chi(k);
        ref.chi_dot = rt.chidot(k);

        statePre = g.state(:);

        % Diagnostic query BEFORE control/plant step.
        % Must match actual filter frame:
        %   UH = current body-u
        %   WH = WH3 anchor
        [Vpre, okPre, frame] = query_filter_frame(g, cfg);

        nQuery = nQuery + 1;
        if okPre
            nOK = nOK + 1;
        else
            nNoCov = nNoCov + 1;
        end

        if frame.insideGrid
            nGrid = nGrid + 1;
        elseif gridExitK == 0
            gridExitK = k;
        end

        if okPre && frame.insideGrid
            Vmax = max(Vmax, Vpre);
            if Vpre > 0 && exitK == 0
                exitK = k;
            end
        end

        % Pre-step linear-rate diagnostics in the same filter frame.
        alpha = NaN;
        beta  = NaN(nu, 1);
        rhs   = NaN;
        b     = NaN;
        lb    = NaN(nu, 1);
        ub    = NaN(nu, 1);
        dV_min = NaN;
        u_min_rate = NaN(nu, 1);

        if okPre && frame.insideGrid
            alpha = frame.gradV(:)' * (frame.A * frame.xF(:));
            beta  = (frame.gradV(:)' * frame.B)';

            rhs = -cfg.gamma * (Vpre + cfg.margin);
            b   = rhs - alpha;

            [lb, ub] = LivenessFilter.input_bounds(frame.U0F, f.spec);

            u_min_rate = (beta >= 0) .* lb + (beta < 0) .* ub;
            dV_min = alpha + beta' * u_min_rate;
        end

        % Advance one full GUAM closed-loop step.
        % The filter is called inside controller.control().
        g.step(ref);

        statePost = g.state(:);
        li = f.last_info;

        % Pull input info from filter diagnostics.
        u_nom = getfield_default(li, 'u_nom', NaN(nu, 1));
        u     = getfield_default(li, 'u',     NaN(nu, 1));
        u0    = getfield_default(li, 'u0',    NaN(nu, 1));

        if any(isnan(u0)) && all(isfinite(u_nom)) && all(isfinite(lb))
            u0 = min(max(u_nom(:), lb), ub);
        end

        u_nom = u_nom(:);
        u0    = u0(:);
        u     = u(:);

        dV_raw = NaN;
        dV_box = NaN;
        dV_flt = NaN;

        rawBoxOK = NaN;
        rawCBFOK = NaN;
        boxCBFOK = NaN;
        feasLin  = NaN;

        if okPre && frame.insideGrid && all(isfinite(beta)) && all(isfinite(u_nom))
            dV_raw = alpha + beta' * u_nom;
            rawBoxOK = double(all(u_nom >= lb - 1e-12) && all(u_nom <= ub + 1e-12));
            rawCBFOK = double(beta' * u_nom <= b + cfg.tol);
        end

        if okPre && frame.insideGrid && all(isfinite(beta)) && all(isfinite(u0))
            dV_box = alpha + beta' * u0;
            boxCBFOK = double(beta' * u0 <= b + cfg.tol);
        end

        if okPre && frame.insideGrid && all(isfinite(beta)) && all(isfinite(u))
            dV_flt = alpha + beta' * u;
        end

        if okPre && frame.insideGrid && isfinite(dV_min)
            feasLin = double(dV_min <= rhs + cfg.tol);
        end

        if onoff && isfield(li, 'active') && li.active
            nAct = nAct + 1;

            if firstActK == 0
                firstActK = k;
            end

            if exitK == 0
                actBeforeExit = actBeforeExit + 1;
            end
        end

        % Consistency check:
        % f.last_info.V should match diagnostic Vpre if actual filter and
        % this script use the same current-u + WH3-anchor frame.
        if onoff && okPre && frame.insideGrid && isfield(li, 'ok') && li.ok
            VmisMax = max(VmisMax, abs(li.V - Vpre));
        end

        % Filter-reported CBF inequality check.
        if onoff && isfield(li, 'ok') && li.ok && isfinite(li.dVdt) && isfinite(li.rhs)
            cbfDen = cbfDen + 1;
            if li.dVdt <= li.rhs + cfg.tol
                cbfNum = cbfNum + 1;
            elseif firstCBFfailK == 0
                firstCBFfailK = k;
            end
        end

        % Fixed-frame finite-difference check using same pre-step UH/WH/X0.
        dVfd = NaN;
        rateErr = NaN;

        if onoff && okPre && frame.insideGrid && isfield(li, 'ok') && li.ok && isfinite(li.rhs)
            [VnextFixed, okNextFixed, insideNextFixed] = query_fixed_frame(g, frame);

            if okNextFixed && insideNextFixed
                dVfd = (VnextFixed - Vpre) / dt;
                fdDen = fdDen + 1;

                if dVfd <= li.rhs + cfg.fdTol
                    fdNum = fdNum + 1;
                elseif firstFDfailK == 0
                    firstFDfailK = k;
                end

                if isfinite(li.dVdt)
                    rateErr = dVfd - li.dVdt;
                    dVerrMax = max(dVerrMax, abs(rateErr));
                    dVerrSum = dVerrSum + rateErr;
                    dVerrDen = dVerrDen + 1;
                end
            end
        end

        % ----------------------------- trace log -----------------------------
        j = trace.n + 1;
        trace.n = j;

        trace.k(j) = k;
        trace.t(j) = (k - 1) * dt;

        trace.statePre(j, :)  = statePre(:)';
        trace.statePost(j, :) = statePost(:)';
        trace.refVel(j, :)    = ref.vel(:)';
        trace.refPos(j, :)    = ref.pos(:)';

        trace.xF(j, :)        = frame.xF(:)';
        trace.uhF(j)          = frame.uhF;
        trace.whF(j)          = frame.whF;

        trace.V(j)            = Vpre;
        trace.gradV(j, :)     = frame.gradV(:)';
        trace.gradNorm(j)     = norm(frame.gradV);
        trace.insideGrid(j)   = double(frame.insideGrid);
        trace.ok(j)           = double(okPre);

        trace.alpha(j)        = alpha;
        trace.beta(j, :)      = beta(:)';
        trace.betaNorm(j)     = norm(beta);
        trace.rhs(j)          = rhs;
        trace.b(j)            = b;

        trace.dV_raw(j)       = dV_raw;
        trace.dV_box(j)       = dV_box;
        trace.dV_flt(j)       = dV_flt;
        trace.dV_min(j)       = dV_min;
        trace.dV_info(j)      = getfield_default(li, 'dVdt', NaN);
        trace.dVfd(j)         = dVfd;
        trace.rateErr(j)      = rateErr;

        trace.slack_raw(j)    = rhs - dV_raw;
        trace.slack_box(j)    = rhs - dV_box;
        trace.slack_flt(j)    = rhs - dV_flt;
        trace.slack_min(j)    = rhs - dV_min;
        trace.slack_info(j)   = getfield_default(li, 'rhs', NaN) - getfield_default(li, 'dVdt', NaN);

        trace.u_nom(j, :)     = u_nom(:)';
        trace.u0(j, :)        = u0(:)';
        trace.u(j, :)         = u(:)';
        trace.lb(j, :)        = lb(:)';
        trace.ub(j, :)        = ub(:)';
        trace.u_min(j, :)     = u_min_rate(:)';

        trace.rawBoxOK(j)     = rawBoxOK;
        trace.rawCBFOK(j)     = rawCBFOK;
        trace.boxCBFOK(j)     = boxCBFOK;
        trace.feasLin(j)      = feasLin;

        trace.active(j)       = double(getfield_default(li, 'active', false));
        trace.du(j)           = getfield_default(li, 'du', norm(u - u0));
        trace.satClip(j)      = getfield_default(li, 'sat_clip', norm(u_nom - u0));
        trace.cmdChange(j)    = getfield_default(li, 'command_changed', norm(u - u_nom));
        trace.exitflag(j)     = getfield_default(li, 'quadprog_exitflag', NaN);
        trace.solver{j}       = getfield_default(li, 'solver', 'none');

        trace.w(j)            = statePre(6);
        trace.thetaDeg(j)     = rad2deg(statePre(8));
        trace.q(j)            = statePre(11);
        trace.uBody(j)        = statePre(4);
        trace.wPost(j)        = statePost(6);
        trace.thetaPostDeg(j) = rad2deg(statePost(8));
        trace.qPost(j)        = statePost(11);
        trace.uBodyPost(j)    = statePost(4);

        trace.gridExitFlag(j) = double(~frame.insideGrid);
        trace.brtExitFlag(j)  = double(okPre && frame.insideGrid && Vpre > 0);
        % ---------------------------------------------------------------------

        % Divergence check.
        if ~all(isfinite(g.state)) || abs(g.state(6)) > cfg.divW
            div = true;

            if all(isfinite(g.state))
                [Vpost, okPost, framePost] = query_filter_frame(g, cfg);

                if ~framePost.insideGrid && gridExitK == 0
                    gridExitK = k + 1;
                end

                if okPost && framePost.insideGrid
                    Vmax = max(Vmax, Vpost);
                    if Vpost > 0 && exitK == 0
                        exitK = k + 1;
                    end
                end
            end
            break;
        end
    end

    if isinf(Vmax)
        Vmax = NaN;
    end

    okRate   = safe_rate(nOK, nQuery);
    gridRate = safe_rate(nGrid, nQuery);
    cbfRate  = safe_rate(cbfNum, cbfDen);
    fdRate   = safe_rate(fdNum, fdDen);

    if dVerrDen == 0
        dVerrMean = NaN;
    else
        dVerrMean = dVerrSum / dVerrDen;
    end

    wfin  = g.state(6);
    thfin = rad2deg(g.state(8));

    verdict = classify_result(div, exitK, gridExitK, Vmax, V0, nOK, nNoCov);

    trace = trim_full_trace(trace);

    R = struct();
    R.V0       = V0;
    R.Vmax     = Vmax;
    R.exitK    = exitK;
    R.gridExitK = gridExitK;
    R.nAct     = nAct;
    R.firstActK = firstActK;
    R.actBeforeExit = actBeforeExit;
    R.firstCBFfailK = firstCBFfailK;
    R.firstFDfailK  = firstFDfailK;
    R.div      = div;
    R.wfin     = wfin;
    R.thfin    = thfin;
    R.okRate   = okRate;
    R.gridRate = gridRate;
    R.cbfRate  = cbfRate;
    R.fdRate   = fdRate;
    R.dVerrMax = dVerrMax;
    R.dVerrMean = dVerrMean;
    R.VmisMax  = VmisMax;
    R.verdict  = verdict;
    R.trace    = trace;
end

% -------------------------------------------------------------------------
function [V, ok, frame] = query_filter_frame(g, cfg)
    % Query the BRT in the same frame intended for the filter:
    %   UH = current body-u
    %   WH = WH3 safety anchor

    uhF = filter_sched_uh(g, cfg);
    whF = cfg.whAnchor;

    [X0F, U0F] = g.controller.interp_xu0(uhF, whF);
    AF = g.controller.interp_mtrx(g.controller.LON.Ap, uhF, whF);
    BF = g.controller.interp_mtrx(g.controller.LON.Bp, uhF, whF);

    xF = [g.state(4)  - X0F(1); ...
          g.state(6)  - X0F(3); ...
          g.state(11) - X0F(5); ...
          g.state(8)  - X0F(11)];

    [V, gradV, ok] = g.controller.liveness_lon.lut.query(xF, uhF, whF);

    spec = g.controller.liveness_lon.spec;
    insideGrid = all(xF(:) >= spec.grid_min(:) - 1e-12) && ...
                 all(xF(:) <= spec.grid_max(:) + 1e-12);

    frame = struct();
    frame.uhF = uhF;
    frame.whF = whF;
    frame.X0F = X0F;
    frame.U0F = U0F;
    frame.AF  = AF;
    frame.BF  = BF;

    % Short aliases used by run_case.
    frame.A = AF;
    frame.B = BF;

    frame.xF  = xF;
    frame.gradV = gradV;
    frame.insideGrid = insideGrid;
end

% -------------------------------------------------------------------------
function [V, ok, insideGrid] = query_fixed_frame(g, frame)
    % Query V after plant step using SAME pre-step UH/WH/X0 frame.

    X0F = frame.X0F;
    uhF = frame.uhF;
    whF = frame.whF;

    xF = [g.state(4)  - X0F(1); ...
          g.state(6)  - X0F(3); ...
          g.state(11) - X0F(5); ...
          g.state(8)  - X0F(11)];

    [V, ~, ok] = g.controller.liveness_lon.lut.query(xF, uhF, whF);

    spec = g.controller.liveness_lon.spec;
    insideGrid = all(xF(:) >= spec.grid_min(:) - 1e-12) && ...
                 all(xF(:) <= spec.grid_max(:) + 1e-12);
end

% -------------------------------------------------------------------------
function uhF = filter_sched_uh(g, cfg)
    % Current body-axis u used for the BRT/filter UH schedule.
    %
    % This must match the RSLQR anchor implementation exactly.

    uhRaw = g.state(4);

    if cfg.clipUH
        uhF = max(g.controller.UH(1), min(g.controller.UH(end), uhRaw));
    else
        uhF = uhRaw;
    end
end

% -------------------------------------------------------------------------
function r = safe_rate(num, den)
    if den == 0
        r = NaN;
    else
        r = num / den;
    end
end

% -------------------------------------------------------------------------
function verdict = classify_result(div, exitK, gridExitK, Vmax, V0, nOK, nNoCov)
    if nOK == 0
        verdict = 'NO-SCHEDULE-COVERAGE';
        return;
    end

    if isnan(V0)
        verdict = 'NO-INITIAL-BRT-VALUE';
        return;
    end

    if V0 > 0
        verdict = 'INIT-OUTSIDE-BRT';
        return;
    end

    if div
        if gridExitK > 0 && exitK > 0
            verdict = 'DIVERGED-after-GRID-and-BRT-EXIT';
        elseif gridExitK > 0
            verdict = 'DIVERGED-after-GRID-EXIT';
        elseif exitK > 0
            verdict = 'DIVERGED-after-BRT-EXIT';
        else
            verdict = 'DIVERGED';
        end
        return;
    end

    if gridExitK > 0
        if exitK > 0 || (~isnan(Vmax) && Vmax > 0)
            verdict = 'GRID-EXIT-and-BRT-EXIT-bounded';
        else
            verdict = 'GRID-EXIT-bounded';
        end
        return;
    end

    if exitK > 0 || (~isnan(Vmax) && Vmax > 0)
        verdict = 'BRT-EXIT-bounded';
        return;
    end

    if nNoCov > 0
        verdict = 'BRT-PRESERVED-partial-schedule';
    else
        verdict = 'BRT-PRESERVED';
    end
end

% -------------------------------------------------------------------------
function trace = init_full_trace(M, nu)
    trace.n = 0;

    trace.k = NaN(M, 1);
    trace.t = NaN(M, 1);

    trace.statePre  = NaN(M, 12);
    trace.statePost = NaN(M, 12);
    trace.refVel    = NaN(M, 3);
    trace.refPos    = NaN(M, 3);

    trace.xF = NaN(M, 4);
    trace.uhF = NaN(M, 1);
    trace.whF = NaN(M, 1);

    trace.V        = NaN(M, 1);
    trace.gradV    = NaN(M, 4);
    trace.gradNorm = NaN(M, 1);

    trace.insideGrid = NaN(M, 1);
    trace.ok         = NaN(M, 1);

    trace.alpha    = NaN(M, 1);
    trace.beta     = NaN(M, nu);
    trace.betaNorm = NaN(M, 1);

    trace.rhs = NaN(M, 1);
    trace.b   = NaN(M, 1);

    trace.dV_raw  = NaN(M, 1);
    trace.dV_box  = NaN(M, 1);
    trace.dV_flt  = NaN(M, 1);
    trace.dV_min  = NaN(M, 1);
    trace.dV_info = NaN(M, 1);
    trace.dVfd    = NaN(M, 1);
    trace.rateErr = NaN(M, 1);

    trace.slack_raw  = NaN(M, 1);
    trace.slack_box  = NaN(M, 1);
    trace.slack_flt  = NaN(M, 1);
    trace.slack_min  = NaN(M, 1);
    trace.slack_info = NaN(M, 1);

    trace.u_nom = NaN(M, nu);
    trace.u0    = NaN(M, nu);
    trace.u     = NaN(M, nu);
    trace.lb    = NaN(M, nu);
    trace.ub    = NaN(M, nu);
    trace.u_min = NaN(M, nu);

    trace.rawBoxOK = NaN(M, 1);
    trace.rawCBFOK = NaN(M, 1);
    trace.boxCBFOK = NaN(M, 1);
    trace.feasLin  = NaN(M, 1);

    trace.active    = NaN(M, 1);
    trace.du        = NaN(M, 1);
    trace.satClip   = NaN(M, 1);
    trace.cmdChange = NaN(M, 1);
    trace.exitflag  = NaN(M, 1);

    trace.solver = cell(M, 1);

    trace.w            = NaN(M, 1);
    trace.thetaDeg     = NaN(M, 1);
    trace.q            = NaN(M, 1);
    trace.uBody        = NaN(M, 1);
    trace.wPost        = NaN(M, 1);
    trace.thetaPostDeg = NaN(M, 1);
    trace.qPost        = NaN(M, 1);
    trace.uBodyPost    = NaN(M, 1);

    trace.gridExitFlag = NaN(M, 1);
    trace.brtExitFlag  = NaN(M, 1);
end

% -------------------------------------------------------------------------
function trace = trim_full_trace(trace)
    n = trace.n;
    fields = fieldnames(trace);

    for i = 1:numel(fields)
        fn = fields{i};
        val = trace.(fn);

        if strcmp(fn, 'n')
            continue;
        end

        if isnumeric(val) || islogical(val)
            if size(val, 1) >= n && numel(val) > 1
                trace.(fn) = val(1:n, :);
            end
        elseif iscell(val)
            if size(val, 1) >= n
                trace.(fn) = val(1:n, :);
            end
        end
    end
end

% -------------------------------------------------------------------------
function val = getfield_default(s, field, defaultVal)
    if isstruct(s) && isfield(s, field)
        val = s.(field);
    else
        val = defaultVal;
    end
end

% -------------------------------------------------------------------------
function plot_full_guam_trace(trace, cfg, pertName, modeName)
    if isempty(trace) || trace.n == 0
        return;
    end

    if cfg.showPlots
        figVisible = 'on';
    else
        figVisible = 'off';
    end

    tag = sanitize_filename(sprintf('%s_%s_gamma%.2f_margin%.3g', ...
        pertName, modeName, cfg.gamma, cfg.margin));

    k = trace.k;

    stateNames = {'du', 'dw', 'dq', 'dtheta'};
    inputNames = {'Pi1','Pi2','Pi3','Pi4','Pi5','Pi6','Pi7','Pi8','Pi9','de','df'};

    % ---------------------------------------------------------------------
    % Figure 1: BRT value, value rates, slack, flags
    % ---------------------------------------------------------------------
    fig1 = figure('Name', ['full_value_rate_' tag], ...
                  'Color', 'w', 'Visible', figVisible);
    tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    plot(k, trace.V, 'LineWidth', 1.2); hold on;
    yline(0, '--', 'V=0');
    yline(-cfg.margin, ':', 'V=-margin');
    xlabel('step k');
    ylabel('V');
    title('Scheduled-anchor BRT value');
    grid on;

    nexttile;
    plot(k, [trace.dV_raw, trace.dV_box, trace.dV_flt, trace.dV_min, trace.dV_info, trace.dVfd, trace.rhs], ...
         'LineWidth', 1.0);
    yline(0, '--');
    xlabel('step k');
    ylabel('rate');
    title('Value-function rates');
    legend({'dV raw nom','dV box nom','dV filtered','dV min-box','dV info','finite-diff dV','rhs'}, ...
           'Location', 'best');
    grid on;

    nexttile;
    plot(k, [trace.slack_raw, trace.slack_box, trace.slack_flt, trace.slack_min, trace.slack_info], ...
         'LineWidth', 1.0);
    yline(0, '--', 'constraint boundary');
    xlabel('step k');
    ylabel('rhs - dV');
    title('CBF slack');
    legend({'raw nominal','box nominal','filtered','min-rate','info'}, ...
           'Location', 'best');
    grid on;

    nexttile;
    plot(k, trace.xF, 'LineWidth', 1.0);
    xlabel('step k');
    ylabel('x_F');
    title('Anchor-frame state x_F = [du,dw,dq,dtheta]');
    legend(stateNames, 'Location', 'best');
    grid on;

    nexttile;
    plot(k, trace.rateErr, 'LineWidth', 1.0); hold on;
    yline(0, '--');
    xlabel('step k');
    ylabel('dVfd - dVdt');
    title('Full GUAM rate prediction error');
    grid on;

    nexttile;
    stairs(k, [trace.active, trace.rawBoxOK, trace.rawCBFOK, trace.boxCBFOK, trace.feasLin, trace.insideGrid], ...
           'LineWidth', 1.0);
    ylim([-0.1, 1.1]);
    xlabel('step k');
    ylabel('logical');
    title('Filter / feasibility / grid flags');
    legend({'active','raw box OK','raw CBF OK','box CBF OK','linear feasible','inside grid'}, ...
           'Location', 'best');
    grid on;

    saveas(fig1, fullfile(cfg.plotDir, [tag '_01_value_rates.png']));
    savefig(fig1, fullfile(cfg.plotDir, [tag '_01_value_rates.fig']));

    % ---------------------------------------------------------------------
    % Figure 2: full-state behavior
    % ---------------------------------------------------------------------
    fig2 = figure('Name', ['full_state_' tag], ...
                  'Color', 'w', 'Visible', figVisible);
    tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    plot(k, [trace.uBody, trace.w, trace.q, trace.thetaDeg], 'LineWidth', 1.0);
    xlabel('step k');
    ylabel('state');
    title('Full GUAM longitudinal states before step');
    legend({'u body [ft/s]','w [ft/s]','q [rad/s]','theta [deg]'}, ...
           'Location', 'best');
    grid on;

    nexttile;
    plot(k, [trace.uBodyPost, trace.wPost, trace.qPost, trace.thetaPostDeg], 'LineWidth', 1.0);
    xlabel('step k');
    ylabel('state');
    title('Full GUAM longitudinal states after step');
    legend({'u body post','w post','q post','theta post [deg]'}, ...
           'Location', 'best');
    grid on;

    nexttile;
    plot(k, trace.refVel, 'LineWidth', 1.0);
    xlabel('step k');
    ylabel('ref vel');
    title('Reference velocity');
    legend({'v_N ref','v_E ref','v_D ref'}, 'Location', 'best');
    grid on;

    saveas(fig2, fullfile(cfg.plotDir, [tag '_02_full_state.png']));
    savefig(fig2, fullfile(cfg.plotDir, [tag '_02_full_state.fig']));

    % ---------------------------------------------------------------------
    % Figure 3: value gradient and beta
    % ---------------------------------------------------------------------
    fig3 = figure('Name', ['full_grad_beta_' tag], ...
                  'Color', 'w', 'Visible', figVisible);
    tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    plot(k, trace.gradV, 'LineWidth', 1.0);
    xlabel('step k');
    ylabel('grad V');
    title('Value-function gradient components');
    legend({'dV/ddu','dV/ddw','dV/ddq','dV/ddtheta'}, ...
           'Location', 'best');
    grid on;

    nexttile;
    plot(k, [trace.gradNorm, trace.betaNorm, trace.alpha], 'LineWidth', 1.0);
    xlabel('step k');
    ylabel('value');
    title('Gradient norm, beta norm, alpha');
    legend({'||grad V||','||beta||','alpha'}, 'Location', 'best');
    grid on;

    nexttile;
    plot(k, trace.beta, 'LineWidth', 1.0);
    xlabel('step k');
    ylabel('beta_i');
    title('Input sensitivity beta = B^T gradV');
    legend(inputNames, 'Location', 'eastoutside');
    grid on;

    saveas(fig3, fullfile(cfg.plotDir, [tag '_03_grad_beta.png']));
    savefig(fig3, fullfile(cfg.plotDir, [tag '_03_grad_beta.fig']));

    % ---------------------------------------------------------------------
    % Figure 4: input modification summary
    % ---------------------------------------------------------------------
    fig4 = figure('Name', ['full_input_summary_' tag], ...
                  'Color', 'w', 'Visible', figVisible);
    tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    plot(k, [trace.satClip, trace.du, trace.cmdChange], 'LineWidth', 1.0);
    xlabel('step k');
    ylabel('norm');
    title('Input modification norms');
    legend({'||u_{nom}-u0|| clipping','||u-u0|| CBF move','||u-u_{nom}|| total'}, ...
           'Location', 'best');
    grid on;

    nexttile;
    plot(k, trace.u_nom - trace.u0, 'LineWidth', 1.0);
    xlabel('step k');
    ylabel('u_{nom} - u0');
    title('Raw nominal bound violation after clipping');
    legend(inputNames, 'Location', 'eastoutside');
    grid on;

    nexttile;
    plot(k, trace.u - trace.u0, 'LineWidth', 1.0);
    xlabel('step k');
    ylabel('u - u0');
    title('Pure CBF projection relative to clipped nominal');
    legend(inputNames, 'Location', 'eastoutside');
    grid on;

    saveas(fig4, fullfile(cfg.plotDir, [tag '_04_input_summary.png']));
    savefig(fig4, fullfile(cfg.plotDir, [tag '_04_input_summary.fig']));

    % ---------------------------------------------------------------------
    % Figure 5: per-channel input with bounds
    % ---------------------------------------------------------------------
    fig5 = figure('Name', ['full_input_channels_' tag], ...
                  'Color', 'w', 'Visible', figVisible);
    tiledlayout(4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    nu = size(trace.u, 2);

    for i = 1:nu
        nexttile;
        plot(k, trace.u_nom(:, i), '--', 'LineWidth', 0.9); hold on;
        plot(k, trace.u0(:, i), ':', 'LineWidth', 1.1);
        plot(k, trace.u(:, i), '-', 'LineWidth', 1.2);
        plot(k, trace.lb(:, i), 'k--', 'LineWidth', 0.8);
        plot(k, trace.ub(:, i), 'k--', 'LineWidth', 0.8);

        xlabel('k');
        ylabel(inputNames{i});
        title(inputNames{i});
        grid on;

        if i == 1
            legend({'u_{nom}','u0 clipped','u filtered','lb','ub'}, ...
                   'Location', 'best');
        end
    end

    saveas(fig5, fullfile(cfg.plotDir, [tag '_05_input_channels.png']));
    savefig(fig5, fullfile(cfg.plotDir, [tag '_05_input_channels.fig']));

    % ---------------------------------------------------------------------
    % Figure 6: exit indicators and scheduling
    % ---------------------------------------------------------------------
    fig6 = figure('Name', ['full_exit_schedule_' tag], ...
                  'Color', 'w', 'Visible', figVisible);
    tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    stairs(k, [trace.gridExitFlag, trace.brtExitFlag, trace.insideGrid, trace.ok], ...
           'LineWidth', 1.0);
    ylim([-0.1, 1.1]);
    xlabel('step k');
    ylabel('logical');
    title('Grid/BRT coverage flags');
    legend({'grid exit','BRT exit','inside grid','LUT ok'}, 'Location', 'best');
    grid on;

    nexttile;
    plot(k, [trace.uhF, trace.whF], 'LineWidth', 1.0);
    xlabel('step k');
    ylabel('schedule');
    title('BRT scheduling variables');
    legend({'UH current-u schedule','WH anchor'}, 'Location', 'best');
    grid on;

    nexttile;
    plot(k, trace.exitflag, 'LineWidth', 1.0);
    xlabel('step k');
    ylabel('quadprog exitflag');
    title('QP exit flag');
    grid on;

    saveas(fig6, fullfile(cfg.plotDir, [tag '_06_exit_schedule.png']));
    savefig(fig6, fullfile(cfg.plotDir, [tag '_06_exit_schedule.fig']));
end

% -------------------------------------------------------------------------
function s = sanitize_filename(s)
    s = regexprep(s, '[^\w\-.]', '_');
end
