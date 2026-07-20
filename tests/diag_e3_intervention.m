% diag_e3_full_guam_caseA_currentU_linear.m

% Purpose:
%   Compare filter-OFF baseline and filter-ON anchored linear CBF-like filter
%   under off-trim initial perturbations.
close all; clear all;

% cd('C:\Users\grape\OneDrive\CAU\AISL\LpC-GUAM');
addpath(genpath(pwd));

% Trim Table setting
S = load('trim_table_Poly_ConcatVer4p0.mat');
% Environment Parameters
cfg.M        = 10000;        % Total simulation steps
cfg.whAnchor = S.WH(3);     % WH for calling trim table in perspective of BRT
cfg.target_vel = 150.0;       % Target velocity related to reference generator
cfg.gamma    = 5.0;         % Sensitivity of the Class-K function
cfg.margin   = 0.001;         % CBF margin for zero-hold effect
cfg.decimation = 1;         % Reference gervenor decimation
cfg.makePlots = true;       % Save plot  
cfg.showPlots = true;       % Show plots
cfg.plotDir   = fullfile(pwd, 'diag_plots', 'full_guam_caseA_currentU_linear');

if cfg.makePlots && ~exist(cfg.plotDir, 'dir')
    mkdir(cfg.plotDir);
end

% This must match the RSLQR implementation.
% If RSLQR passes raw state(4) to the filter schedule, set false.
% If RSLQR clips current-u to trim-table UH range, set true.
cfg.clipUH   = true;

% Initial perturbation condition (off-trim state)
perts = {11, 0.0,        'q=0.0'};

fprintf('\nFull GUAM Case A: current-u + WH3-anchor, linear rate\n');
fprintf('WH3 anchor = %.6f ft/s, gamma=%.2f, margin=%.3g, clipUH=%d\n\n', cfg.whAnchor, cfg.gamma, cfg.margin, cfg.clipUH);

fprintf('%-9s %-4s %-8s %-8s %-6s %-7s %-8s %-8s %-9s %s\n', ...
        'pert','mode','V0','Vmax','exitK','gridEx','nAct','cbf%','dVerrMx','verdict');

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
        % Runtime logging
        fprintf('%-9s %-4s %-8.3f %-8.3f %-6d %-7d %-8d %-8.1f %-9.2e %s\n', ...
                name, md, R.V0, R.Vmax, R.exitK, R.gridExitK, R.nAct, 100*R.cbfRate, R.dVerrMax, R.verdict);
        if cfg.makePlots
            % Save log file
            plot_full_guam_trace(R.trace, cfg, name, md);
            traceFile = fullfile(cfg.plotDir, ...
                sanitize_filename(sprintf('trace_%s_%s_ref_%.1f_dec_%d.mat', name, md, cfg.target_vel, cfg.decimation)));
            save(traceFile, 'R', 'cfg');
        end
    end
    fprintf('\n');
end

% -------------------------------------------------------------------------
function R = run_case(idx, val, onoff, cfg)
    % Build the config hub from the diagnostic struct's mission fields; the
    % rest of `cfg` (gamma/margin/whAnchor/tol/divW/...) stays a plain struct
    % used only inside diagnostics.
    hub = Config('lon_brt_verify', struct('steps', cfg.M, 'target_vel', cfg.target_vel));
    c   = Controller(hub.controller, hub.sim.dt);
    g   = LpC_GUAM(hub);
    f   = c.safety_filter;

    % Config defaults used only inside diagnostics.
    tol   = getfield_default(cfg, 'tol',   1e-6);
    divW  = getfield_default(cfg, 'divW',  200);

    % Set anchor for both OFF and ON so diagnostics and controller share config.
    c.safety_filter_wh_anchor = cfg.whAnchor;

    if onoff
        f.mode        = 'blend';
        f.gamma       = cfg.gamma;
        f.live_margin = cfg.margin;
    else
        f.mode = 'off';
    end

    rt = hub.controller.getReferenceTrajectory();
    [x0, engine0, surface0] = c.initial_condition(rt);
    g.reset(x0, engine0, surface0);
    c.reset();
    f.reset_counters();

    % Apply off-trim initial perturbation.
    g.state(idx) = g.state(idx) + val;

    dt = g.simConfig.dt;

    nu = c.safety_filter.spec.nu;

    trace = init_full_trace(cfg.M, nu);
    trace = add_transition_trace_fields(trace, cfg.M);

    % Initial transition-envelope query.
    [envInit, ~, ~] = query_transition_frame(g, c, cfg);
    if envInit.ok
        V0 = envInit.V;
    else
        V0 = NaN;
    end

    % Primary metrics: transition-envelope basis.
    Vmax      = -inf;
    exitK     = 0;
    gridExitK = 0;
    div       = false;

    % Filter intervention.
    nAct      = 0;
    firstActK = 0;

    % Envelope coverage.
    nQuery = 0;
    nOK    = 0;
    nNoCov = 0;

    % CBF-model consistency.
    cbfNum = 0;
    cbfDen = 0;

    % Full GUAM finite-difference consistency.
    dVerrMax = 0;
    dVerrDen = 0;

    for k = 1:cfg.M
        if mod(k-1, cfg.decimation) == 0
            ref.pos     = rt.pos(:, k);
            ref.vel     = rt.vel(:, k);
            ref.chi     = rt.chi(k);
            ref.chi_dot = rt.chidot(k);
        end

        statePre = g.state(:);

        % -----------------------------------------------------------------
        % Transition-BRT diagnostic query BEFORE control/plant step.
        %
        % envPre:
        %   primary metric frame, Venv = min(valid Vcur, valid Vnext)
        %
        % selPre:
        %   selected frame that should match the filter-side BRT target.
        % -----------------------------------------------------------------
        [envPre, selPre, candPre] = query_transition_frame(g, c, cfg);

        nQuery = nQuery + 1;

        if envPre.ok
            nOK = nOK + 1;

            Vmax = max(Vmax, envPre.V);

            if envPre.V > 0 && exitK == 0
                exitK = k;
            end
        else
            nNoCov = nNoCov + 1;

            if gridExitK == 0
                gridExitK = k;
            end
        end

        % -----------------------------------------------------------------
        % Selected-frame diagnostics.
        % CBF/rate checks must use the selected BRT frame, not always current.
        % -----------------------------------------------------------------
        alpha = NaN;
        beta  = NaN(nu, 1);
        rhs   = NaN;
        b     = NaN;
        lb    = NaN(nu, 1);
        ub    = NaN(nu, 1);

        if selPre.valid
            alpha = selPre.gradV(:)' * (selPre.A * selPre.xF(:));
            beta  = (selPre.gradV(:)' * selPre.B)';

            rhs = -cfg.gamma * (selPre.V + cfg.margin);
            b   = rhs - alpha;

            [lb, ub] = LivenessFilter.input_bounds(selPre.U0F, f.spec);
        end

        % -----------------------------------------------------------------
        % Advance one full GUAM closed-loop step.
        % The liveness filter is called inside controller.control().
        % -----------------------------------------------------------------
        [eng_cmd, srf_cmd] = c.control(g.state, ref);
        g.step(eng_cmd, srf_cmd);

        statePost = g.state(:);
        li = f.last_info;

        % Pull input info from filter diagnostics.
        u_nom = getfield_default(li, 'u_nom', NaN(nu, 1));
        u0    = getfield_default(li, 'u0',    NaN(nu, 1));
        u     = getfield_default(li, 'u',     NaN(nu, 1));

        u_nom = u_nom(:);
        u0    = u0(:);
        u     = u(:);

        if any(isnan(u0)) && all(isfinite(u_nom)) && all(isfinite(lb))
            u0 = min(max(u_nom, lb), ub);
        end

        % Selected-frame model rate after actual filter output.
        dV_flt = NaN;
        if selPre.valid && all(isfinite(beta)) && all(isfinite(u))
            dV_flt = alpha + beta' * u;
        end

        % Actual filter intervention count.
        if onoff && isfield(li, 'active') && li.active
            nAct = nAct + 1;

            if firstActK == 0
                firstActK = k;
            end
        end

        % Filter-reported CBF inequality check.
        % This uses filter-side li.dVdt and li.rhs. Since the filter now
        % performs BRT transition internally, these should already correspond
        % to the selected BRT frame.
        if onoff && isfield(li, 'ok') && li.ok && ...
                isfield(li, 'dVdt') && isfield(li, 'rhs') && ...
                isfinite(li.dVdt) && isfinite(li.rhs)

            cbfDen = cbfDen + 1;

            if li.dVdt <= li.rhs + tol
                cbfNum = cbfNum + 1;
            end
        end

        % Full GUAM finite-difference check in the selected pre-step frame.
        dVfd = NaN;
        rateErr = NaN;

        if onoff && selPre.valid && ...
                isfield(li, 'ok') && li.ok && ...
                isfield(li, 'rhs') && isfinite(li.rhs)

            [VnextFixed, okNextFixed, insideNextFixed] = query_fixed_candidate(g, c, selPre);

            if okNextFixed && insideNextFixed
                dVfd = (VnextFixed - selPre.V) / dt;

                if isfield(li, 'dVdt') && isfinite(li.dVdt)
                    rateErr = dVfd - li.dVdt;
                    dVerrMax = max(dVerrMax, abs(rateErr));
                    dVerrDen = dVerrDen + 1;
                end
            end
        end

        % -----------------------------------------------------------------
        % Trace log.
        % Keep trace.V as primary transition-envelope value.
        % Keep trace.xF / uhF / whF as selected-frame state.
        % -----------------------------------------------------------------
        j = trace.n + 1;
        trace.n = j;

        trace.k(j) = k;
        trace.t(j) = (k - 1) * dt;

        trace.statePre(j, :)  = statePre(:)';
        trace.statePost(j, :) = statePost(:)';
        trace.refVel(j, :)    = ref.vel(:)';
        trace.refPos(j, :)    = ref.pos(:)';

        % Primary envelope values.
        trace.V(j)      = envPre.V;
        trace.ok(j)     = double(envPre.ok);
        trace.insideGrid(j) = double(envPre.ok);

        trace.Vcur(j)   = candPre.current.V;
        trace.Vnext(j)  = candPre.next.V;
        trace.Venv(j)   = envPre.V;
        trace.Vsel(j)   = selPre.V;

        trace.curIdx(j)  = candPre.current.idx;
        trace.nextIdx(j) = candPre.next.idx;
        trace.selIdx(j)  = selPre.idx;

        trace.curValid(j)  = double(candPre.current.valid);
        trace.nextValid(j) = double(candPre.next.valid);
        trace.envValid(j)  = double(envPre.ok);
        trace.selValid(j)  = double(selPre.valid);

        trace.selectedIsNext(j)  = double(strcmp(selPre.name, 'next'));
        trace.transitionReady(j) = double(candPre.next.valid && candPre.next.V < 0);

        % Selected BRT frame.
        trace.xF(j, :)    = selPre.xF(:)';
        trace.uhF(j)      = selPre.uhF;
        trace.whF(j)      = selPre.whF;
        trace.gradV(j, :) = selPre.gradV(:)';
        trace.gradNorm(j) = norm(selPre.gradV);

        trace.xCur(j, :)  = candPre.current.xF(:)';
        trace.xNext(j, :) = candPre.next.xF(:)';
        trace.xSel(j, :)  = selPre.xF(:)';

        trace.alpha(j)    = alpha;
        trace.beta(j, :)  = beta(:)';
        trace.betaNorm(j) = norm(beta);
        trace.rhs(j)      = rhs;
        trace.b(j)        = b;

        trace.dV_flt(j)  = dV_flt;
        trace.dV_info(j) = getfield_default(li, 'dVdt', NaN);
        trace.dVfd(j)    = dVfd;
        trace.rateErr(j) = rateErr;

        trace.slack_flt(j)  = rhs - dV_flt;
        trace.slack_info(j) = getfield_default(li, 'rhs', NaN) - ...
                              getfield_default(li, 'dVdt', NaN);

        trace.u_nom(j, :) = u_nom(:)';
        trace.u0(j, :)    = u0(:)';
        trace.u(j, :)     = u(:)';
        trace.lb(j, :)    = lb(:)';
        trace.ub(j, :)    = ub(:)';

        trace.active(j)    = double(getfield_default(li, 'active', false));
        trace.du(j)        = getfield_default(li, 'du', norm(u - u0));
        trace.satClip(j)   = getfield_default(li, 'sat_clip', norm(u_nom - u0));
        trace.cmdChange(j) = getfield_default(li, 'command_changed', norm(u - u_nom));
        trace.exitflag(j)  = getfield_default(li, 'quadprog_exitflag', NaN);
        trace.solver{j}    = getfield_default(li, 'solver', 'none');

        trace.w(j)            = statePre(6);
        trace.thetaDeg(j)     = rad2deg(statePre(8));
        trace.q(j)            = statePre(11);
        trace.uBody(j)        = statePre(4);

        trace.wPost(j)        = statePost(6);
        trace.thetaPostDeg(j) = rad2deg(statePost(8));
        trace.qPost(j)        = statePost(11);
        trace.uBodyPost(j)    = statePost(4);

        trace.gridExitFlag(j) = double(~envPre.ok);
        trace.brtExitFlag(j)  = double(envPre.ok && envPre.V > 0);
        trace.envExitFlag(j)  = double(envPre.ok && envPre.V > 0);
        trace.selExitFlag(j)  = double(selPre.valid && selPre.V > 0);

        % -----------------------------------------------------------------
        % Divergence check.
        % -----------------------------------------------------------------
        if ~all(isfinite(g.state)) || abs(g.state(6)) > divW
            div = true;

            if all(isfinite(g.state))
                [envPost, ~, ~] = query_transition_frame(g, c, cfg);

                if ~envPost.ok && gridExitK == 0
                    gridExitK = k + 1;
                end

                if envPost.ok
                    Vmax = max(Vmax, envPost.V);

                    if envPost.V > 0 && exitK == 0
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

    cbfRate = safe_rate(cbfNum, cbfDen);

    wfin  = g.state(6);
    thfin = rad2deg(g.state(8));

    verdict = classify_result(div, exitK, gridExitK, Vmax, V0, nOK, nNoCov);

    trace = trim_full_trace(trace);

    R = struct();
    R.V0        = V0;
    R.Vmax      = Vmax;
    R.exitK     = exitK;
    R.gridExitK = gridExitK;

    R.nAct      = nAct;
    R.firstActK = firstActK;

    R.div       = div;
    R.wfin      = wfin;
    R.thfin     = thfin;

    R.okRate    = safe_rate(nOK, nQuery);
    R.cbfRate   = cbfRate;
    R.dVerrMax  = dVerrMax;

    R.verdict   = verdict;
    R.trace     = trace;
end

% -------------------------------------------------------------------------
function [env, selected, cand] = query_transition_frame(g, c, cfg)
    cand.current = make_brt_candidate(g, c, cfg, 'current');
    cand.next    = make_brt_candidate(g, c, cfg, 'next');

    c = cand.current;
    n = cand.next;

    vals  = [];
    names = {};

    if c.valid
        vals(end+1) = c.V; %#ok<AGROW>
        names{end+1} = 'current'; %#ok<AGROW>
    end

    if n.valid
        vals(end+1) = n.V; %#ok<AGROW>
        names{end+1} = 'next'; %#ok<AGROW>
    end

    env = struct();

    if isempty(vals)
        env.V = NaN;
        env.ok = false;
        env.insideGrid = false;
        env.best = 'none';
    else
        [env.V, ii] = min(vals);
        env.ok = true;
        env.insideGrid = true;
        env.best = names{ii};
    end

    % This mirrors the intended transition logic:
    %   if next BRT is valid and already contains the state, move to next.
    %   otherwise keep current.
    %
    % If current and next are the same final index, selected = next/current
    % is equivalent, but we keep next for consistency with filter code.
    if c.idx == n.idx
        selected = n;
    elseif n.valid && n.V < 0
        selected = n;
    else
        selected = c;
    end
end

% -------------------------------------------------------------------------
function brt = make_brt_candidate(g, c, cfg, which)
    bc   = c.baseline_controller;
    sf   = c.safety_filter;
    spec = sf.spec;

    uhRaw = g.state(4);

    if cfg.clipUH
        uhRaw = max(bc.UH(1), min(bc.UH(end), uhRaw));
    end

    [~, current_idx] = min(abs(bc.UH(:) - uhRaw));
    next_idx = min(numel(bc.UH), current_idx + 1);

    switch lower(which)
        case 'current'
            idx = current_idx;
        case 'next'
            idx = next_idx;
        otherwise
            error('Unknown BRT candidate type: %s', which);
    end

    uh = bc.UH(idx);
    wh = cfg.whAnchor;

    [X0, U0] = bc.interp_xu0(uh, wh);
    A = bc.interp_mtrx(bc.LON.Ap, uh, wh);
    B = bc.interp_mtrx(bc.LON.Bp, uh, wh);

    x = [g.state(4)  - X0(1); ...
         g.state(6)  - X0(3); ...
         g.state(11) - X0(5); ...
         g.state(8)  - X0(11)];

    insideGrid = all(x(:) >= spec.grid_min(:) - 1e-12) && ...
                 all(x(:) <= spec.grid_max(:) + 1e-12);

    [V, gradV, ok] = sf.lut.query(x, uh, wh);

    if isempty(gradV) || numel(gradV) ~= 4
        gradV = NaN(4, 1);
    end

    brt = struct();
    brt.name = lower(which);

    brt.idx  = idx;
    brt.uhF  = uh;
    brt.whF  = wh;

    brt.X0F  = X0;
    brt.U0F  = U0;
    brt.A    = A;
    brt.B    = B;

    brt.xF    = x;
    brt.V     = V;
    brt.gradV = gradV(:);

    brt.ok         = ok;
    brt.insideGrid = insideGrid;
    brt.valid      = ok && insideGrid && isfinite(V);
end

% -------------------------------------------------------------------------
function [V, ok, insideGrid] = query_fixed_candidate(g, c, brt)
    X0F = brt.X0F;
    uhF = brt.uhF;
    whF = brt.whF;

    xF = [g.state(4)  - X0F(1); ...
          g.state(6)  - X0F(3); ...
          g.state(11) - X0F(5); ...
          g.state(8)  - X0F(11)];

    [V, ~, ok] = c.safety_filter.lut.query(xF, uhF, whF);

    spec = c.safety_filter.spec;
    insideGrid = all(xF(:) >= spec.grid_min(:) - 1e-12) && ...
                 all(xF(:) <= spec.grid_max(:) + 1e-12);
end

% -------------------------------------------------------------------------
function trace = add_transition_trace_fields(trace, M)
    trace.Vcur  = NaN(M, 1);
    trace.Vnext = NaN(M, 1);
    trace.Venv  = NaN(M, 1);
    trace.Vsel  = NaN(M, 1);

    trace.curIdx  = NaN(M, 1);
    trace.nextIdx = NaN(M, 1);
    trace.selIdx  = NaN(M, 1);

    trace.curValid  = NaN(M, 1);
    trace.nextValid = NaN(M, 1);
    trace.envValid  = NaN(M, 1);
    trace.selValid  = NaN(M, 1);

    trace.selectedIsNext  = NaN(M, 1);
    trace.transitionReady = NaN(M, 1);

    trace.xCur  = NaN(M, 4);
    trace.xNext = NaN(M, 4);
    trace.xSel  = NaN(M, 4);

    trace.envExitFlag = NaN(M, 1);
    trace.selExitFlag = NaN(M, 1);
end

% -------------------------------------------------------------------------
function verdict = classify_result(div, exitK, gridExitK, Vmax, V0, nOK, nNoCov)
    if nOK == 0
        verdict = 'NO-TRANSITION-BRT-COVERAGE';
        return;
    end

    if isnan(V0)
        verdict = 'NO-INITIAL-TRANSITION-BRT-VALUE';
        return;
    end

    if V0 > 0
        verdict = 'INIT-OUTSIDE-TRANSITION-BRT';
        return;
    end

    if div
        if gridExitK > 0 && exitK > 0
            verdict = 'DIVERGED-after-TRANSITION-GRID-and-BRT-EXIT';
        elseif gridExitK > 0
            verdict = 'DIVERGED-after-TRANSITION-GRID-EXIT';
        elseif exitK > 0
            verdict = 'DIVERGED-after-TRANSITION-BRT-EXIT';
        else
            verdict = 'DIVERGED-inside-TRANSITION-BRT';
        end
        return;
    end

    if gridExitK > 0
        if exitK > 0 || (~isnan(Vmax) && Vmax > 0)
            verdict = 'TRANSITION-GRID-and-BRT-EXIT-bounded';
        else
            verdict = 'TRANSITION-GRID-EXIT-bounded';
        end
        return;
    end

    if exitK > 0 || (~isnan(Vmax) && Vmax > 0)
        verdict = 'TRANSITION-BRT-EXIT-bounded';
        return;
    end

    if nNoCov > 0
        verdict = 'TRANSITION-BRT-PRESERVED-partial-coverage';
    else
        verdict = 'TRANSITION-BRT-PRESERVED';
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
end

% -------------------------------------------------------------------------
function s = sanitize_filename(s)
    s = regexprep(s, '[^\w\-.]', '_');
end
