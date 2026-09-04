% check_trimref_liveness - build the trim-schedule reference trajectory and
% verify the LON BRT / CBF liveness prediction along it.
%
% For every sample the LON state [u;w;q;theta] is expressed as a perturbation
% from (a) the current segment anchor trim k and (b) the next anchor k+1, and
% the corresponding BRT value functions are interpolated. V<0 means the point
% is inside that BRT (the target set of that anchor is reachable within the
% BRT horizon under bounded perturbation inputs).
%
% Prints per-segment statistics (entry/exit progress fractions, handover
% window) and draws the three verification figures (plot_trimref_brt).
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

%% 1) Reference trajectory (20 trims @ WH3, 5 s holds, uniform T_seg = 2 s)
% NOTE: V along the path is a geometric quantity - changing T_seg rescales
% only the TIME spent outside a BRT, not the entry/exit progress fractions.
opts = struct();
opts.T_seg = 2.0;            % uniform segment time [s] (or opts.accel = 4.0)
traj = TrimScheduleTrajectory.build(opts);
fprintf('trajectory: N = %d samples, t = [0, %.2f] s, %d segments\n', ...
        numel(traj.time), traj.time(end), traj.n_trim - 1);

%% 2) Load LON BRT tables (data/GUAM_LON_BRT_HJIR_UH{k}_WH3.mat)
spec = FilterConfig.channelSpec('lon');
gv = cell(1, 4);
for d = 1:4
    gv{d} = linspace(spec.grid_min(d), spec.grid_max(d), spec.grid_num(d));
end
brt = cell(1, traj.n_trim);
for k = 1:traj.n_trim
    f = fullfile(root, 'data', sprintf('GUAM_LON_BRT_HJIR_UH%d_WH%d.mat', ...
                                       k, traj.wh_idx));
    S = load(f, 'data');
    brt{k} = S.data;
end

%% 3) Evaluate V along the trajectory against current / next anchors
N = numel(traj.time);
ev = struct('V_cur', nan(1, N), 'V_next', nan(1, N));
for i = 1:N
    k  = traj.seg(i);
    kn = min(k + 1, traj.n_trim);
    ev.V_cur(i)  = brt_value(brt{k},  gv, traj.lon(:, i) - traj.trim_lon(:, k));
    ev.V_next(i) = brt_value(brt{kn}, gv, traj.lon(:, i) - traj.trim_lon(:, kn));
end

%% 4) Report
gap = ev.V_cur >= 0 & ev.V_next >= 0;
fprintf('\nsamples outside current BRT (V_cur>=0) : %d / %d\n', nnz(ev.V_cur >= 0), N);
fprintf('samples outside next    BRT (V_next>=0): %d / %d\n', nnz(ev.V_next >= 0), N);
fprintf('samples in NEITHER BRT (coverage gap)  : %d / %d\n', nnz(gap), N);

fprintf('\n k   UH_k -> UH_k+1   T_seg   enter-next@   leave-cur@   handover window\n');
for k = 1:traj.n_trim - 1
    m  = traj.seg == k & traj.progress < k + 1;
    fr = traj.progress(m) - k;
    vc = ev.V_cur(m);  vn = ev.V_next(m);
    i_in  = find(vn < 0, 1, 'first');
    i_out = find(vc >= 0, 1, 'first');
    f_in  = fr(max(i_in, 1));  if isempty(i_in),  f_in  = NaN; end
    f_out = NaN;               if ~isempty(i_out), f_out = fr(i_out); end
    fprintf('%2d  %6.1f -> %6.1f  %5.2f s     %5.1f %%       %5.1f %%     %5.1f %%\n', ...
            k, traj.UH(k), traj.UH(k + 1), traj.T_seg(k), ...
            100 * f_in, 100 * f_out, 100 * (min(f_out, 1) - f_in));
end

if any(gap)
    warning('check_trimref_liveness:gap', ...
        'Reference leaves BOTH BRTs in segments: %s', ...
        mat2str(unique(traj.seg(gap))));
else
    fprintf('\nOK: every sample lies inside at least one BRT (no coverage gap).\n');
end

%% 5) Closed-loop GUAM run on the same mission (nominal RSLQR, filter off)
% Flies the 'trim_schedule' scenario through the full LpC_GUAM plant so the
% figures can overlay how the vehicle actually moved (red) against the
% reference (black) and the BRT tube.
params = struct('filter_mode', 'off');
cfg  = Config('trim_schedule', params);
guam = LpC_GUAM(cfg);
rt   = guam.refTraj;
M    = size(rt.pos, 2);
st   = zeros(12, M);
guam.reset();
for i = 1:M
    ref = struct('pos', rt.pos(:, i), 'vel', rt.vel(:, i), ...
                 'chi', rt.chi(i), 'chi_dot', rt.chidot(i));
    st(:, i) = guam.state;
    guam.step(ref);
end
guam_trace = struct('time', rt.time, 'lon', st([4 6 11 8], :));
fprintf('\nclosed-loop run: %d steps, final u = %.1f ft/s (ref %.1f)\n', ...
        M, st(4, end), rt.vel(1, end));

%% 6) Figures
figs = plot_trimref_brt(traj, ev, brt, gv, struct('guam', guam_trace));   %#ok<NASGU>

%% 7) (u, t, theta) tube zoom for ONE chosen segment (edit tube_seg, run this section)
% BRT_k / BRT_{k+1} zero-contour slices + reference (black) + flown GUAM
% (red) over that segment's scheduled time window.
tube_seg = 13;                     % <- segment to inspect (1..19)
plot_trimref_tube_zoom(traj, brt, gv, tube_seg, struct('guam', guam_trace));

%% 8) Handover slices for ONE chosen segment (edit hand_seg, run this section)
hand_seg = 1;                     % <- segment to inspect (1..19)
plot_trimref_handover(traj, ev, brt, gv, hand_seg, guam_trace);

%% 9) CBF risk detection along the FLOWN trajectory
% The runtime CBF quantity is V evaluated at the vehicle's OWN state
% (its own (w,q) included) against the scheduled anchors - exactly what
% LivenessFilter.transition_ready = (V_next < 0) uses. A segment is
% "detected as risky" when the vehicle is NOT inside BRT_{k+1} at the
% scheduled handover time t_node(k+1).
tg = guam_trace.time;
Mg = numel(tg);
evg = struct('V_cur', nan(1, Mg), 'V_next', nan(1, Mg), 'seg', zeros(1, Mg));
for i = 1:Mg
    k = find(traj.t_node <= tg(i) + 1e-9, 1, 'last');
    if isempty(k), k = 1; end
    k = min(k, traj.n_trim - 1);
    evg.seg(i)    = k;
    evg.V_cur(i)  = brt_value(brt{k},   gv, guam_trace.lon(:, i) - traj.trim_lon(:, k));
    evg.V_next(i) = brt_value(brt{k+1}, gv, guam_trace.lon(:, i) - traj.trim_lon(:, k+1));
end

gapg = evg.V_cur >= 0 & evg.V_next >= 0;
fprintf('\n=== CBF detection on the flown trajectory ===\n');
fprintf('flown samples outside BOTH scheduled BRTs: %d / %d\n', nnz(gapg), Mg);
fprintf(' k   handover t   V_next@handover   ready?   actual entry t   delay [s]   gap smpl\n');
n_ready = 0;
for k = 1:traj.n_trim - 1
    t_h = traj.t_node(k + 1);
    [~, ih] = min(abs(tg - t_h));
    % evaluate against BRT_{k+1} explicitly (evg.V_next at exactly t_h
    % already belongs to the NEXT scheduled segment - off-by-one trap)
    Vh = brt_value(brt{k+1}, gv, guam_trace.lon(:, ih) - traj.trim_lon(:, k+1));
    ready = Vh < 0;
    n_ready = n_ready + ready;
    % first time the flown state enters BRT_{k+1} (search from segment
    % start onward, past the scheduled handover if needed)
    m0 = find(tg >= traj.t_node(k), 1, 'first');
    t_enter = NaN;
    for ii = m0:Mg
        Vn = brt_value(brt{k+1}, gv, guam_trace.lon(:, ii) - traj.trim_lon(:, k+1));
        if Vn < 0, t_enter = tg(ii); break; end
    end
    ngap = nnz(gapg & evg.seg == k);
    fprintf('%2d   %7.2f      %+8.3f        %s     %8.2f      %+7.2f      %5d\n', ...
            k, t_h, Vh, ternary(ready, ' yes ', ' NO  '), t_enter, ...
            t_enter - t_h, ngap);
end
fprintf('segments transition-ready at scheduled handover: %d / %d\n', ...
        n_ready, traj.n_trim - 1);

%% 10) Per-segment symptom table + detection figure
% Restates the same data in the vocabulary used when reading the tube
% figures: did the vehicle start outside BRT_k, did it ever reach
% BRT_{k+1}, did it fall back out before the handover.
fprintf('\n k  | start in BRT_k | reached BRT_k+1 | re-exit before handover | ready@handover\n');
for k = 1:traj.n_trim - 1
    m  = find(tg >= traj.t_node(k) & tg <= traj.t_node(k + 1));
    Vc = zeros(1, numel(m));  Vn = zeros(1, numel(m));
    for j = 1:numel(m)
        Vc(j) = brt_value(brt{k},     gv, guam_trace.lon(:, m(j)) - traj.trim_lon(:, k));
        Vn(j) = brt_value(brt{k + 1}, gv, guam_trace.lon(:, m(j)) - traj.trim_lon(:, k + 1));
    end
    ien = find(Vn < 0, 1, 'first');
    reached = ~isempty(ien);
    reexit  = reached && any(Vn(ien:end) >= 0);
    fprintf('%2d  |      %s       |       %s        |          %s            |      %s\n', ...
        k, ternary(Vc(1) < 0, 'yes', 'NO '), ternary(reached, 'yes', 'NO '), ...
        ternary(reexit, 'YES', 'no '), ternary(Vn(end) < 0, 'yes', 'NO '));
end

f9 = plot_cbf_detection(traj, brt, gv, guam_trace, struct('mark', [1 2 5 6]));

%% 11) PREDICTIVE gate validation (agreed CBF-style forecast)
% Gate: V_hat = V_next + dV/dt * (t_rem - t_mar) must be < 0, with dV/dt a
% trailing finite difference of the measured V (LivenessPredictor). This
% section replays the flown baseline and asks: how much EARLIER does the
% forecast warn, compared to when the problem actually materializes?
fd_window = 0.2;    % [s] finite-difference window
t_margin  = 0.3;    % [s] arrive-early margin t_mar
fprintf('\n=== predictive gate (fd=%.2fs, t_mar=%.2fs) vs actual events ===\n', ...
        fd_window, t_margin);
fprintf(' k | first warn t | warn mode   | actual event t | event        | lead [s]\n');
for k = 1:traj.n_trim - 1
    t0 = traj.t_node(k);  t1 = traj.t_node(k + 1);
    m  = find(tg >= t0 & tg <= t1);
    P  = LivenessPredictor(fd_window, t_margin);
    t_warn = NaN;  warn_mode = '     -     ';
    Vn = zeros(1, numel(m));
    for j = 1:numel(m)
        Vn(j) = brt_value(brt{k+1}, gv, guam_trace.lon(:, m(j)) - traj.trim_lon(:, k+1));
        % 3rd argument is the REMAINING time to the node, not the absolute
        % node time (passing t1 inflates the horizon ~20x and floods the
        % late segments with spurious re-exit warnings).
        pr = P.update(tg(m(j)), Vn(j), t1 - tg(m(j)));
        if isnan(t_warn) && ~pr.ready_pred && tg(m(j)) > t0 + fd_window
            t_warn = tg(m(j));
            warn_mode = ternary(Vn(j) >= 0, 'late-entry ', 're-exit    ');
        end
    end
    % actual event: not inside at handover (late entry never fixed) or
    % re-exit moment (V back above 0 after having been inside)
    ien = find(Vn < 0, 1, 'first');
    t_evt = NaN;  evt = 'none        ';
    if isempty(ien) || Vn(end) >= 0
        if isempty(ien)
            t_evt = t1;  evt = 'never enters';
        else
            ire = find(Vn(ien:end) >= 0, 1, 'first');
            t_evt = tg(m(ien + ire - 1));  evt = 're-exit     ';
        end
    end
    fprintf('%2d |    %7.2f   | %s |    %7.2f     | %s |  %6.2f\n', ...
            k, t_warn, warn_mode, t_evt, evt, t_evt - t_warn);
end

% -------------------------------------------------------------------------
function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end

function V = brt_value(data, gv, x_pert)
% Interpolate one BRT value function at a perturbation state (clamped to
% the grid box, matching ValueFunctionLUT's clamping).
xc = x_pert(:)';
for d = 1:4
    xc(d) = min(max(xc(d), gv{d}(1)), gv{d}(end));
end
V = interpn(gv{1}, gv{2}, gv{3}, gv{4}, data, xc(1), xc(2), xc(3), xc(4), 'linear');
end
