% compare_cbf_on_off - fly the searched schedule with the liveness filter off
% and on, and measure what the filter costs.
%
% ---------------------------------------------------------------------------
% WHAT THIS IS AND IS NOT
% ---------------------------------------------------------------------------
% This is a compatibility check, not a joint design. The schedule comes from
% optimize_schedule_sim, which searched with the filter OFF: it is the fastest
% ramp the aircraft can follow while staying inside a reachable set on its own.
% The question here is only whether adding the filter behind it changes
% anything, and by how much.
%
% If the filter never fires, the two architectures agree - the schedule already
% keeps the aircraft where the filter would have kept it. If it fires, the
% aircraft stops tracking the ramp and the transition stretches; how much it
% stretches is the number this script reports.
%
% ---------------------------------------------------------------------------
% THE FILTER WAS NOT RUNNING BEFORE THIS
% ---------------------------------------------------------------------------
% FilterConfig pointed at 'tables/BRT', which does not exist here, with the
% channel prefix 'GUAM_LON_HJIR' while the files are named
% 'GUAM_LON_BRT_HJIR_UH*_WH*.mat' in data/. The lookup matched zero files and
% reported available = 0, and nothing downstream checks that, so every run
% with mode 'blend' silently filtered nothing. Both were corrected; the filter
% now loads all 20 longitudinal tables. Any earlier result that assumed an
% active filter needs re-examining.
%
% ---------------------------------------------------------------------------
% WHAT IS MEASURED
% ---------------------------------------------------------------------------
% The reference is identical in both runs, so "reference reaches cruise" is
% identical too and says nothing. What changes is when the AIRCRAFT gets
% there, so that is the number compared, along with how often the filter was
% active and whether either run left the sets.
%
% Output: logger/cbf_compare.mat, logger/cbf_compare.png
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt       = 0.01;
WDOT     = 11.667;
T_SETTLE = 45.0;
T_CRUISE = 25.0;
ALT0     = 1600;
LON      = [1 3 5 11];

S  = load(fullfile(root,'logger','schedule_sim.mat'));
a  = S.a;   UH = S.UH;   du_seg = S.du_seg;
T   = load('trim_table_Poly_ConcatVer4p0.mat');
trim_lon = T.XU0_interp(LON, 1:20, 3);

spec = FilterConfig.channelSpec('lon');
gv = cell(1,4);
for d = 1:4
    gv{d} = linspace(spec.grid_min(d), spec.grid_max(d), spec.grid_num(d));
end
V = cell(1,20);
for k = 1:20
    Sd = load(fullfile(root,'data', ...
        sprintf('GUAM_LON_BRT_HJIR_UH%d_WH%d.mat',k,3)), 'data');
    V{k} = griddedInterpolant(gv, Sd.data, 'linear', 'nearest');
end

line = @(c) fprintf('%s\n', repmat(c, 1, 88));
fprintf('schedule from optimize_schedule_sim, reference time %.2f s\n\n', sum(du_seg./a));

R = struct();
modes = {'off', 'blend'};
for im = 1:2
    md = modes{im};
    params = struct('filter_mode', md);  params.T_seg = 2.0;
    guam = LpC_GUAM(Config('trim_schedule', params));
    guam.reset();
    st = guam.saveState();  st.state(3) = -ALT0;  guam.restoreState(st);
    p = [0; 0; -ALT0];

    for i = 1:round(T_SETTLE/dt)
        p = p + [0; 0; WDOT*dt];
        guam.step(struct('pos',p,'vel',[0;0;WDOT],'chi',0,'chi_dot',0));
    end
    f = guam.controller.safety_filter;
    fprintf('  mode %-6s : filter tables loaded = %d\n', md, f.lut.available);
    f.reset_counters();

    N_tr = sum(max(round(du_seg./a/dt), 1));
    N_cr = round(T_CRUISE/dt);
    N = N_tr + N_cr;
    t = zeros(1,N);  uc = zeros(1,N);  ua = zeros(1,N);
    Vb = nan(1,N);   seg = zeros(1,N);
    tt = 0;  n = 0;

    for k = 1:19
        Nk = max(round(du_seg(k)/a(k)/dt), 1);
        for i = 1:Nk
            n = n + 1;
            s = guam.state;  th = s(8);
            uh =  s(4)*cos(th) + s(6)*sin(th);
            wh = -s(4)*sin(th) + s(6)*cos(th);
            xa = [uh; wh; s(11); th];
            Vb(n) = set_value(V, trim_lon, spec, xa, k, k+1);
            v = UH(k) + a(k)*(i*dt);
            t(n) = tt;  uc(n) = v;  ua(n) = uh;  tt = tt + dt;  seg(n) = k;
            p = p + [v*dt; 0; WDOT*dt];
            guam.step(struct('pos',p,'vel',[v;0;WDOT],'chi',0,'chi_dot',0));
        end
    end
    for i = 1:N_cr
        n = n + 1;
        s = guam.state;  th = s(8);
        uh =  s(4)*cos(th) + s(6)*sin(th);
        wh = -s(4)*sin(th) + s(6)*cos(th);
        Vb(n) = set_value(V, trim_lon, spec, [uh; wh; s(11); th], 20, 20);
        t(n) = tt;  uc(n) = UH(end);  ua(n) = uh;  tt = tt + dt;  seg(n) = 20;
        p = p + [UH(end)*dt; 0; WDOT*dt];
        guam.step(struct('pos',p,'vel',[UH(end);0;WDOT],'chi',0,'chi_dot',0));
    end

    R(im).mode = md;  R(im).t = t(1:n);  R(im).uc = uc(1:n);  R(im).ua = ua(1:n);
    R(im).Vb = Vb(1:n);  R(im).seg = seg(1:n);
    R(im).n_calls = f.n_calls;  R(im).n_active = f.n_active;
end

%% report --------------------------------------------------------------
line('='); fprintf('LIVENESS FILTER OFF vs ON\n'); line('=');
Tref = sum(du_seg ./ a);
fprintf('  reference reaches cruise at %.2f s in both runs (same command)\n\n', Tref);
fprintf('  %-8s | filter active | aircraft 95%% | 99%% | 100%% | worst V | left sets\n', 'mode');
for im = 1:2
    r = R(im);
    i95  = first_at(r.ua, 0.95*UH(end));
    i99  = first_at(r.ua, 0.99*UH(end));
    i100 = first_at(r.ua, UH(end));
    if r.n_calls > 0, act = sprintf('%5.1f%%', 100*r.n_active/r.n_calls);
    else,             act = ' n/a  '; end
    fprintf('  %-8s | %13s | %11s | %5s | %5s | %7.4f | %s\n', r.mode, act, ...
        tstr(r,i95), tstr(r,i99), tstr(r,i100), max(r.Vb), ...
        ternary(any(r.Vb > 0), 'YES', 'no'));
end

i100off = first_at(R(1).ua, UH(end));
i100on  = first_at(R(2).ua, UH(end));
if ~isnan(i100off) && ~isnan(i100on)
    fprintf('\n  the filter costs %+.2f s to reach cruise\n', R(2).t(i100on) - R(1).t(i100off));
end
fprintf('  filter fired on %.1f%% of steps (%d of %d)\n', ...
        100*R(2).n_active/max(R(2).n_calls,1), R(2).n_active, R(2).n_calls);

save(fullfile(root,'logger','cbf_compare.mat'), 'R', 'a', 'UH', 'du_seg', 'Tref');

f1 = figure('Position',[80 80 1100 400], 'Color','w');
subplot(1,2,1); hold on; grid on;
plot(R(1).t, R(1).ua, 'LineWidth',1.6);
plot(R(2).t, R(2).ua, '--', 'LineWidth',1.6);
plot(R(1).t, R(1).uc, ':', 'LineWidth',1.2);
xlabel('t [s]'); ylabel('u [ft/s]'); title('flown speed');
legend({'filter off','filter on','command'}, 'Location','southeast');
subplot(1,2,2); hold on; grid on;
plot(R(1).t, R(1).Vb, 'LineWidth',1.6);
plot(R(2).t, R(2).Vb, '--', 'LineWidth',1.6);
yline(0,'r--'); xlabel('t [s]'); ylabel('V'); title('set value (V<0 = inside)');
legend({'filter off','filter on'}, 'Location','best');
saveas(f1, fullfile(root,'logger','cbf_compare.png'));
fprintf('\nsaved logger/cbf_compare.mat and logger/cbf_compare.png\n');

%% ---------------------------------------------------------------------
function v = set_value(V, TL, spec, xa, k1, k2)
% Best of the two bracketing anchors, each read at the deviation from its own
% trim point. Outside the stored grid is not certifiable, so it reads as
% unsafe rather than being clamped to a boundary value.
v = inf;
for k = unique([k1 k2])
    e = xa - TL(:,k);
    if all(e >= spec.grid_min & e <= spec.grid_max)
        v = min(v, V{k}(e(1), e(2), e(3), e(4)));
    end
end
if ~isfinite(v), v = 1; end
end

function i = first_at(u, target)
i = find(u >= target, 1, 'first');
if isempty(i), i = NaN; end
end

function s = tstr(r, i)
if isnan(i), s = '  --  '; else, s = sprintf('%6.2f', r.t(i)); end
end

function s = ternary(c,a,b), if c, s=a; else, s=b; end, end
