% run_option1_mission_fig - whole-mission 3-D path for the dV/dt governor.
%
% Runs option 1 (TrimProgressGovernor with the LivenessPredictor gate, i.e.
% the prediction made by extrapolating the BRT value with its time derivative,
% V_hat = V_next + dV/dt*(t_rem - t_mar)) and draws the same (u, t, theta)
% figure used for the reference governor:
%
%   black dashed  scheduled trajectory  x_e(r)
%   blue          governed command      x_e(v)
%   red           flown state           (u, theta)
%   faint tubes   BRT of each trim point, at the instant the governed command
%                 first reached that trim speed (so each is one stored table,
%                 not a blend of two)
%
% Saves the figure as .fig (MATLAB's own format, re-openable and rotatable)
% next to a .png copy.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

%% Setup
dt_sim = 0.01;
params = struct('filter_mode', 'off');
params.T_seg = 2.0;

traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt_sim, params));
brtV = BRTValue(fullfile(root, 'data'), traj.trim_lon, traj.wh_idx);

% TrimProgressGovernor wants the raw tables and grid vectors
spec = FilterConfig.channelSpec('lon');
gv = cell(1, 4);
for d = 1:4
    gv{d} = linspace(spec.grid_min(d), spec.grid_max(d), spec.grid_num(d));
end
brt = cell(1, traj.n_trim);
for k = 1:traj.n_trim
    S = load(fullfile(root, 'data', ...
        sprintf('GUAM_LON_BRT_HJIR_UH%d_WH%d.mat', k, traj.wh_idx)), 'data');
    brt{k} = S.data;
end

guam = LpC_GUAM(Config('trim_schedule', params));
gov  = TrimProgressGovernor(traj, brt, gv, dt_sim, 0.05, ...
                            struct('use_predictor', true));
guam.reset();  gov.reset();

%% Governed run
T_max = 80;  Ng = round(T_max/dt_sim);
G = struct('t', zeros(1,Ng), 'st', zeros(12,Ng), 'v', zeros(1,Ng), ...
           'progress', zeros(1,Ng));
for i = 1:Ng
    tt = (i-1)*dt_sim;
    [ref, gi] = gov.step(guam.state, tt);
    G.t(i) = tt;  G.st(:,i) = guam.state;
    G.v(i) = gi.u_ref;  G.progress(i) = gi.progress;
    guam.step(ref);
    if gi.progress >= traj.n_trim && tt > traj.t_node(end) + 3
        fn = fieldnames(G);
        for f = 1:numel(fn), G.(fn{f}) = G.(fn{f})(:,1:i); end
        break;
    end
end
n = numel(G.t);
G.r = interp1(traj.time, traj.lon(1,:), min(G.t, traj.time(end)));
i_fin = find(G.progress >= traj.n_trim, 1, 'first');
t_done = G.t(i_fin);

fprintf('\n=== option 1 (dV/dt predictive gate), T_seg = %.2f s ===\n', params.T_seg);
fprintf('final trim reached : %.2f s   (plan wanted %.2f s, delay %+.2f s)\n', ...
        t_done, traj.t_node(end), t_done - traj.t_node(end));
fprintf('frozen             : %.2f s\n', gov.n_hold * dt_sim);
Vnow = zeros(1,n);
for i = 1:n, Vnow(i) = brtV.value(G.st([4 6 11 8],i), G.v(i)); end
fprintf('samples with V >= 0: %d / %d (%.1f %%)\n', ...
        nnz(Vnow >= 0), n, 100*nnz(Vnow >= 0)/n);
fprintf('worst V            : %+.3f\n', max(Vnow));
fprintf('altitude deviation : %.1f ft\n', max(abs(-G.st(3,:) - 80)));

%% Whole-mission 3-D figure
f = figure('Name', 'option 1 (dV/dt): full mission', 'Position', [70 50 1150 800]);
ax = axes(f); hold(ax, 'on'); grid(ax, 'on');

% one tube per trim point, at the instant the command first reached it
u_a = brtV.u_anchor;
t_sl = nan(1,numel(u_a));  i_sl = nan(1,numel(u_a));
for k = 1:numel(u_a)
    ik = find(G.v >= u_a(k) - 1e-9, 1, 'first');
    if ~isempty(ik), t_sl(k) = G.t(ik);  i_sl(k) = ik; end
end
keep = ~isnan(t_sl);  t_sl = t_sl(keep);  i_sl = i_sl(keep);

cmap = parula(256);
for k = 1:numel(t_sl)
    col = cmap(round((k-1)/max(numel(t_sl)-1,1)*255)+1, :);
    C = brt_zero_contour(brtV, G.v(i_sl(k)), G.st([4 6 11 8], i_sl(k)));
    for c = 1:numel(C)
        p = C{c};
        plot3(ax, p(1,:), t_sl(k)*ones(1,size(p,2)), p(2,:), ...
              'Color', [col 0.45], 'LineWidth', 1.0, 'HandleVisibility', 'off');
    end
end
plot3(ax, nan, nan, nan, '-', 'Color', [0.35 0.55 0.55], 'LineWidth', 1.0, ...
      'DisplayName', sprintf('BRT of each trim point (%d of %d reached)', ...
                             numel(t_sl), numel(u_a)));

xr = zeros(2,n);  xv = zeros(2,n);
for i = 1:n
    er = brtV.trim_at(G.r(i));  xr(:,i) = [er(1); rad2deg(er(4))];
    ev = brtV.trim_at(G.v(i));  xv(:,i) = [ev(1); rad2deg(ev(4))];
end
plot3(ax, xr(1,:), G.t, xr(2,:), '--', 'Color', [0.15 0.15 0.15], ...
      'LineWidth', 1.8, 'DisplayName', 'scheduled trajectory  x_e(r)');
plot3(ax, xv(1,:), G.t, xv(2,:), '-', 'Color', [0 0.2 0.85], ...
      'LineWidth', 2.4, 'DisplayName', 'governed command  x_e(v)');
plot3(ax, G.st(4,:), G.t, rad2deg(G.st(8,:)), '-', 'Color', [0.85 0.2 0.2], ...
      'LineWidth', 2.0, 'DisplayName', 'flown state  (u, \theta)');

xlabel(ax, 'u [ft/s]'); ylabel(ax, 'time [s]'); zlabel(ax, '\theta [deg]');
title(ax, sprintf(['Option 1, dV/dt predictive gate: plan reached the final ' ...
      'trim at %.2f s, the governor at %.2f s'], traj.t_node(end), t_done));
legend(ax, 'Location', 'northeast');  view(ax, -60, 22);
ylim(ax, [G.t(1) G.t(end)]);

%% Save
out = fullfile(root, 'logger');
savefig(f, fullfile(out, 'option1_dvdt_mission.fig'));
exportgraphics(f, fullfile(out, 'option1_dvdt_mission.png'), 'Resolution', 150);
fprintf('\nsaved %s(.fig/.png)\n', fullfile(out, 'option1_dvdt_mission'));
