% diag_reexit_ungoverned - why does the UNGOVERNED run enter BRT(k+1) during
% segment k and then leave it again near the end?
%
% Runs the plain reference trajectory through RSLQR (no governor) and, for one
% segment, reports
%   - V of BRT(k+1) at the flown state, so entry/exit instants are visible
%   - the speed the autopilot is really chasing: ref_vel + 0.1*e_pos_body
%     (RSLQR.m:227), against the plan
%   - dV/dt split per state, so the motion responsible for the exit is named
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

SEG = 5;                      % segment k: trim k -> k+1

dt = 0.01;
params = struct('filter_mode', 'off');  params.T_seg = 2.0;
traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt, params));
brtV = BRTValue(fullfile(root, 'data'), traj.trim_lon, traj.wh_idx);
guam = LpC_GUAM(Config('trim_schedule', params));

%% Ungoverned run: the schedule fed straight to RSLQR
rt = guam.refTraj;  N = size(rt.pos, 2);
guam.reset();
L = struct('t', rt.time, 'st', zeros(12, N), 'r', zeros(1, N), ...
           'e_pos', zeros(1, N), 'u_eff', zeros(1, N));
for i = 1:N
    s = guam.state;
    ref = struct('pos', rt.pos(:, i), 'vel', rt.vel(:, i), ...
                 'chi', rt.chi(i), 'chi_dot', rt.chidot(i));
    R = RSLQR.rotm_i2b(s(7), s(8), s(9));
    e = R * (ref.pos - s(1:3));
    L.st(:, i) = s;
    L.r(i)     = ref.vel(1);
    L.e_pos(i) = e(1);
    L.u_eff(i) = ref.vel(1) + 0.1 * e(1);      % what the inner loop chases
    guam.step(ref);
end

%% V of BRT(k+1) along the flown trajectory
% Calling value() with v = u_anchor(k+1) gives weight 1 to brt{k+1} at its own
% trim, i.e. exactly the "next BRT" test.
v_next = brtV.u_anchor(SEG + 1);
V = zeros(1, N);
for i = 1:N, V(i) = brtV.value(L.st([4 6 11 8], i), v_next); end

ta = traj.t_node(SEG);  tb = traj.t_node(SEG + 1);
w  = L.t >= ta - 1.5 & L.t <= tb + 3.0;
iw = find(w);

inside = V < 0;
i_ent = iw(find(inside(iw), 1, 'first'));
i_ext = [];
if ~isempty(i_ent)
    j = find(~inside(i_ent:iw(end)), 1, 'first');
    if ~isempty(j), i_ext = i_ent + j - 1; end
end

fprintf('\n=== ungoverned, segment %d (trim %d -> %d, %.1f -> %.1f ft/s) ===\n', ...
        SEG, SEG, SEG+1, brtV.u_anchor(SEG), v_next);
fprintf('plan window : %.2f - %.2f s\n', ta, tb);
if isempty(i_ent)
    fprintf('never entered BRT(%d) in this window\n', SEG+1);  return;
end
fprintf('entered BRT(%d) at %.2f s   (V = %+.3f)\n', SEG+1, L.t(i_ent), V(i_ent));
if isempty(i_ext)
    fprintf('did not leave again in this window\n');
else
    fprintf('LEFT     BRT(%d) at %.2f s   (V = %+.3f)   after %.2f s inside\n', ...
            SEG+1, L.t(i_ext), V(i_ext), L.t(i_ext) - L.t(i_ent));
end

%% What was the autopilot actually chasing?
fprintf('\n   t   |  plan r | +0.1*e_pos | eff cmd |  flown u | u - eff\n');
for tt = ta : 0.4 : min(tb + 2.4, L.t(end))
    [~, i] = min(abs(L.t - tt));
    fprintf('%6.2f | %7.1f | %10.2f | %7.1f | %8.1f | %+7.2f\n', ...
            L.t(i), L.r(i), 0.1*L.e_pos(i), L.u_eff(i), L.st(4,i), ...
            L.st(4,i) - L.u_eff(i));
end

%% dV/dt split per state at the exit
if ~isempty(i_ext)
    i = i_ext;
    x = L.st([4 6 11 8], i);
    h = [0.05; 0.05; 0.002; 0.002];
    g = zeros(4,1);
    for d = 1:4
        xp = x; xp(d) = xp(d) + h(d);
        xm = x; xm(d) = xm(d) - h(d);
        g(d) = (brtV.value(xp, v_next) - brtV.value(xm, v_next)) / (2*h(d));
    end
    xd = (L.st([4 6 11 8], i+1) - L.st([4 6 11 8], i-1)) / (2*dt);
    contrib = g .* xd;
    nm = {'u    ', 'w    ', 'q    ', 'theta'};
    fprintf('\ndV/dt split at the exit (t = %.2f s), total %+.3f /s\n', ...
            L.t(i), sum(contrib));
    for d = 1:4
        fprintf('  %s  dV/dx = %+8.4f   dx/dt = %+8.4f   -> %+7.3f  (%5.1f %%)\n', ...
                nm{d}, g(d), xd(d), contrib(d), ...
                100*contrib(d)/sum(abs(contrib)));
    end
    fprintf('\nstate at exit, relative to trim %d:\n', SEG+1);
    dxe = x - traj.trim_lon(:, SEG+1);
    fprintf('  du = %+7.2f ft/s   dw = %+7.2f ft/s   dq = %+7.3f rad/s   dtheta = %+7.2f deg\n', ...
            dxe(1), dxe(2), dxe(3), rad2deg(dxe(4)));
end
