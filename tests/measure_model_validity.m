% measure_model_validity - how large may the tracking error get before the
% linear closed-loop model stops describing the aircraft?
%
% ---------------------------------------------------------------------------
% WHY THIS IS THE BLOCKING QUESTION
% ---------------------------------------------------------------------------
% The schedule search is now limited by a hand-picked cap on how far the linear
% model is trusted, not by the reachable sets. With no cap it returns 23.6 s but
% propagates Phi to a 43.9 ft/s error - Phi was measured with 0.1 ft/s nudges,
% so that is extrapolation by a factor of 400. With a 4 ft/s cap it returns
% 133.8 s and every segment sits exactly on the cap. Neither number is a
% property of the vehicle; both are properties of a guess.
%
% So the guess is replaced with a measurement: fly the real nonlinear simulator
% and the linear model side by side on the same command, and watch where they
% part company.
%
% ---------------------------------------------------------------------------
% THE EXPERIMENT
% ---------------------------------------------------------------------------
% Uniform acceleration, several values, each flown twice:
%
%   real   : LpC_GUAM stepped with a descending reference, full nonlinear aero
%   linear : d(i+1) = Phi d(i) - xi_e'(v) a dt, from the descending model
%
% Both are driven by the identical command, both start trimmed at hover, and
% both are read out as the heading-frame [u w q theta] the reachable sets use.
% What is reported is the disagreement between them AS A FUNCTION OF how large
% the error has grown - which is the curve the cap should come from.
%
% The mission descends at WH3 (11.667 ft/s), the condition the trim points and
% their sets were built at; flown level the aircraft sits a standing 11.7 ft/s
% off every anchor and none of this means anything.
%
% Output: logger/model_validity.mat, logger/model_validity.png
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

%% 0. setup ------------------------------------------------------------
dt     = 0.01;
WDOT   = 11.667;          % WH3 descent [ft/s, + down]
A_LIST = [2 3 5 8];       % uniform accelerations to probe [ft/s^2]
ALT0   = 1800;            % start altitude [ft]; the longest case costs ~1500
% The hover trim point itself descends at 11.667 ft/s, so "trimmed at hover"
% is not "at rest". Starting the simulator at rest and holding 6 s left it
% still 11.7 ft/s short of the trim condition when the ramp began, and that
% showed up as a spurious model error of exactly 11.67. The loop needs tens of
% seconds to acquire the descent - measured 32 s to reach 11.4 - so the hold is
% long, and the statistics below start only once the ramp does.
T_HOLD = 45.0;            % settle at the descending hover trim [s]
SEL    = [4 6 11 8];
LON    = [1 3 5 11];

T   = load('trim_table_Poly_ConcatVer4p0.mat');
UH  = T.UH(1:20);  UH = UH(:)';
XU0 = T.XU0_interp(:, 1:20, 3);
trim_lon = XU0(LON, :);

Sm  = load(fullfile(root,'logger','closed_loop_model_descending.mat'));
M   = Sm.M;   dXE = Sm.dXE;
params = struct('filter_mode','off');  params.T_seg = 2.0;

line = @(c) fprintf('%s\n', repmat(c, 1, 92));
b2h  = @(x) [ x(1,:).*cos(x(4,:)) + x(2,:).*sin(x(4,:)) ;
             -x(1,:).*sin(x(4,:)) + x(2,:).*cos(x(4,:)) ;
              x(3,:) ; x(4,:) ];

R = struct();
for ia = 1:numel(A_LIST)
    a  = A_LIST(ia);
    Tr = UH(end)/a;
    N  = round((T_HOLD + Tr)/dt);

    %% command profile: hold hover, then ramp at a, descending throughout
    t  = (0:N-1)*dt;
    v  = min(max(a*(t - T_HOLD), 0), UH(end));
    z  = -ALT0 + WDOT*t;
    x  = cumtrapz(t, v);

    %% --- real simulator ---------------------------------------------
    guam = LpC_GUAM(Config('trim_schedule', params));  guam.reset();
    s = guam.saveState();  s.state(3) = -ALT0;  guam.restoreState(s);
    XR = zeros(4, N);
    for i = 1:N
        st = guam.state;
        XR(:, i) = [st(4); st(6); st(11); st(8)];   % u, w (body), q, theta
        guam.step(struct('pos',[x(i);0;z(i)], 'vel',[v(i);0;WDOT], ...
                         'chi',0, 'chi_dot',0));
    end

    %% --- linear model on the identical command ----------------------
    XL = zeros(4, N);
    d  = zeros(34, 1);
    for i = 1:N
        vv = v(i);
        [P, xe] = interp_model(M, dXE, vv);
        XL(:, i) = xe(SEL) + d(SEL);
        adot = 0;
        if i < N, adot = (v(i+1) - v(i))/dt; end
        d = P*d - dXE(:, anchor_of(M.u_anchor, vv))*adot*dt;
    end

    %% --- compare in the heading frame the sets use ------------------
    HR = b2h(XR);   HL = b2h(XL);
    % error relative to the interpolated trim point, which is what a lookup sees
    ER = HR - interp_trim(trim_lon, UH, v);
    EL = HL - interp_trim(trim_lon, UH, v);
    DIS = abs(EL - ER);                       % model minus truth

    R(ia).a = a;  R(ia).t = t;  R(ia).v = v;
    R(ia).ER = ER;  R(ia).EL = EL;  R(ia).DIS = DIS;
    R(ia).diverged = any(~isfinite(XR(:))) || max(abs(XR(1,:))) > 400;
end

%% report --------------------------------------------------------------
line('='); fprintf('LINEAR MODEL vs SIMULATOR, uniform acceleration\n'); line('=');
fprintf('  "true error" is how far the simulator strays from the trim point.\n');
fprintf('  "model error" is how far the linear prediction is from the simulator.\n');
fprintf('  The cap belongs where the second becomes a large fraction of the first.\n\n');
fprintf('    a   |  true error  (u, w, q[deg/s], th[deg])  |  model error (same)      | ratio\n');
for ia = 1:numel(R)
    if R(ia).diverged
        fprintf(' %5.1f |  *** simulator diverged ***\n', R(ia).a);  continue;
    end
    m  = R(ia).t >= R(ia).t(1) + 45.0;        % ramp only, settle excluded
    te = max(abs(R(ia).ER(:,m)), [], 2);
    me = max(R(ia).DIS(:,m), [], 2);
    fprintf(' %5.1f | %6.2f %6.2f %7.2f %7.2f | %6.2f %6.2f %7.2f %7.2f | %5.0f%%\n', ...
        R(ia).a, te(1), te(2), rad2deg(te(3)), rad2deg(te(4)), ...
        me(1), me(2), rad2deg(me(3)), rad2deg(me(4)), ...
        100*max(me(1:2)./max(te(1:2), 1e-6)));
end

save(fullfile(root,'logger','model_validity.mat'), 'R', 'A_LIST', 'WDOT', 'ALT0');

f = figure('Position',[80 80 1100 700], 'Color','w');
nm = {'u [ft/s]','w [ft/s]','q [deg/s]','\theta [deg]'};
sc = [1 1 180/pi 180/pi];
for ch = 1:4
    subplot(2,2,ch); hold on; grid on;
    for ia = 1:numel(R)
        if R(ia).diverged, continue; end
        plot(R(ia).t, sc(ch)*R(ia).ER(ch,:), 'LineWidth', 1.5);
        plot(R(ia).t, sc(ch)*R(ia).EL(ch,:), '--', 'LineWidth', 1.2);
    end
    xlabel('t [s]'); ylabel(nm{ch});
    title([nm{ch} ' : solid = simulator, dashed = linear model']);
end
saveas(f, fullfile(root,'logger','model_validity.png'));
fprintf('\nsaved logger/model_validity.mat and logger/model_validity.png\n');

%% helpers -------------------------------------------------------------
function k = anchor_of(ua, v)
k = find(ua <= v, 1, 'last');
if isempty(k), k = 1; end
k = min(k, numel(ua));
end

function [P, xe] = interp_model(M, dXE, v) %#ok<INUSD>
ua = M.u_anchor;
if v <= ua(1),        k = 1;              al = 0;
elseif v >= ua(end),  k = numel(ua)-1;    al = 1;
else
    k = find(ua <= v, 1, 'last');  k = min(k, numel(ua)-1);
    al = (v - ua(k)) / (ua(k+1) - ua(k));
end
P  = (1-al)*M.Phi{k}    + al*M.Phi{k+1};
xe = (1-al)*M.xi_e(:,k) + al*M.xi_e(:,k+1);
end

function XT = interp_trim(trim_lon, UH, v)
XT = zeros(4, numel(v));
for i = 1:numel(v)
    XT(:, i) = interp1(UH, trim_lon', min(max(v(i), UH(1)), UH(end)))';
end
end
