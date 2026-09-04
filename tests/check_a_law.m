% check_a_law - fly the schedule produced by the closed-form acceleration law.
%
% ---------------------------------------------------------------------------
% THE LAW
% ---------------------------------------------------------------------------
% The proposal is that the admissible acceleration at each anchor follows from
% geometry alone - no prediction, no iteration, no simulation loop:
%
%     a_k = eta * min(  (P_th - Gam*Dth) / (2*mu_th) ,  min(rho_th)/mu_th )
%     P_th = rho_th(k) + rho_th(k+1)          tube pitch room, both ends
%     Gam  = 1/sqrt(1 - (Du/P_u)^2)           supposed speed-axis coupling
%     mu_th(k) = lambda(k)/g                  pitch tilt needed per unit accel
%
% then T_k = Du / a_k. The physical claim behind mu_th = lambda/g is that
% accelerating forward at a requires tilting by a/g, times whatever fraction
% lambda of the acceleration the tilt has to supply - the pusher supplies the
% rest, so lambda falls from 1.0 at hover to 0.475 by u = 127.
%
% ---------------------------------------------------------------------------
% WHAT THIS SCRIPT FINDS
% ---------------------------------------------------------------------------
%  1. At eta = 1 the law returns 15.5 s, far outside anything that has flown.
%     Calibrating eta on the low-speed band alone (2.5 s/segment, already
%     certified) gives eta = 0.458 and a total of 33.9 s - within 1 % of the
%     33.7 s the flight search found. So the TOTAL is right.
%
%  2. The DISTRIBUTION is not. At equal total time the law's schedule leaves
%     the tube (6.5 %, worst V = +0.34 at u = 141) where the band schedule
%     does not. Giving the law 40 s instead of 34 does not repair it.
%
%  3. Swapping halves with the band schedule isolates the fault: segments
%     1-12 (u = 0-101) from the law are SAFE, segments 13-19 are not. The law
%     is right over the first two thirds and wrong over the last third.
%
%  4. Repairing the suspected lambda artifact at anchors 17-20 makes it far
%     worse (V = +2.0), which identifies the cause: at high speed the binding
%     limit is not pitch attitude at all - it is thrust minus drag, and the
%     law has no term for it. The artifact was accidentally masking that.
%
% Output: logger/a_law_check.mat / .png
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt     = 0.01;
WH_IDX = 3;
N_TRIM = 20;
T_HOLD = 5.0;
G      = 32.174;
DTH    = 0.01047;                      % anchor-to-anchor trim pitch step [rad]
AX     = {'u (fwd speed)', 'w (vert speed)', 'q (pitch rate)', 'theta (pitch)'};

T = load(fullfile(root,'logger','tube_dimensions.mat'));   % half : 4 x 20
O = load(fullfile(root,'logger','segment_time_opt.mat'));
half = T.half;  u_anchor = O.u_anchor(1:N_TRIM);

params = struct('filter_mode','off');  params.T_seg = 2.0;
traj0 = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt, params));
brtV  = BRTValue(fullfile(root,'data'), traj0.trim_lon, WH_IDX);

% lambda(u) from the closed-loop DC gain, as reported. NOT verified here; the
% value at anchor 17 jumping 0.475 -> 1.045 is a known suspect (table splice).
uh_l  = [1 5 8 10 13 16 17 20];
lam   = interp1(uh_l, [1.000 1.004 0.972 0.849 0.614 0.475 1.045 1.069], 1:N_TRIM);
lam_r = lam;  lam_r(17:20) = 0.475;    % repair: lambda cannot rise again at high speed

%% ---------------------------------------------------------------------
%  1. the law at eta = 1, and eta calibrated on the low band alone
%% ---------------------------------------------------------------------
[Traw, tab] = law_schedule(half, u_anchor, lam/G, DTH, N_TRIM);
fprintf('--- the law at eta = 1 ---\n');
fprintf('seg |    u    |  P_u  | Du/P_u | Gam  | P_th  | a [ft/s2] | T [s]\n');
for k = 1:N_TRIM-1
    fprintf('%3d | %7.1f | %5.2f | %6.3f | %4.2f | %5.3f | %9.2f | %5.2f\n', ...
            k, u_anchor(k), tab(k,1), tab(k,2), tab(k,3), tab(k,4), tab(k,5), Traw(k));
end
fprintf('total at eta = 1 : %.2f s   (nothing this short has ever flown)\n', sum(Traw));

LOW = 1:4;  MID = 5:8;  HIGH = 9:19;
eta = mean(Traw(LOW)) / 2.5;
Tcal = Traw / eta;
fprintf(['\neta calibrated on the low band alone : %.3f\n' ...
         '  band means with that eta   low %.2f   mid %.2f   high %.2f  (s)\n' ...
         '  flight-certified bands     low 2.50   mid 1.80   high 1.50  (s)\n' ...
         '  total %.1f s   vs the 33.7 s found by flight search\n'], ...
        eta, mean(Tcal(LOW)), mean(Tcal(MID)), mean(Tcal(HIGH)), sum(Tcal));

Trep = law_schedule(half, u_anchor, lam_r/G, DTH, N_TRIM);
Trep = Trep / (mean(Trep(LOW))/2.5);

%% ---------------------------------------------------------------------
%  2. fly it, with the band schedule as the control
%% ---------------------------------------------------------------------
band = [repmat(2.5,1,4) repmat(1.8,1,4) repmat(1.5,1,11)];
cases = { ...
  'band 2.5/1.8/1.5 (control)',  band; ...
  'law, eta from low band',      Tcal; ...
  'law, total matched to 33.7',  Traw * (33.7/sum(Traw)); ...
  'law, scaled 2.30',            Traw * 2.30; ...
  'law, scaled 2.60',            Traw * 2.60; ...
  'law 1-12 + band 1.5 for 13-19', [Tcal(1:12) repmat(1.5,1,7)]; ...
  'band 1-12 + law for 13-19',     [band(1:12) Tcal(13:19)]; ...
  'law, lambda repaired at 17-20', Trep; ...
  'law, floored at 1.5 s',         max(Tcal, 1.5); ...
};

fprintf('\n--- flight ---\n');
fprintf('%-32s | %6s | %7s | %8s | %7s | %s\n', ...
        'schedule', 'time', 'V>=0 %', 'worst V', 'at u', 'verdict');
fprintf('%s\n', repmat('-', 1, 82));
res = struct([]);
for i = 1:size(cases,1)
    [tt, uu] = profile_from(u_anchor, cases{i,2}, dt, T_HOLD);
    R = fly(params, tt, uu, dt, brtV, T_HOLD, half, u_anchor);
    if R.viol == 0, vd = 'SAFE'; else, vd = 'unsafe'; end
    fprintf('%-32s | %6.1f | %7.1f | %+8.3f | %7.1f | %s\n', ...
            cases{i,1}, tt(end)-2*T_HOLD, R.viol, R.Vw, R.u_worst, vd);
    res(i).name = cases{i,1};   res(i).Ts   = cases{i,2};
    res(i).t    = tt(end)-2*T_HOLD;  res(i).viol = R.viol;  res(i).Vw = R.Vw;
    res(i).u_worst = R.u_worst;  res(i).V = R.V;  res(i).tt = tt;  res(i).uu = uu;
    res(i).use  = R.use;
end

%% ---------------------------------------------------------------------
%  3. which axis is being spent
%% ---------------------------------------------------------------------
fprintf('\n--- axis usage (peak deviation / tube half-width, %%) ---\n');
fprintf('%-32s |', 'schedule');  fprintf(' %14s', AX{:});  fprintf('\n');
for i = 1:numel(res)
    fprintf('%-32s |', res(i).name);  fprintf(' %13.1f%%', res(i).use);  fprintf('\n');
end
fprintf(['\nThe law does what it advertises: it drives pitch attitude to its\n' ...
         'budget (86-95 %%) where the band schedule stops at 74 %%. The tube is\n' ...
         'still left, because the tube is 4-D and vertical speed is holding\n' ...
         '69 %% at the same time - budget the law assumes is free.\n']);

%% ---------------------------------------------------------------------
%  4. figure
%% ---------------------------------------------------------------------
SHOW = [1 2 6 7 8];
f = figure('Position',[100 100 980 760]);
subplot(3,1,1); hold on;
for i = SHOW, plot(res(i).tt, res(i).uu, 'LineWidth', 1.1); end
ylabel('commanded speed [ft/s]'); grid on;
legend({res(SHOW).name}, 'Location','southeast', 'FontSize', 7);
title('the closed-form acceleration law vs the flight-found band schedule');

subplot(3,1,2); hold on;
plot(u_anchor(1:end-1), band, 'k-o', 'LineWidth', 1.4);
plot(u_anchor(1:end-1), Tcal, 'r-s', 'LineWidth', 1.4);
plot(u_anchor(1:end-1), Trep, 'b-^', 'LineWidth', 1.0);
xlabel('anchor speed [ft/s]'); ylabel('segment time [s]'); grid on;
legend({'band schedule','the law','law, lambda repaired'}, 'Location','northeast', 'FontSize', 8);

subplot(3,1,3); hold on;
for i = SHOW, plot(res(i).tt, res(i).V, 'LineWidth', 1.1); end
yline(0,'k-','LineWidth',1.4); ylim([-0.5 0.5]);
xlabel('time [s]'); ylabel('BRT value (V<0 = inside)'); grid on;
saveas(f, fullfile(root,'logger','a_law_check.png'));
save(fullfile(root,'logger','a_law_check.mat'), 'res','Traw','Tcal','Trep','eta','half','lam');
fprintf('\nsaved logger/a_law_check.{mat,png}\n');

% =========================================================================
function [Tk, tab] = law_schedule(half, ua, mu_th, dth, N)
% 법칙 본체. 구간마다 나눗셈 한 번, 반복도 예측도 없음.
%
% 입력
%   half  : 4 x 20. 안전 영역의 반폭. half(1,k)=속도 방향, half(4,k)=자세각 방향.
%           check_tube_dimensions.m 이 트림점에서 한 축씩 걸어 나가며 잰 값.
%   ua    : 1 x 20. 트림점 20개의 전진속도 [ft/s]. 0 부터 160.3 까지 8.44 간격.
%   mu_th : 1 x 20. "가속도 1 ft/s^2 당 숙여야 하는 각도" [rad]. = lambda/g.
%           호버에서 1/32.174 = 0.0311 rad (1.78도).
%   dth   : 트림점 하나 건널 때 목표 자세각 자체가 움직이는 양 [rad]. 0.01047.
%   N     : 트림점 개수 (20)
% 출력
%   Tk    : 1 x 19. 각 구간에 줄 시간 [s].
%   tab   : 중간값 기록 (표로 찍기 위한 것, 계산에는 안 쓰임)

Tk = zeros(1, N-1);  tab = zeros(N-1, 5);
for k = 1:N-1
    % 이번 구간에서 올려야 할 속도. 항상 8.44 ft/s.
    Du   = ua(k+1) - ua(k);

    % 양 끝 트림점의 여유를 합친 값. 두 안전 영역이 겹쳐 있으니 건너가는 동안
    % 양쪽 방을 다 쓸 수 있다는 뜻.
    P_u  = half(1,k) + half(1,k+1);      % 속도 방향으로 쓸 수 있는 총 폭
    P_th = half(4,k) + half(4,k+1);      % 자세각 방향으로 쓸 수 있는 총 폭

    % 속도 방향 관문. 걸음(8.44)이 합친 폭보다 크면 두 영역이 속도 방향으로
    % 아예 안 겹친다는 뜻이라 시간을 아무리 줘도 안 됨 -> 트림점을 더 촘촘히
    % 놓거나 안전 영역을 다시 계산해야 함. (실제로는 20구간 전부 통과함)
    rat  = Du / P_u;
    if rat >= 1
        error('check_a_law:gate', 'speed-axis gate fails at segment %d', k);
    end

    % 속도 방향이 빠듯할수록 커지는 값. 자세각 예산을 깎으라는 신호.
    % 다만 곱해지는 dth 가 0.01 로 작아서 실제로는 P_th 의 4 %밖에 못 깎음.
    Gam = 1/sqrt(1 - rat^2);

    % 허용 가속도. 두 가지 읽기 중 더 엄한 쪽을 택함.
    %   앞 : 합친 자세각 방을 양 끝이 나눠 쓴다고 보고 (그래서 2로 나눔)
    %   뒤 : 더 좁은 쪽 영역 하나를 아예 못 벗어난다고 보고
    % 어느 쪽이든 나눗셈의 뜻은 같음 -- "쓸 수 있는 각도"를 "가속도 1당
    % 숙이는 각도"로 나누면 쓸 수 있는 가속도가 나온다.
    a   = min( (P_th - Gam*dth)/(2*mu_th(k)), ...
               min(half(4,k), half(4,k+1))/mu_th(k) );

    % 그 가속도로 8.44 ft/s 를 올리는 데 걸리는 시간.
    Tk(k)    = Du / a;
    tab(k,:) = [P_u, rat, Gam, P_th, a];
end
end

function [t, u] = profile_from(us, Ts, dt, T_hold)
u = repmat(us(1), 1, round(T_hold/dt));
for k = 1:numel(Ts)
    n = max(round(Ts(k)/dt), 1);
    u = [u, us(k) + (us(k+1)-us(k))*(0:n-1)/n]; %#ok<AGROW>
end
u = [u, repmat(us(end), 1, round(T_hold/dt)+1)];
t = (0:numel(u)-1)*dt;
end

function R = fly(params, tt, uu, dt, brtV, T_hold, half, u_anchor)
guam = LpC_GUAM(Config('trim_schedule', params));  guam.reset();
N = numel(tt);  pos = cumtrapz(tt, uu);
V = zeros(1,N);  dev = zeros(4,N);
for k = 1:N
    x  = guam.state([4 6 11 8]);
    dev(:,k) = x(:) - reshape(brtV.trim_at(uu(k)), [], 1);
    V(k) = brtV.value(x, uu(k));
    guam.step(struct('pos',[pos(k);0;-80], 'vel',[uu(k);0;0], 'chi',0, 'chi_dot',0));
end
w = round(T_hold/dt)+1 : N-round(T_hold/dt);
R.V = V;  [R.Vw, iw] = max(V(w));  R.u_worst = uu(w(iw));
R.viol = 100*sum(V(w) >= 0)/numel(w);
R.use  = zeros(1,4);
for d = 1:4
    [pk, ip] = max(abs(dev(d,w)));
    allow = interp1(u_anchor, half(d,:), ...
                    min(max(uu(w(ip)), u_anchor(1)), u_anchor(end)));
    R.use(d) = 100*pk/allow;
end
end
