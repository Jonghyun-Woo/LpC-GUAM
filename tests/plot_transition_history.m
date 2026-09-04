function plot_transition_history(varargin)
% PLOT_TRANSITION_HISTORY  28.93초 천이를 시간축으로 보여 준다.
%
%   제자리 비행에서 순항 비행으로 넘어가는 동안 무엇이 어떻게 바뀌는지를
%   한 장에 담는다.  로터 사용률과 피치각이 순항 트림값에 안착하는 것이
%   천이가 실제로 끝났다는 증거가 된다.
%
%   실행
%       matlab -batch "addpath(genpath('.')); plot_transition_history"
%
%   선택 인자 (이름-값 쌍)
%       'save'    true   logger/transition_history.png 로 저장할지
%       'hover'   20     천이 전 제자리 비행을 몇 초 보여 줄지
%       'cruise'  40     천이 후 순항을 몇 초 보여 줄지
%
%   그림 구성 (세로 5칸, 가로축은 천이 시작을 0 으로 잡은 시각)
%       1  전진속도        지령과 실제
%       2  피치각          실제와 트림값
%       3  리프트 로터     8개 각각, 최대 회전수 대비 사용률
%       4  푸셔 로터       사용률
%       5  구동기 여유     한계까지 얼마나 남았나
%
%   읽는 법
%       세로 점선 두 개가 천이의 시작과 끝이다.  그 사이에서 리프트 로터가
%       내려가고 푸셔가 올라가며, 끝난 뒤에는 세 값이 모두 순항 트림값
%       (가로 점선) 에 붙어 움직이지 않는다.

OPT = struct('save', true, 'hover', 20, 'cruise', 40);
for i = 1:2:numel(varargin), OPT.(varargin{i}) = varargin{i+1}; end

root = fileparts(fileparts(mfilename('fullpath')));

%% ---- 준비물 ---------------------------------------------------------
T  = load('trim_table_Poly_ConcatVer4p0.mat');
UH = T.UH(1:20).';
du = diff(UH);
TH_TRIM = T.XU0_interp(11, 1:20, 3);             % 트림 피치각 [rad]

S = load(fullfile(root,'logger','overlap_schedule.mat'));
a = S.a;

RPM = 60/(2*pi);
LIFT_TRIM_20 = mean(T.XU0_interp(17:24, 20, 3)) * RPM;   % 순항 리프트 트림
PUSH_TRIM_20 =      T.XU0_interp(25,    20, 3)  * RPM;   % 순항 푸셔 트림
TH_TRIM_20   =      TH_TRIM(20) * 180/pi;

cfg    = RSLQRConfig();
ENGMAX = cfg.eng_max(:).' * RPM;                 % 로터 최대 회전수 [RPM]

%% ---- 비행 ----------------------------------------------------------
[tt, LOG, Ttr] = fly(UH, du, a, OPT.hover, OPT.cruise);
%   LOG 열 : 1 u_body, 2 theta[deg], 3 u_cmd, 4..11 리프트, 12 푸셔  [RPM]

u    = LOG(:,1);      th   = LOG(:,2);   ucmd = LOG(:,3);
lift = LOG(:,4:11);   push = LOG(:,12);
liftPct = 100 * lift ./ ENGMAX(1:8);
pushPct = 100 * push ./ ENGMAX(9);

%% ---- 그림 ----------------------------------------------------------
f = figure('Position',[60 40 1000 940],'Color','w');
CT = [0.06 0.43 0.39];   CO = [0.85 0.45 0.10];  CG = [0.55 0.55 0.55];

ax1 = subplot(5,1,1);  hold(ax1,'on');  grid(ax1,'on');
hc = plot(ax1, tt, ucmd, '--', 'LineWidth',1.6, 'Color',CG);
hf = plot(ax1, tt, u,    '-',  'LineWidth',2.0, 'Color',CT);
ylabel(ax1,'u  [ft/s]');
legend(ax1,[hc hf],{'commanded','flown'},'Location','southeast');
title(ax1, sprintf('Transition history   (transition = %.2f s)', Ttr),'FontWeight','normal');

ax2 = subplot(5,1,2);  hold(ax2,'on');  grid(ax2,'on');
plot(ax2, tt, th, '-', 'LineWidth',2.0, 'Color',CT);
yline(ax2, TH_TRIM_20, ':', 'LineWidth',1.4, 'Color',CO);
ylabel(ax2,'\theta  [deg]');
text(ax2, tt(end)-2, TH_TRIM_20+1.5, 'cruise trim', 'Color',CO,'FontSize',8, ...
     'HorizontalAlignment','right');

ax3 = subplot(5,1,3);  hold(ax3,'on');  grid(ax3,'on');
plot(ax3, tt, liftPct, 'LineWidth',1.1);
yline(ax3, 100*LIFT_TRIM_20/mean(ENGMAX(1:8)), ':', 'LineWidth',1.4, 'Color',CO);
ylabel(ax3,'lift rotors  [% of max]');  ylim(ax3,[0 100]);

ax4 = subplot(5,1,4);  hold(ax4,'on');  grid(ax4,'on');
plot(ax4, tt, pushPct, '-', 'LineWidth',2.0, 'Color',CT);
yline(ax4, 100*PUSH_TRIM_20/ENGMAX(9), ':', 'LineWidth',1.4, 'Color',CO);
ylabel(ax4,'pusher  [% of max]');  ylim(ax4,[0 100]);

ax5 = subplot(5,1,5);  hold(ax5,'on');  grid(ax5,'on');
mgn = 100 - max([liftPct pushPct], [], 2);       % 한계까지 남은 여유
plot(ax5, tt, mgn, '-', 'LineWidth',2.0, 'Color',CT);
ylabel(ax5,'actuator margin  [%]');  xlabel(ax5,'time from transition start  [s]');
ylim(ax5,[0 100]);

for ax = [ax1 ax2 ax3 ax4 ax5]
    x1 = xline(ax, 0,   'k--','LineWidth',1.2);
    x2 = xline(ax, Ttr, 'k--','LineWidth',1.2);
    x1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    x2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    xlim(ax, [tt(1) tt(end)]);
end

if OPT.save
    out = fullfile(root,'logger','transition_history.png');
    exportgraphics(f, out, 'Resolution', 150);
    fprintf('  saved %s\n', out);
end

%% ---- 요약 표 -------------------------------------------------------
seg = {'hover', tt < -1; ...
       'transition start', tt > 1 & tt < 3; ...
       'transition mid',   abs(tt - Ttr/2) < 1; ...
       'transition end',   tt > Ttr-2 & tt < Ttr; ...
       'cruise (>20 s)',   tt > Ttr + 20};
fprintf('\n  phase              lift[%%]  push[%%]  theta[deg]  u[ft/s]\n');
for i = 1:size(seg,1)
    m = seg{i,2};
    fprintf('  %-18s %7.0f %8.0f %11.2f %8.1f\n', seg{i,1}, ...
        mean(mean(liftPct(m,:))), mean(pushPct(m)), mean(th(m)), mean(u(m)));
end
fprintf('\n  cruise trim :  lift %.0f%%   push %.0f%%   theta %.2f deg\n', ...
    100*LIFT_TRIM_20/mean(ENGMAX(1:8)), 100*PUSH_TRIM_20/ENGMAX(9), TH_TRIM_20);
fprintf('  actuator margin, minimum over the whole run : %.0f%%\n', min(mgn));
end


%% ========================================================================
function [tt, LOG, Ttr] = fly(UH, du, a, T_HOVER, T_CRUISE)
% 제자리 비행 - 천이 - 순항을 이어서 날리고 기록한다.
dt   = 0.01;   WDOT = 11.667;   ALT0 = 4000;   T_SETTLE = 45.0;
RPM  = 60/(2*pi);

params = struct('filter_mode','off');  params.T_seg = 2.0;
guam = LpC_GUAM(Config('trim_schedule', params));  guam.reset();
st = guam.saveState();  st.state(3) = -ALT0;  guam.restoreState(st);
p = [0; 0; -ALT0];

for i = 1:round(T_SETTLE/dt)                     % 화면에 안 보이는 정착 구간
    p = p + [0; 0; WDOT*dt];
    guam.step(struct('pos',p,'vel',[0;0;WDOT],'chi',0,'chi_dot',0));
end

LOG = [];   tt = [];   t = -T_HOVER;

for i = 1:round(T_HOVER/dt)                      % 천이 전 제자리 비행
    p = p + [0; 0; WDOT*dt];
    [eng,~,~,~] = guam.step(struct('pos',p,'vel',[0;0;WDOT],'chi',0,'chi_dot',0));
    s = guam.state;
    LOG(end+1,:) = [s(4), s(8)*180/pi, 0, eng(:).'*RPM];   %#ok<AGROW>
    tt(end+1) = t;  t = t + dt;                            %#ok<AGROW>
end

for k = 1:19                                     % 천이
    N = max(round(du(k)/a(k)/dt), 1);
    for i = 1:N
        v = UH(k) + a(k)*(i*dt);
        p = p + [v*dt; 0; WDOT*dt];
        [eng,~,~,~] = guam.step(struct('pos',p,'vel',[v;0;WDOT],'chi',0,'chi_dot',0));
        s = guam.state;
        LOG(end+1,:) = [s(4), s(8)*180/pi, v, eng(:).'*RPM];   %#ok<AGROW>
        tt(end+1) = t;  t = t + dt;                            %#ok<AGROW>
    end
end
Ttr = t;                                         % 천이가 끝난 시각

for i = 1:round(T_CRUISE/dt)                     % 순항 유지
    p = p + [UH(20)*dt; 0; WDOT*dt];
    [eng,~,~,~] = guam.step(struct('pos',p,'vel',[UH(20);0;WDOT],'chi',0,'chi_dot',0));
    s = guam.state;
    LOG(end+1,:) = [s(4), s(8)*180/pi, UH(20), eng(:).'*RPM];  %#ok<AGROW>
    tt(end+1) = t;  t = t + dt;                                %#ok<AGROW>
end
tt = tt(:);
end
