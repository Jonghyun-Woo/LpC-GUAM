function plot_overlap_segments(varargin)
% PLOT_OVERLAP_SEGMENTS  28.93초 스케줄의 19개 구간을 (du, dth) 평면에 그린다.
%
%   구간 k 마다 BRT_k 와 BRT_k+1 의 단면 경계, 그리고 기체가 실제로 지나간
%   궤적을 함께 그린다.  두 BRT 가 겹치는 곳을 궤적이 통과하는지 눈으로 볼 수
%   있다.
%
%   실행
%       matlab -batch "addpath(genpath('.')); plot_overlap_segments"
%
%   선택 인자 (이름-값 쌍)
%       'save'    true      logger/ 에 png 로 저장할지
%       'refly'   false     궤적을 다시 비행해서 얻을지.  false 면 저장된 것을 쓴다
%       'plane'   'theta'   'theta' | 'w' | 'q'   세로축으로 쓸 성분
%
%   그림 구성
%       그림 1 : 구간 1~6    (2 x 3)
%       그림 2 : 구간 7~12   (2 x 3)
%       그림 3 : 구간 13~18  (2 x 3)
%       그림 4 : 구간 19     (단독)
%
%   그림 읽는 법
%       청록 실선   앵커 k 의 BRT 경계
%       주황 파선   앵커 k+1 의 BRT 경계.  가로축을 Du 만큼 옮겨 놓았으므로
%                   두 곡선이 겹치는 가로 구간이 곧 겹침 구간이다
%       검은 선     실제 비행 궤적.  흰 원이 시작, 검은 삼각형이 끝
%       회색 영역   두 BRT 가 세로축 값을 공유하며 겹치는 가로 범위
%
%       각 BRT 단면은 나머지 두 성분을 그 구간 중간의 실측값으로 고정해서
%       자른 것이다.  4차원 덩어리를 2차원으로 보려면 어딘가를 잘라야 하고,
%       기체가 실제로 지나는 자리에서 자르는 것이 가장 뜻이 있다.

OPT = struct('save', true, 'refly', false, 'plane', 'theta');
for i = 1:2:numel(varargin), OPT.(varargin{i}) = varargin{i+1}; end

root = fileparts(fileparts(mfilename('fullpath')));

% 세로축으로 쓸 성분 고르기.  편차 벡터의 순서는 [du; dw; dq; dth]
switch lower(OPT.plane)
    case 'theta', JY = 4;  YLAB = '\delta\theta  [deg]';   YSC = 180/pi;
    case 'w',     JY = 2;  YLAB = '\deltaw  [ft/s]';       YSC = 1;
    case 'q',     JY = 3;  YLAB = '\deltaq  [deg/s]';      YSC = 180/pi;
    otherwise, error('plane 은 theta | w | q 중 하나여야 합니다');
end

%% ---- 준비물 ---------------------------------------------------------
T   = load('trim_table_Poly_ConcatVer4p0.mat');
UH  = T.UH(1:20).';                              % 앵커 속도 [ft/s]
TL  = T.XU0_interp([1 3 5 11], 1:20, 3);         % [u; w; q; th] 지평선 기준
du  = diff(UH);                                  % 트림 간격

S = load(fullfile(root,'logger','overlap_schedule.mat'));
a = S.a;                                         % 28.93초 스케줄

spec  = FilterConfig.channelSpec('lon');
GMIN  = spec.grid_min;  GMAX = spec.grid_max;
gv = cell(1,4);
for i = 1:4
    gv{i} = linspace(spec.grid_min(i), spec.grid_max(i), spec.grid_num(i));
end
V = cell(1,20);
for k = 1:20
    Sf = load(fullfile(root,'data', ...
        sprintf('GUAM_LON_BRT_HJIR_UH%d_WH%d.mat', k, 3)), 'data');
    V{k} = griddedInterpolant(gv, Sf.data, 'linear', 'nearest');
end

%% ---- 궤적 ----------------------------------------------------------
trjfile = fullfile(root,'logger','overlap_traj.mat');
if OPT.refly || ~isfile(trjfile)
    TRJ = fly_schedule(UH, TL, du, a);
    save(trjfile, 'TRJ');
else
    L = load(trjfile);  TRJ = L.TRJ;
end

%% ---- 그림 ----------------------------------------------------------
GRP = {1:6, 7:12, 13:18, 19};                    % 그림 4장으로 나눈다
for gi = 1:4
    ks = GRP{gi};
    if numel(ks) == 1
        f = figure('Position',[120 120 560 460],'Color','w');
        nr = 1;  nc = 1;
    else
        f = figure('Position',[60 60 1280 720],'Color','w');
        nr = 2;  nc = 3;
    end
    for ii = 1:numel(ks)
        ax = subplot(nr, nc, ii);
        h = draw_segment(ax, ks(ii), TRJ, V, TL, du, a, GMIN, GMAX, JY, YLAB, YSC);
        if ii == 1
            legend(ax, h, {sprintf('BRT %d', ks(1)), ...
                        sprintf('BRT %d (+\\Deltau)', ks(1)+1), 'flown'}, ...
                   'Location','southwest','FontSize',8);
        end
    end
    sgtitle(sprintf('Overlap along %s  -  segments %s   (total 28.93 s)', ...
        YLAB, mat2str(ks)));
    if OPT.save
        out = fullfile(root,'logger', ...
            sprintf('overlap_%s_%d.png', lower(OPT.plane), gi));
        exportgraphics(f, out, 'Resolution', 150);
        fprintf('  saved %s\n', out);
    end
end
end


%% ========================================================================
function h = draw_segment(ax, k, TRJ, V, TL, du, a, GMIN, GMAX, JY, YLAB, YSC)
% 구간 k 하나를 그린다.
X  = TRJ{k};                                     % N x 4, 앵커 k 기준 편차
im = max(round(size(X,1)/2), 1);                 % 구간 중간

% 나머지 두 성분을 구간 중간의 실측값으로 고정해 단면을 만든다
oth = setdiff([2 3 4], JY);
fx  = [X(im,oth(1))  X(im,oth(2))];

ug = linspace(GMIN(1), GMAX(1), 181);
yg = linspace(GMIN(JY), GMAX(JY), 181);
[UU, YY] = meshgrid(ug, yg);

q = zeros(4, numel(UU));
q(1,:) = UU(:).';   q(JY,:) = YY(:).';
q(oth(1),:) = fx(1);  q(oth(2),:) = fx(2);
Z1 = reshape(V{k}(q(1,:).', q(2,:).', q(3,:).', q(4,:).'), size(UU));

% 앵커 k+1 : 트림점 차이만큼 편차를 옮겨서 조회한다
sh = TL(2:4,k) - TL(2:4,k+1);
q2 = q;
q2(2,:) = q2(2,:) + sh(1);
q2(3,:) = q2(3,:) + sh(2);
q2(4,:) = q2(4,:) + sh(3);
Z2 = reshape(V{k+1}(q2(1,:).', q2(2,:).', q2(3,:).', q2(4,:).'), size(UU));

hold(ax,'on');  grid(ax,'on');

% 겹치는 가로 범위를 회색으로 칠한다
in1 = any(Z1 <= 0, 1);   in2 = any(Z2 <= 0, 1);
ovl = in1 & circshift_ok(in2, du(k), ug);
if any(ovl)
    xo = ug(ovl);
    yl = [GMIN(JY) GMAX(JY)]*YSC;
    patch(ax, [min(xo) max(xo) max(xo) min(xo)], [yl(1) yl(1) yl(2) yl(2)], ...
        [0.85 0.85 0.85], 'EdgeColor','none', 'FaceAlpha',0.45);
end

[~,h1] = contour(ax, UU,       YY*YSC, Z1, [0 0], 'LineWidth',2, 'LineColor',[0.06 0.43 0.39]);
[~,h2] = contour(ax, UU+du(k), YY*YSC, Z2, [0 0], 'LineWidth',2, 'LineColor',[0.85 0.45 0.10], ...
        'LineStyle','--');
h3 = plot(ax, X(:,1), X(:,JY)*YSC, 'k-', 'LineWidth',2.2);
h = [h1 h2 h3];
plot(ax, X(1,1),   X(1,JY)*YSC,   'ko','MarkerFaceColor','w','MarkerSize',7);
plot(ax, X(end,1), X(end,JY)*YSC, 'k^','MarkerFaceColor','k','MarkerSize',7);

xlabel(ax,'\deltau  [ft/s]');  ylabel(ax, YLAB);
xlim(ax, [-16 22]);
title(ax, sprintf('%d \\rightarrow %d   (a = %.2f, T = %.2f s)', ...
    k, k+1, a(k), du(k)/a(k)), 'FontWeight','normal');
end


function m2 = circshift_ok(in2, dus, ug)
% 앵커 k+1 의 가로 범위를 Du 만큼 옮겨 앵커 k 의 격자에 맞춘다.
x2 = ug + dus;
m2 = interp1(x2, double(in2), ug, 'nearest', 0) > 0.5;
end


%% ========================================================================
function TRJ = fly_schedule(UH, TL, du, a)
% 28.93초 스케줄로 실제 비행하고 구간마다 편차 궤적을 남긴다.
dt = 0.01;  WDOT = 11.667;  ALT0 = 1600;  T_SETTLE = 45.0;

params = struct('filter_mode','off');  params.T_seg = 2.0;
guam = LpC_GUAM(Config('trim_schedule', params));  guam.reset();
st = guam.saveState();  st.state(3) = -ALT0;  guam.restoreState(st);
p = [0; 0; -ALT0];

for i = 1:round(T_SETTLE/dt)                     % 강하 호버에서 정착
    p = p + [0; 0; WDOT*dt];
    guam.step(struct('pos',p,'vel',[0;0;WDOT],'chi',0,'chi_dot',0));
end

TRJ = cell(1,19);
for k = 1:19
    N = max(round(du(k)/a(k)/dt), 1);
    X = zeros(N,4);
    for i = 1:N
        s  = guam.state;  th = s(8);
        uh =  s(4)*cos(th) + s(6)*sin(th);       % 동체 -> 지평선
        wh = -s(4)*sin(th) + s(6)*cos(th);
        X(i,:) = ([uh; wh; s(11); th] - TL(:,k)).';
        v = UH(k) + a(k)*(i*dt);
        p = p + [v*dt; 0; WDOT*dt];
        guam.step(struct('pos',p,'vel',[v;0;WDOT],'chi',0,'chi_dot',0));
    end
    TRJ{k} = X;
end
end
