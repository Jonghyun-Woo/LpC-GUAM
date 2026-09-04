% run_wff_compare_fig - the mission with and without the level-flight
% vertical-speed term.
%
% Left  : w_ff off  (shipped)
% Right : w_ff = u*tan(theta) - w_trim, from the flown attitude
%
% Both panels use the tuned weights q1 = 0.08, q4 = 2.0 and the shipped
% position gain k_pos = 0.1, so the only difference is the new term.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt = 0.01;
MODE = [0 2];
Q    = [0.08 0.01 1000 2.0 0 0]';
params = struct('filter_mode','off');  params.T_seg = 2.0;

traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt, params));
brtV = BRTValue(fullfile(root,'data'), traj.trim_lon, traj.wh_idx);
g0   = LpC_GUAM(Config('trim_schedule', params));
bc0  = g0.controller.baseline_controller;
cfg  = RSLQRConfig;

N = cfg.N_trim;  M = cfg.M_trim;
Ki = zeros(3,3,N,M);  Kx = zeros(3,4,N,M);
for jj = 1:M
    for ii = 1:N
        lon = bc0.ctrl_lon(bc0.trim_xu_eq(bc0.XU0(:,ii,jj)), Q, cfg.Rlon0, cfg.Wlon0);
        Ki(:,:,ii,jj) = lon.Ki;  Kx(:,:,ii,jj) = lon.Kx;
    end
end

f = figure('Name','level-flight w term','Position',[40 40 1500 780]);
ax = gobjects(1,2);
fprintf('\nq1 = 0.08, q4 = 2.0, k_pos = 0.1\n\n');
fprintf('%-22s %9s %9s %9s %10s %8s\n', 'w feedforward', ...
        'max|u-r|','e_pos','alt dev','alt swing','BRT>=0');
fprintf('%s\n', repmat('-',1,72));

for p = 1:2
    guam = LpC_GUAM(Config('trim_schedule', params));
    bc = guam.controller.baseline_controller;
    bc.LON.Ki = Ki;  bc.LON.Kx = Kx;  bc.w_ff_mode = MODE(p);
    guam.reset();

    rt = guam.refTraj;  n = size(rt.pos,2);
    ST = zeros(12,n);  V = zeros(1,n);  EP = zeros(1,n);  R = rt.vel(1,:);
    for i = 1:n
        s = guam.state;
        Rm = RSLQR.rotm_i2b(s(7), s(8), s(9));
        e  = Rm * (rt.pos(:,i) - s(1:3));
        ST(:,i)=s;  EP(i)=e(1);
        V(i) = brtV.value(s([4 6 11 8]), R(i));
        guam.step(struct('pos',rt.pos(:,i),'vel',rt.vel(:,i), ...
                         'chi',rt.chi(i),'chi_dot',rt.chidot(i)));
    end
    t = rt.time;  keep = t <= traj.t_node(end) + 3;
    t=t(keep); ST=ST(:,keep); V=V(keep); R=R(keep); EP=EP(keep);

    alt  = -ST(3,:);   late = t >= 15;
    eu   = max(abs(ST(4,:) - R));
    ep   = max(abs(EP));
    ad   = max(abs(alt - 80));
    sw   = max(alt(late)) - min(alt(late));
    viol = 100*nnz(V >= 0)/numel(V);
    if MODE(p)==0, nm = 'off (shipped)'; else, nm = 'u*tan(theta) - w_trim'; end
    fprintf('%-22s %9.2f %9.1f %9.1f %10.1f %8.1f\n', nm, eu, ep, ad, sw, viol);

    ax(p) = subplot(1,2,p); hold on; grid on;
    cmap = parula(256);
    for k = 1:traj.n_trim
        tk = traj.t_node(k);
        if tk > t(end), continue; end
        [~, i] = min(abs(t - tk));
        col = cmap(round((k-1)/(traj.n_trim-1)*255)+1, :);
        C = brt_zero_contour(brtV, R(i), ST([4 6 11 8], i));
        for c = 1:numel(C)
            q = C{c};
            plot3(ax(p), q(1,:), tk*ones(1,size(q,2)), q(2,:), ...
                  'Color',[col 0.35],'LineWidth',0.9,'HandleVisibility','off');
        end
    end
    plot3(ax(p), nan,nan,nan,'-','Color',[.35 .55 .55],'LineWidth',0.9, ...
          'DisplayName','BRT of each trim point');

    xr = zeros(2,numel(t));
    for i = 1:numel(t)
        e = brtV.trim_at(R(i));  xr(:,i) = [e(1); rad2deg(e(4))];
    end
    plot3(ax(p), xr(1,:), t, xr(2,:), '--','Color',[.1 .1 .1], ...
          'LineWidth',1.8,'DisplayName','plan  x_e(r)');
    inside = V < 0;
    plot_runs(ax(p), t, ST,  inside, [0 .6 .25], 2.2, 'flown - inside BRT');
    plot_runs(ax(p), t, ST, ~inside, [.9 .1 .1], 2.6, 'flown - OUTSIDE BRT');

    xlabel(ax(p),'u [ft/s]'); ylabel(ax(p),'time [s]'); zlabel(ax(p),'\theta [deg]');
    title(ax(p), sprintf(['w feedforward: %s\n' ...
          'speed %.2f ft/s | position %.1f ft | altitude %.1f ft | ' ...
          'alt swing %.1f ft | outside BRT %.1f %%'], nm, eu, ep, ad, sw, viol));
    view(ax(p), -60, 22);  ylim(ax(p), [0 t(end)]);
    if p==1, legend(ax(p),'Location','northeast'); end
end
xl = [min(arrayfun(@(a) a.XLim(1),ax)), max(arrayfun(@(a) a.XLim(2),ax))];
zl = [min(arrayfun(@(a) a.ZLim(1),ax)), max(arrayfun(@(a) a.ZLim(2),ax))];
set(ax,'XLim',xl,'ZLim',zl);

out = fullfile(root,'logger','wff_compare');
savefig(f,[out '.fig']);  exportgraphics(f,[out '.png'],'Resolution',140);
fprintf('\nsaved %s(.fig/.png)\n', out);

% =========================================================================
function plot_runs(ax, t, ST, mask, col, lw, name)
d = diff([false, mask(:)', false]);
s = find(d==1);  e = find(d==-1)-1;  shown = false;
for j = 1:numel(s)
    idx = s(j):e(j);
    if numel(idx) < 2, continue; end
    if shown
        plot3(ax, ST(4,idx), t(idx), rad2deg(ST(8,idx)),'-','Color',col, ...
              'LineWidth',lw,'HandleVisibility','off');
    else
        plot3(ax, ST(4,idx), t(idx), rad2deg(ST(8,idx)),'-','Color',col, ...
              'LineWidth',lw,'DisplayName',name);  shown = true;
    end
end
if ~shown
    plot3(ax,nan,nan,nan,'-','Color',col,'LineWidth',lw,'DisplayName',name);
end
end
