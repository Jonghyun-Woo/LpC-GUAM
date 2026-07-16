function figs = plot_on_state_input_compare(R_or_trace, opts)
% plot_on_state_input_compare
%
% Compare u_nom, u0, and u at the SAME ON-run states.
%
% Meaning:
%   u_nom : nominal longitudinal input before liveness filter
%   u0    : box-clipped nominal input
%   u     : filtered input after QP
%
% This is NOT an OFF-vs-ON trajectory comparison.
% This is a same-state input comparison along the ON trajectory.

    if nargin < 2 || isempty(opts)
        opts = struct();
    end

    opts = fill_input_compare_defaults(opts);

    if isfield(R_or_trace, 'trace')
        tr = R_or_trace.trace;
    else
        tr = R_or_trace;
    end

    k = tr.k(:);
    n = numel(k);

    valid = isfinite(k);

    k = k(valid);

    u_nom = tr.u_nom(valid, :);
    u0    = tr.u0(valid, :);
    u     = tr.u(valid, :);
    lb    = tr.lb(valid, :);
    ub    = tr.ub(valid, :);

    inputNames = {'Pi1','Pi2','Pi3','Pi4','Pi5','Pi6','Pi7','Pi8','Pi9','de','df'};
    unitNames  = {'RPM','RPM','RPM','RPM','RPM','RPM','RPM','RPM','RPM','deg','deg'};

    % Surface channels are assumed to be in rad internally.
    % Convert elevator/flap perturbations to deg for plotting.
    scale = ones(1, size(u, 2));
    scale(10:11) = 180 / pi;

    u_nom_p = u_nom .* scale;
    u0_p    = u0    .* scale;
    u_p     = u     .* scale;
    lb_p    = lb    .* scale;
    ub_p    = ub    .* scale;

    figs = struct();

    % ---------------------------------------------------------------------
    % Figure 1: overall modification and CBF-rate diagnostics
    % ---------------------------------------------------------------------
    figs.summary = figure('Name', 'same_ON_state_input_summary', 'Color', 'w');
    tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    plot(k, tr.satClip(valid), 'LineWidth', 1.2); hold on;
    plot(k, tr.du(valid), 'LineWidth', 1.2);
    plot(k, tr.cmdChange(valid), 'LineWidth', 1.2);
    xlabel('step k');
    ylabel('norm');
    title('Input modification norms at same ON states');
    legend({'||u_{nom}-u0|| clipping', ...
            '||u-u0|| QP move', ...
            '||u-u_{nom}|| total'}, ...
            'Location', 'best');
    grid on;

    nexttile;
    plot(k, tr.dV_raw(valid), 'LineWidth', 1.0); hold on;
    plot(k, tr.dV_box(valid), 'LineWidth', 1.0);
    plot(k, tr.dV_flt(valid), 'LineWidth', 1.0);
    plot(k, tr.rhs(valid), 'k--', 'LineWidth', 1.0);
    yline(0, ':');
    xlabel('step k');
    ylabel('value rate');
    title('Value-rate comparison at same ON states');
    legend({'dV raw nominal', 'dV box nominal', 'dV filtered', 'rhs'}, ...
           'Location', 'best');
    grid on;

    nexttile;
    stairs(k, tr.active(valid), 'LineWidth', 1.1); hold on;
    stairs(k, tr.rawBoxOK(valid), 'LineWidth', 1.1);
    stairs(k, tr.boxCBFOK(valid), 'LineWidth', 1.1);
    stairs(k, tr.feasLin(valid), 'LineWidth', 1.1);
    stairs(k, tr.insideGrid(valid), 'LineWidth', 1.1);
    ylim([-0.1, 1.1]);
    xlabel('step k');
    ylabel('logical');
    title('Filter / feasibility flags');
    legend({'active', 'raw box OK', 'box CBF OK', 'linear feasible', 'inside grid'}, ...
           'Location', 'best');
    grid on;

    % ---------------------------------------------------------------------
    % Figure 2: per-channel input comparison
    % ---------------------------------------------------------------------
    figs.channels = figure('Name', 'same_ON_state_input_channels', 'Color', 'w');
    tiledlayout(4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    nu = size(u, 2);

    for i = 1:nu
        nexttile;

        plot(k, u_nom_p(:, i), '--', 'LineWidth', 0.9); hold on;
        plot(k, u0_p(:, i), ':', 'LineWidth', 1.2);
        plot(k, u_p(:, i), '-', 'LineWidth', 1.2);
        plot(k, lb_p(:, i), 'k--', 'LineWidth', 0.7);
        plot(k, ub_p(:, i), 'k--', 'LineWidth', 0.7);

        xlabel('k');
        ylabel(sprintf('%s [%s]', inputNames{i}, unitNames{i}));
        title(inputNames{i});
        grid on;

        if i == 1
            legend({'u_{nom}', 'u0 clipped', 'u filtered', 'lb', 'ub'}, ...
                   'Location', 'best');
        end
    end

    % ---------------------------------------------------------------------
    % Figure 3: pure differences
    % ---------------------------------------------------------------------
    figs.delta = figure('Name', 'same_ON_state_input_deltas', 'Color', 'w');
    tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    plot(k, u_nom_p(:, 1:9) - u0_p(:, 1:9), 'LineWidth', 1.0);
    xlabel('step k');
    ylabel('RPM');
    title('Rotor clipping: u_{nom} - u0');
    legend(inputNames(1:9), 'Location', 'eastoutside');
    grid on;

    nexttile;
    plot(k, u_p(:, 1:9) - u0_p(:, 1:9), 'LineWidth', 1.0);
    xlabel('step k');
    ylabel('RPM');
    title('Rotor QP correction: u - u0');
    legend(inputNames(1:9), 'Location', 'eastoutside');
    grid on;

    nexttile;
    plot(k, [u_nom_p(:,10:11) - u0_p(:,10:11), ...
             u_p(:,10:11)     - u0_p(:,10:11)], ...
         'LineWidth', 1.1);
    xlabel('step k');
    ylabel('deg');
    title('Surface clipping/QP correction');
    legend({'de nom-u0', 'df nom-u0', 'de u-u0', 'df u-u0'}, ...
           'Location', 'best');
    grid on;

    % ---------------------------------------------------------------------
    % Console summary
    % ---------------------------------------------------------------------
    fprintf('\n=== Same-ON-state input comparison ===\n');
    fprintf('n steps             : %d\n', n);
    fprintf('active count         : %d\n', sum(tr.active(valid) == 1, 'omitnan'));
    fprintf('max ||u_nom-u0||     : %.6g\n', max(tr.satClip(valid), [], 'omitnan'));
    fprintf('max ||u-u0||         : %.6g\n', max(tr.du(valid), [], 'omitnan'));
    fprintf('max ||u-u_nom||      : %.6g\n', max(tr.cmdChange(valid), [], 'omitnan'));

    idxFirst = find((tr.active(valid) == 1) | ...
                    (tr.du(valid) > 1e-9) | ...
                    (tr.satClip(valid) > 1e-9), 1, 'first');

    if ~isempty(idxFirst)
        fprintf('first modification k : %d\n', k(idxFirst));
        fprintf('  V        = %.6g\n', tr.V(valid_index(valid, idxFirst)));
        fprintf('  rhs      = %.6g\n', tr.rhs(valid_index(valid, idxFirst)));
        fprintf('  dV_box   = %.6g\n', tr.dV_box(valid_index(valid, idxFirst)));
        fprintf('  dV_flt   = %.6g\n', tr.dV_flt(valid_index(valid, idxFirst)));
    end
end