classdef SimLogger < handle
    % Buffers a single closed-loop run and produces plots / a tube-overlay
    % trace. Read-only w.r.t. the simulation: it never advances the plant.
    %
    % Two-call logging contract (per step k, paired, in this order):
    %   addState(k, guam, ref)         % BEFORE guam.step(): pre-step state/aero/ref
    %   addInputs(guam, engine, surface)% AFTER  guam.step(): effectors + filter last_info
    % The filter diagnostics (last_info) are populated during step(), so they
    % can only be read afterwards; the trajectory coordinates use the pre-step
    % state (matching the diagnostic trace convention: statePre + pre-step BRT).
    %
    % Buffers are preallocated to the step count at construction (NaN-filled)
    % and trimmed to the actual count by finalize().
    %
    % Units follow the plant: position ft, velocity ft/s, angles rad, rotor
    % speed rad/s. exportTrace() emits the fields consumed by the LON tube
    % overlay helpers (tests/visualize_lon_tube_overlay_trace.m).

    properties
        cfg         % LoggerConfig
        dt          % [s] time step (from SimConfig)
        N           % buffer capacity (step count)
        count       % number of steps actually recorded
        nu          % filter effector count (LON = 11)
        buf         % struct of preallocated arrays (rows = steps)
        saveDirCached = '';  % resolved output dir (set on first save)
    end

    methods
        function obj = SimLogger(loggerCfg, simConfig)
            obj.cfg   = loggerCfg;
            obj.dt    = simConfig.dt;
            obj.N     = simConfig.steps;
            obj.count = 0;
            obj.nu    = 11;   % LON liveness channel effector count
            obj.init_buffer();
        end

        function init_buffer(obj)
            N_  = obj.N;
            nu_ = obj.nu;
            b = struct();
            % --- basic ---
            b.t       = nan(N_, 1);
            b.k       = nan(N_, 1);
            b.state   = nan(N_, 12);
            b.engine  = nan(N_, 9);
            b.surface = nan(N_, 5);
            b.alpha   = nan(N_, 1);
            b.beta    = nan(N_, 1);
            b.V       = nan(N_, 1);   % airspeed
            b.ref_pos = nan(N_, 3);
            b.ref_vel = nan(N_, 3);
            % --- filter / BRT ---
            b.brtV       = nan(N_, 1);
            b.u_nom      = nan(N_, nu_);
            b.u0         = nan(N_, nu_);
            b.u          = nan(N_, nu_);
            b.lb         = nan(N_, nu_);
            b.ub         = nan(N_, nu_);
            b.active     = nan(N_, 1);
            b.du         = nan(N_, 1);
            b.satClip    = nan(N_, 1);
            b.cmdChange  = nan(N_, 1);
            b.dVdt       = nan(N_, 1);
            b.rhs        = nan(N_, 1);
            b.insideGrid = nan(N_, 1);
            b.ok         = nan(N_, 1);
            obj.buf = b;
        end

        function addState(obj, guam, k, ref)
            % Record pre-step state, aero frame, and reference. Advances the
            % record counter; must be called once per step BEFORE guam.step().
            if ~obj.cfg.enable, return; end
            j = obj.count + 1;
            obj.count = j;

            obj.buf.k(j) = k;
            obj.buf.t(j) = (k - 1) * obj.dt;

            x = guam.state;
            obj.buf.state(j, :) = x(:)';

            [a, bta, Vair] = guam.aeroFrame.compute(x(4), x(5), x(6));
            obj.buf.alpha(j) = a;
            obj.buf.beta(j)  = bta;
            obj.buf.V(j)     = Vair;

            obj.buf.ref_pos(j, :) = ref.pos(:)';
            obj.buf.ref_vel(j, :) = ref.vel(:)';
        end

        function addInputs(obj, guam, engine, surface)
            % Record post-step effectors and filter diagnostics into the record
            % opened by the matching addState(); must be called AFTER guam.step().
            if ~obj.cfg.enable, return; end
            j = obj.count;
            if j < 1, return; end

            obj.buf.engine(j, :)  = engine(:)';
            obj.buf.surface(j, :) = surface(:)';

            if ~obj.cfg.logFilter, return; end
            li = guam.controller.liveness_lon.last_info;
            if isempty(li), return; end

            nu_ = obj.nu;
            obj.buf.brtV(j)       = getfield_default(li, 'V', NaN);
            obj.buf.u_nom(j, :)   = row(getfield_default(li, 'u_nom', nan(nu_, 1)), nu_);
            obj.buf.u0(j, :)      = row(getfield_default(li, 'u0',    nan(nu_, 1)), nu_);
            obj.buf.u(j, :)       = row(getfield_default(li, 'u',     nan(nu_, 1)), nu_);
            obj.buf.lb(j, :)      = row(getfield_default(li, 'lb',    nan(nu_, 1)), nu_);
            obj.buf.ub(j, :)      = row(getfield_default(li, 'ub',    nan(nu_, 1)), nu_);
            obj.buf.active(j)     = double(getfield_default(li, 'active', NaN));
            obj.buf.du(j)         = getfield_default(li, 'du', NaN);
            obj.buf.satClip(j)    = getfield_default(li, 'sat_clip', NaN);
            obj.buf.cmdChange(j)  = getfield_default(li, 'command_changed', NaN);
            obj.buf.dVdt(j)       = getfield_default(li, 'dVdt', NaN);
            obj.buf.rhs(j)        = getfield_default(li, 'rhs', NaN);
            obj.buf.insideGrid(j) = double(getfield_default(li, 'inside_grid', NaN));
            obj.buf.ok(j)         = double(getfield_default(li, 'ok', NaN));
        end

        function finalize(obj)
            % Trim every buffer array to the number of recorded steps.
            n = obj.count;
            fn = fieldnames(obj.buf);
            for i = 1:numel(fn)
                v = obj.buf.(fn{i});
                if size(v, 1) >= n
                    obj.buf.(fn{i}) = v(1:n, :);
                end
            end
        end

        function tr = exportTrace(obj)
            % Emit a struct with the fields consumed by the LON tube-overlay
            % helpers (absolute coord mode) and the input-comparison plot.
            st = obj.buf.state;
            tr = struct();
            tr.uBody    = st(:, 4);
            tr.w        = st(:, 6);
            tr.q        = st(:, 11);
            tr.thetaDeg = rad2deg(st(:, 8));
            tr.k        = obj.buf.k;
            tr.brtExitFlag  = double(obj.buf.ok == 1 & obj.buf.brtV > 0);
            tr.gridExitFlag = double(obj.buf.insideGrid == 0);
            % Input-comparison fields.
            tr.u_nom     = obj.buf.u_nom;
            tr.u0        = obj.buf.u0;
            tr.u         = obj.buf.u;
            tr.lb        = obj.buf.lb;
            tr.ub        = obj.buf.ub;
            tr.active    = obj.buf.active;
            tr.du        = obj.buf.du;
            tr.satClip   = obj.buf.satClip;
            tr.cmdChange = obj.buf.cmdChange;
            tr.V         = obj.buf.brtV;
            tr.rhs       = obj.buf.rhs;
        end

        % ----------------------------------------------------------------
        % Plotting / saving
        % ----------------------------------------------------------------
        function plot(obj)
            % Draw all enabled figures and save them if configured.
            if ~obj.cfg.enable, return; end
            if obj.cfg.plotBasic
                figs = obj.plotBasic();
                obj.saveIfNeeded(figs);
            end
            if obj.cfg.plotFilter && obj.cfg.logFilter
                figs = obj.plotFilter();
                obj.saveIfNeeded(figs);
            end
        end

        function figs = plotBasic(obj)
            % Position / velocity / attitude / effector figures (ported from
            % the original run_transition_sim.m inline plots; data from buf).
            t   = obj.buf.t;
            st  = obj.buf.state;          % N x 12
            rp  = obj.buf.ref_pos;        % N x 3
            rv  = obj.buf.ref_vel;        % N x 3
            figs = struct();

            % --- Inertial position (NED) ---
            figs.position = figure('Name', 'Position (NED)');
            lbl = {'North [ft]', 'East [ft]', 'Down [ft]'};
            for i = 1:3
                subplot(3, 1, i);
                plot(t, st(:, i), 'b', t, rp(:, i), 'r--');
                ylabel(lbl{i}); grid on;
                if i == 1, legend('sim', 'ref'); title('Inertial position'); end
            end
            xlabel('Time [s]');

            % --- Body velocity ---
            figs.velocity = figure('Name', 'Body velocity');
            lbl = {'u [ft/s]', 'v [ft/s]', 'w [ft/s]'};
            for i = 1:3
                subplot(3, 1, i);
                plot(t, st(:, 3 + i), 'b', t, rv(:, i), 'r--');
                ylabel(lbl{i}); grid on;
                if i == 1, legend('sim', 'ref (heading frame)'); title('Velocity'); end
            end
            xlabel('Time [s]');

            % --- Attitude ---
            figs.attitude = figure('Name', 'Attitude');
            lbl = {'\phi [deg]', '\theta [deg]', '\psi [deg]'};
            for i = 1:3
                subplot(3, 1, i);
                plot(t, rad2deg(st(:, 6 + i)), 'b');
                ylabel(lbl{i}); grid on;
                if i == 1, title('Euler angles'); end
            end
            xlabel('Time [s]');

            % --- Effectors ---
            figs.effectors = figure('Name', 'Effectors');
            subplot(3, 1, 1);
            plot(t, obj.buf.engine(:, 1:8) .* (60 / (2 * pi)));
            ylabel('Lift rotors [RPM]'); grid on; title('Effectors');
            subplot(3, 1, 2);
            plot(t, obj.buf.engine(:, 9) .* (60 / (2 * pi)));
            ylabel('Pusher [RPM]'); grid on;
            subplot(3, 1, 3);
            plot(t, rad2deg(obj.buf.surface));
            ylabel('Surfaces [deg]'); grid on;
            legend('LA', 'RA', 'LE', 'RE', 'RUD');
            xlabel('Time [s]');
        end

        function figs = plotFilter(obj)
            % Single-run filter/BRT diagnostics:
            %   A) BRT value V(t)
            %   B) per-channel nominal (u_nom) vs filtered (u) input, with lb/ub.
            % Returns an empty struct (no figures) when the run carries no BRT
            % data (e.g. filter OFF -> brtV all NaN).
            figs = struct();
            k    = obj.buf.k;
            brtV = obj.buf.brtV;
            if all(isnan(brtV))
                return;
            end

            % --- A) BRT value over time ---
            figs.brt_value = figure('Name', 'BRT value', 'Color', 'w');
            plot(k, brtV, 'LineWidth', 1.2); hold on;
            yline(0, '--', 'V=0');
            xlabel('step k'); ylabel('V');
            title('Selected-anchor BRT value V(k)'); grid on;

            % --- B) per-channel nominal vs filtered input ---
            inputNames = {'Pi1','Pi2','Pi3','Pi4','Pi5','Pi6','Pi7','Pi8','Pi9','de','df'};
            unitNames  = {'RPM','RPM','RPM','RPM','RPM','RPM','RPM','RPM','RPM','deg','deg'};
            % Surface channels (elevator/flap) are stored in rad; show in deg.
            scale = ones(1, obj.nu);
            scale(10:11) = 180 / pi;

            u_nom = obj.buf.u_nom .* scale;
            u_f   = obj.buf.u     .* scale;
            lb    = obj.buf.lb    .* scale;
            ub    = obj.buf.ub    .* scale;

            figs.input_compare = figure('Name', 'Nominal vs filtered input', 'Color', 'w');
            tiledlayout(4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
            for i = 1:obj.nu
                nexttile;
                plot(k, u_nom(:, i), '--', 'LineWidth', 0.9); hold on;
                plot(k, u_f(:, i),   '-',  'LineWidth', 1.2);
                plot(k, lb(:, i), 'k:', 'LineWidth', 0.7);
                plot(k, ub(:, i), 'k:', 'LineWidth', 0.7);
                xlabel('k'); ylabel(sprintf('%s [%s]', inputNames{i}, unitNames{i}));
                title(inputNames{i}); grid on;
                if i == 1
                    legend({'u_{nom}', 'u filtered', 'lb', 'ub'}, 'Location', 'best');
                end
            end
        end

        function d = resolveSaveDir(obj)
            % Resolve (and create, once) the output directory. Cached so all
            % figures of one run land in the same folder.
            if ~isempty(obj.saveDirCached)
                d = obj.saveDirCached;
                return;
            end
            if ~isempty(obj.cfg.saveDir)
                d = obj.cfg.saveDir;
            else
                d = fullfile(pwd, 'logs', datestr(now, 'yyyymmdd_HHMMSS'));
            end
            if ~isfolder(d), mkdir(d); end
            obj.saveDirCached = d;
        end

        function saveFigure(obj, figHandle, name)
            % Save one figure as PNG into the resolved output directory.
            if ~obj.cfg.saveFigures, return; end
            d = obj.resolveSaveDir();
            exportgraphics(figHandle, fullfile(d, [name '.png']));
        end

        function saveIfNeeded(obj, figs)
            % Save every figure in a name->handle struct (field name = file name).
            if ~obj.cfg.saveFigures, return; end
            fn = fieldnames(figs);
            for i = 1:numel(fn)
                h = figs.(fn{i});
                if isgraphics(h)
                    obj.saveFigure(h, fn{i});
                end
            end
        end
    end
end

% -------------------------------------------------------------------------
function r = row(v, nu)
% Coerce a value to a 1 x nu row, padding/truncating defensively.
v = v(:)';
if numel(v) == nu
    r = v;
else
    r = nan(1, nu);
    m = min(numel(v), nu);
    r(1:m) = v(1:m);
end
end
