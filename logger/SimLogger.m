classdef SimLogger < handle
    % Data logger for one closed-loop run (no plotting; that lives in visualize/).
    % Per step, call addState BEFORE guam.step() and addInputs AFTER it: the
    % filter diagnostics are only populated during step(). Units follow the
    % plant (ft, ft/s, rad, rad/s).

    properties
        cfg         % LoggerConfig
        dt          % [s] time step (from SimConfig)
        N           % buffer capacity (step count)
        count       % number of steps actually recorded
        nu          % filter effector count (LON = 11)
        buf         % struct of preallocated arrays (rows = steps)
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

            [alpha, beta, Vair] = guam.aeroFrame.compute(x(4), x(5), x(6));
            obj.buf.alpha(j) = alpha;
            obj.buf.beta(j)  = beta;
            obj.buf.V(j)     = Vair;

            obj.buf.ref_pos(j, :) = ref.pos(:)';
            obj.buf.ref_vel(j, :) = ref.vel(:)';
        end

        function addInputs(obj, controller, engine, surface)
            % Record post-step effectors and filter diagnostics into the record
            % opened by the matching addState(); must be called AFTER guam.step().
            % engine/surface are the actual actuator positions (guam.engine /
            % guam.surface); the filter diagnostics come from the controller.
            if ~obj.cfg.enable, return; end
            j = obj.count;
            if j < 1, return; end

            obj.buf.engine(j, :)  = engine(:)';
            obj.buf.surface(j, :) = surface(:)';

            if ~obj.cfg.logFilter, return; end
            filter_info = controller.safety_filter.last_info;
            if isempty(filter_info), return; end

            nu_ = obj.nu;
            obj.buf.brtV(j)       = getfield_default(filter_info, 'V', NaN);
            obj.buf.u_nom(j, :)   = row(getfield_default(filter_info, 'u_nom', nan(nu_, 1)), nu_);
            obj.buf.u0(j, :)      = row(getfield_default(filter_info, 'u0',    nan(nu_, 1)), nu_);
            obj.buf.u(j, :)       = row(getfield_default(filter_info, 'u',     nan(nu_, 1)), nu_);
            obj.buf.lb(j, :)      = row(getfield_default(filter_info, 'lb',    nan(nu_, 1)), nu_);
            obj.buf.ub(j, :)      = row(getfield_default(filter_info, 'ub',    nan(nu_, 1)), nu_);
            obj.buf.active(j)     = double(getfield_default(filter_info, 'active', NaN));
            obj.buf.du(j)         = getfield_default(filter_info, 'du', NaN);
            obj.buf.satClip(j)    = getfield_default(filter_info, 'sat_clip', NaN);
            obj.buf.cmdChange(j)  = getfield_default(filter_info, 'command_changed', NaN);
            obj.buf.dVdt(j)       = getfield_default(filter_info, 'dVdt', NaN);
            obj.buf.rhs(j)        = getfield_default(filter_info, 'rhs', NaN);
            obj.buf.insideGrid(j) = double(getfield_default(filter_info, 'inside_grid', NaN));
            obj.buf.ok(j)         = double(getfield_default(filter_info, 'ok', NaN));
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
            % Emit a struct with the fields consumed by the tube-overlay
            % helpers (absolute coord mode) and the input-comparison plot.
            % Both lon [u,w,q,theta] and lat [v,p,r,phi] states are provided.
            st = obj.buf.state;
            tr = struct();
            tr.uBody    = st(:, 4);
            tr.w        = st(:, 6);
            tr.q        = st(:, 11);
            tr.thetaDeg = rad2deg(st(:, 8));
            tr.v        = st(:, 5);
            tr.p        = st(:, 10);
            tr.r        = st(:, 12);
            tr.phiDeg   = rad2deg(st(:, 7));
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

        function data = exportData(obj)
            % Full logged buffer plus dt/nu, for the visualization layer
            % (visualize/plot_sim_diagnostics.m). Contains every field the basic
            % and filter diagnostic plots need.
            data    = obj.buf;
            data.dt = obj.dt;
            data.nu = obj.nu;
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
