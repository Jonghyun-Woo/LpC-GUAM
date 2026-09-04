classdef ClosedLoopModel < handle
    % CLOSEDLOOPMODEL linear one-step model of the WHOLE governed loop.
    %
    % The governor asks "if I hold this command, is the vehicle safe over the
    % next T seconds?" and currently answers by stepping the real plant 600
    % times (~0.33 s per question). Everything inside that loop - RSLQR, the
    % control allocation, the actuator lags, the aero polynomial, the rigid
    % body and the state feedback - is either already linear or is linearised
    % here, so the whole lot collapses into one matrix per trim point:
    %
    %     xi(k+1) - xi_e(v)  =  Phi(v) * ( xi(k) - xi_e(v) )
    %
    % Phi is measured, not assembled: each state is nudged by +/-h and the
    % simulator is evaluated ONE step, so actuator dynamics, the allocation
    % and the hard-coded outer position gain are all captured without having
    % to model them by hand.
    %
    % State vector (34):
    %    1: 3   e_pos      reference position minus plant position [ft]
    %    4: 6   u, v, w    body velocity [ft/s]
    %    7: 9   phi,th,psi euler angles [rad]
    %   10:12   p, q, r    body rates [rad/s]
    %   13:21   eng_pos    9 propeller speeds [rad/s]
    %   22:26   srf_pos    5 surface deflections [rad]
    %   27:29   ctrl_lon_i longitudinal integrator
    %   30:32   ctrl_lat_i lateral integrator
    %   33:34   phi_v, theta_v
    %
    % Position enters as the ERROR, not the absolute value: holding a command
    % makes the reference ramp, so the absolute north position has no
    % equilibrium, while the tracking error settles to a constant (the ~38 ft
    % that drives RSLQR's +0.1*e_pos speed bias).

    properties (Constant)
        N_XI  = 34;
        Z_REF = -80;            % altitude reference used by the governor [ft]
    end

    properties
        w_ref = 0       % vertical profile the model is built at [ft/s, + down]
                        %   0      : level, the shipped behaviour
                        %   11.667 : the WH3 descent the trim points and their
                        %            reachable sets were computed for
        dt
        u_anchor        % 1xN trim speeds the model is built at [ft/s]
        Phi             % 1xN cell of N_XI x N_XI
        xi_e            % N_XI x N linearisation points
        c               % N_XI x N one-step drift, Psi(xi_e) - xi_e
        h               % N_XI x 1 perturbation sizes actually used
        resid           % 1xN  ||c||, how close xi_e is to a true equilibrium

        % xi_e is not required to be an exact fixed point. The closed loop has
        % a position-error mode with a ~10 s time constant (RSLQR's +0.1*e_pos
        % outer gain), so settling it exactly would take minutes. Storing the
        % leftover drift c and predicting with the affine form
        %     d(k+1) = c + Phi*d(k),      d = xi - xi_e
        % is exact to first order wherever xi_e happens to sit, and removes
        % the need to settle at all.
    end

    methods
        function obj = ClosedLoopModel(dt)
            obj.dt = dt;
        end

        % ---------------------------------------------------------------
        % state <-> vector
        % ---------------------------------------------------------------
        function xi = pack(obj, guam, ref_pos) %#ok<INUSL>
            s  = guam.saveState();
            xi = [ref_pos(:) - s.state(1:3); ...
                  s.state(4:12); ...
                  s.eng_pos(:); ...
                  s.srf_pos(:); ...
                  s.ctrl_lon_i(:); ...
                  s.ctrl_lat_i(:); ...
                  s.phi_v; s.theta_v];
        end

        function unpack(obj, guam, xi, ref_pos) %#ok<INUSL>
            s = guam.saveState();                      % keeps time etc.
            s.state(1:3)  = ref_pos(:) - xi(1:3);
            s.state(4:12) = xi(4:12);
            s.eng_pos     = xi(13:21);
            s.srf_pos     = xi(22:26);
            s.ctrl_lon_i  = xi(27:29);
            s.ctrl_lat_i  = xi(30:32);
            s.phi_v       = xi(33);
            s.theta_v     = xi(34);
            guam.restoreState(s);
        end

        % ---------------------------------------------------------------
        % the one-step map whose Jacobian we want
        % ---------------------------------------------------------------
        function xi2 = step_map(obj, guam, xi, v)
            % One closed-loop step with the command held at v. The absolute
            % reference position is arbitrary because only the error is a
            % state, so anchor it at the origin and advance it by v*dt.
            %
            % w_ref sets the vertical profile the model is linearised about.
            % It matters: the 20 trim points are WH3 points, trimmed at an
            % 11.667 ft/s descent, and settling the loop against a LEVEL
            % reference leaves the aircraft a standing -11.667 ft/s off every
            % anchor in heading-frame w and +2 to +4 deg in pitch. Measured
            % against the simulator with w_ref = 11.667 those deviations fall
            % to 0.02 ft/s and 0.01 deg. Phi and xi_e describe whichever
            % profile is set here, so it must match the mission being flown.
            P0 = [0; 0; obj.Z_REF];
            P1 = P0 + [v*obj.dt; 0; obj.w_ref*obj.dt];
            obj.unpack(guam, xi, P0);
            guam.step(struct('pos', P1, 'vel', [v; 0; obj.w_ref], ...
                             'chi', 0, 'chi_dot', 0));
            xi2 = obj.pack(guam, P1);
        end

        % ---------------------------------------------------------------
        % equilibrium for a held command
        % ---------------------------------------------------------------
        function [xi_eq, res] = find_equilibrium(obj, guam, v, T_settle, xi_0)
            % Hold the command and let the loop settle. Starting from a
            % neighbouring trim (xi_0) converges far faster than from hover.
            N = round(T_settle / obj.dt);
            if nargin >= 5 && ~isempty(xi_0)
                xi = xi_0;
            else
                xi = obj.pack(guam, [0; 0; obj.Z_REF]);
            end
            for i = 1:N
                xi = obj.step_map(guam, xi, v);
            end
            xi_eq = xi;
            res   = norm(obj.step_map(guam, xi_eq, v) - xi_eq);
        end

        % ---------------------------------------------------------------
        % Jacobian by central difference on the one-step map
        % ---------------------------------------------------------------
        function P = jacobian(obj, guam, xi_eq, v, hvec)
            n = obj.N_XI;
            P = zeros(n, n);
            for j = 1:n
                xp = xi_eq;  xp(j) = xp(j) + hvec(j);
                xm = xi_eq;  xm(j) = xm(j) - hvec(j);
                fp = obj.step_map(guam, xp, v);
                fm = obj.step_map(guam, xm, v);
                P(:, j) = (fp - fm) / (2*hvec(j));
            end
        end

        function [hv, tab] = tune_h(obj, guam, xi_eq, v, h0, n_halve)
            % Halve each perturbation until the column stops changing, which
            % is the balance Stevens 3.7 describes: too large mixes in
            % nonlinearity, too small drowns in round-off.
            n = obj.N_XI;  hv = h0(:);  tab = zeros(n, n_halve);
            for j = 1:n
                prev = [];  best = h0(j);
                for a = 1:n_halve
                    hj = h0(j) / 2^(a-1);
                    xp = xi_eq;  xp(j) = xp(j) + hj;
                    xm = xi_eq;  xm(j) = xm(j) - hj;
                    col = (obj.step_map(guam, xp, v) - ...
                           obj.step_map(guam, xm, v)) / (2*hj);
                    if ~isempty(prev)
                        d = norm(col - prev) / max(norm(prev), 1e-12);
                        tab(j, a) = d;
                        if d < 1e-3, best = hj; break; end
                    end
                    prev = col;  best = hj;
                end
                hv(j) = best;
            end
        end

        % ---------------------------------------------------------------
        % build the whole family
        % ---------------------------------------------------------------
        function build(obj, guam, u_anchor, h0, T_settle, verbose)
            if nargin < 6, verbose = true; end
            n = numel(u_anchor);
            obj.u_anchor = u_anchor(:)';
            obj.Phi      = cell(1, n);
            obj.xi_e     = zeros(obj.N_XI, n);
            obj.c        = zeros(obj.N_XI, n);
            obj.resid    = zeros(1, n);
            obj.h        = h0(:);

            guam.reset();
            xi_prev = [];
            for k = 1:n
                v = u_anchor(k);
                [xe, res] = obj.find_equilibrium(guam, v, T_settle, xi_prev);
                obj.xi_e(:, k) = xe;
                obj.c(:, k)    = obj.step_map(guam, xe, v) - xe;
                obj.resid(k)   = res;
                obj.Phi{k}     = obj.jacobian(guam, xe, v, obj.h);
                xi_prev        = xe;
                if verbose
                    fprintf(['trim %2d  v = %6.1f ft/s | drift %.2e | ' ...
                             'max|eig| = %.5f\n'], k, v, res, ...
                            max(abs(eig(obj.Phi{k}))));
                end
            end
        end

        % ---------------------------------------------------------------
        % interpolate, exactly like BRTValue does for the tubes
        % ---------------------------------------------------------------
        function [P, xe, cc] = at(obj, v)
            u = obj.u_anchor;
            if v <= u(1)
                k = 1;                 a = 0;
            elseif v >= u(end)
                k = numel(u) - 1;      a = 1;
            else
                k = find(u <= v, 1, 'last');  k = min(k, numel(u)-1);
                a = (v - u(k)) / (u(k+1) - u(k));
            end
            P  = (1-a)*obj.Phi{k}    + a*obj.Phi{k+1};
            xe = (1-a)*obj.xi_e(:,k) + a*obj.xi_e(:,k+1);
            cc = (1-a)*obj.c(:,k)    + a*obj.c(:,k+1);
        end

        % ---------------------------------------------------------------
        % linear rollout: the thing that replaces TrimRefGovernor.predict
        % ---------------------------------------------------------------
        function X = rollout(obj, xi_now, v, N)
            % Returns 4 x (N+1) of [u; w; q; theta] along the linear
            % prediction, i.e. exactly what worst_violation() needs.
            [P, xe, cc] = obj.at(v);
            d = xi_now - xe;
            X = zeros(4, N+1);
            % within xi: u=4, w=6, q=11, theta=8
            sel = [4 6 11 8];
            for k = 1:N+1
                X(:, k) = xe(sel) + d(sel);
                d = cc + P * d;                 % affine, see note on c above
            end
        end
    end
end
