classdef LivenessPredictor < handle
    % LIVENESSPREDICTOR CBF-style finite-difference reachability forecast.
    %
    % Predictive gate (HJ-navigation paper Eq.(14) with the sign flipped for
    % liveness and alpha replaced by the schedule deadline):
    %
    %   reach OK   :  dV/dt <= -V / (t_rem - t_mar)
    %   governor   :  dV/dt >  -V / (t_rem - t_mar)
    %
    % evaluated in the equivalent "predicted value at the deadline" form,
    % which covers BOTH failure modes with one formula:
    %
    %   V_hat = V + (dV/dt) * t_eff        (t_eff = deadline-adjusted horizon)
    %   ready_pred = (V_hat < -eps_v)
    %
    %   - outside (V>=0), closing too slowly -> V_hat > 0 : LATE ENTRY
    %   - inside  (V<0),  V rising           -> V_hat > 0 : RE-EXIT
    %
    % dV/dt is a trailing finite difference of the MEASURED V along the
    % actual trajectory (paper Eq.(15)) - no linear model, so the WIKI 7.6
    % linearization problem does not apply.
    %
    % kappa_star answers the scalar-reference-governor question directly:
    % advancing the schedule at rate kappa stretches the time to the node to
    % t_rem/kappa, so the largest admissible advance rate is
    %
    %   kappa* = t_rem*dVdt / (dVdt*t_mar - V - eps_v)      (V >= 0, dVdt < 0)
    %
    % i.e. "how much must I slow the schedule down so that the forecast just
    % barely makes the deadline". kappa* = 1 means no intervention needed.
    %
    % reset() MUST be called whenever the anchor changes: V jumps
    % discontinuously there and the finite difference would be meaningless.

    properties
        fd_window = 0.2;    % [s] trailing window for the finite difference
        t_margin  = 0.3;    % [s] arrive-early margin t_mar
        eps_v     = 0.0;    % value margin at the deadline (grid-noise guard)

        buf_t = [];
        buf_V = [];
    end

    methods
        function obj = LivenessPredictor(fd_window, t_margin, eps_v)
            if nargin >= 1 && ~isempty(fd_window), obj.fd_window = fd_window; end
            if nargin >= 2 && ~isempty(t_margin),  obj.t_margin  = t_margin;  end
            if nargin >= 3 && ~isempty(eps_v),     obj.eps_v     = eps_v;     end
        end

        function reset(obj)
            obj.buf_t = [];
            obj.buf_V = [];
        end

        function out = update(obj, t, V, t_rem)
            % Push a sample and evaluate the predictive gate.
            %   t     : current time [s]
            %   V     : value of the NEXT-anchor BRT at the current state
            %   t_rem : time left until the schedule wants to hand over [s]
            obj.buf_t(end + 1) = t;
            obj.buf_V(end + 1) = V;
            keep = obj.buf_t >= t - obj.fd_window - 1e-9;
            obj.buf_t = obj.buf_t(keep);
            obj.buf_V = obj.buf_V(keep);

            % Trailing slope; during warm-up fall back to dVdt = 0, which
            % degenerates to the present-tense sign test.
            span = obj.buf_t(end) - obj.buf_t(1);
            warm = numel(obj.buf_t) >= 2 && span >= 0.5 * obj.fd_window;
            if warm
                dVdt = (obj.buf_V(end) - obj.buf_V(1)) / span;
            else
                dVdt = 0;
            end

            % The margin is applied in the conservative direction of each mode:
            %   entry  (V>=0): shorten the horizon -> must arrive t_mar early
            %   re-exit (V<0): extend  the horizon -> must still be inside
            %                  t_mar AFTER the deadline, otherwise an exit in
            %                  the last t_mar seconds goes unseen
            if V >= 0
                t_eff = max(t_rem - obj.t_margin, 0);
            else
                t_eff = max(t_rem + obj.t_margin, 0);
            end
            V_hat = V + dVdt * t_eff;

            % Largest schedule advance rate that still meets the deadline.
            if V < 0
                % Re-exit mode: stretching the deadline does not help (V is
                % rising), so there is no useful kappa - freeze and settle.
                kappa_star = 0;
            elseif dVdt >= -1e-9
                kappa_star = 0;                 % not closing at all
            else
                den = dVdt * obj.t_margin - V - obj.eps_v;
                kappa_star = t_rem * dVdt / den;
            end
            kappa_star = min(max(kappa_star, 0), 1);

            out = struct('dVdt', dVdt, 't_eff', t_eff, 'V_hat', V_hat, ...
                         'kappa_star', kappa_star, 'warm', warm, ...
                         'ready_pred', V_hat < -obj.eps_v, ...
                         'ready_now',  V < 0);
        end
    end
end
