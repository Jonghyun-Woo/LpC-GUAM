classdef BRTSlice
    % Static helpers to slice a 4D LON BRT value function (u,w,q,theta) into
    % the (u,theta) plane at a given (w,q) perturbation, in ABSOLUTE
    % coordinates of the anchor trim. Shared by plot_trimref_brt and
    % plot_trimref_handover.
    %
    % Dimension-squish rule: always slice at the (w,q) the trajectory
    % actually passes through (perturbation wrt the anchor trim), so
    % inside/outside reads directly off the 2D plot.

    methods (Static)
        function s = uth(data, gv, w_pert, q_pert)
            % 2D V(u,theta) slice at fixed (w,q) perturbation, bilinear in (w,q).
            iw = interp1(gv{2}, 1:numel(gv{2}), ...
                         min(max(w_pert, gv{2}(1)), gv{2}(end)));
            iq = interp1(gv{3}, 1:numel(gv{3}), ...
                         min(max(q_pert, gv{3}(1)), gv{3}(end)));
            iw0 = floor(iw); iw1 = min(iw0 + 1, numel(gv{2})); fw = iw - iw0;
            iq0 = floor(iq); iq1 = min(iq0 + 1, numel(gv{3})); fq = iq - iq0;
            s = (1-fw)*(1-fq)*squeeze(data(:, iw0, iq0, :)) + ...
                   fw *(1-fq)*squeeze(data(:, iw1, iq0, :)) + ...
                (1-fw)*   fq *squeeze(data(:, iw0, iq1, :)) + ...
                   fw *   fq *squeeze(data(:, iw1, iq1, :));
        end

        function C = zero_contour(data, gv, lon_pt, trim_pt)
            % Zero-level contour of the (u,theta) slice of one BRT, taken at
            % the trajectory point's (w,q) perturbation wrt this BRT's anchor.
            % Returns cell of 2 x n polylines in ABSOLUTE (u [ft/s], theta [deg]).
            s  = BRTSlice.uth(data, gv, lon_pt(2) - trim_pt(2), ...
                                        lon_pt(3) - trim_pt(3));
            cc = contourc(gv{1}, rad2deg(gv{4}), s', [0 0]);
            C  = {};
            idx = 1;
            while idx < size(cc, 2)
                n = cc(2, idx);
                v = cc(:, idx + 1 : idx + n);
                v(1, :) = v(1, :) + trim_pt(1);
                v(2, :) = v(2, :) + rad2deg(trim_pt(4));
                C{end + 1} = v; %#ok<AGROW>
                idx = idx + n + 1;
            end
        end

        function draw(axp, data, gv, lon_pt, trim_pt, col)
            % Filled inside-region (V<0) + zero contour of one BRT slice,
            % absolute coords, plus the anchor trim marker.
            C = BRTSlice.zero_contour(data, gv, lon_pt, trim_pt);
            for c = 1:numel(C)
                v = C{c};
                patch(axp, 'XData', v(1, :), 'YData', v(2, :), ...
                      'FaceColor', col, 'FaceAlpha', 0.15, ...
                      'EdgeColor', col, 'LineWidth', 1.8);
            end
            plot(axp, trim_pt(1), rad2deg(trim_pt(4)), '^', 'Color', col, ...
                 'MarkerFaceColor', col, 'MarkerSize', 7);
        end
    end
end
