function C = brt_zero_contour(brtV, v, x_lon)
% BRT_ZERO_CONTOUR zero level of the interpolated V(.,v), sliced at the flown
% (w,q) offset, returned in absolute (u [ft/s], theta [deg]).
%
%   C = brt_zero_contour(brtV, v, x_lon)
%
% brtV  : BRTValue
% v     : command the tube belongs to [ft/s]
% x_lon : flown longitudinal state [u; w; q; theta]
% C     : cell array of 2xM contours, rows = [u; theta_deg]
%
% The tube is 4-D. Squashing it to (u, theta) is only meaningful if the other
% two axes are pinned to the (w, q) the vehicle actually passed through, which
% is what x_lon is for.
u = brtV.u_anchor;
if v <= u(1),        k = 1;               a = 0;
elseif v >= u(end),  k = numel(u) - 1;    a = 1;
else
    k = find(u <= v, 1, 'last');  k = min(k, numel(u)-1);
    a = (v - u(k)) / (u(k+1) - u(k));
end
xe = (1-a)*brtV.trim_lon(:,k) + a*brtV.trim_lon(:,k+1);
sk  = slice_uth(brtV, brtV.brt{k},   x_lon(2)-xe(2), x_lon(3)-xe(3));
sk1 = slice_uth(brtV, brtV.brt{k+1}, x_lon(2)-xe(2), x_lon(3)-xe(3));
s = (1-a)*sk + a*sk1;
cc = contourc(brtV.gv{1}, rad2deg(brtV.gv{4}), s', [0 0]);
C = {};  idx = 1;
while idx < size(cc,2)
    m = cc(2, idx);
    p = cc(:, idx+1 : idx+m);
    p(1,:) = p(1,:) + xe(1);
    p(2,:) = p(2,:) + rad2deg(xe(4));
    C{end+1} = p; %#ok<AGROW>
    idx = idx + m + 1;
end
end

function s = slice_uth(brtV, data, w_off, q_off)
gv = brtV.gv;
iw = interp1(gv{2}, 1:numel(gv{2}), min(max(w_off, gv{2}(1)), gv{2}(end)));
iq = interp1(gv{3}, 1:numel(gv{3}), min(max(q_off, gv{3}(1)), gv{3}(end)));
i0 = floor(iw); i1 = min(i0+1, numel(gv{2})); fw = iw - i0;
j0 = floor(iq); j1 = min(j0+1, numel(gv{3})); fq = iq - j0;
s = (1-fw)*(1-fq)*squeeze(data(:,i0,j0,:)) + fw*(1-fq)*squeeze(data(:,i1,j0,:)) ...
  + (1-fw)*   fq *squeeze(data(:,i0,j1,:)) + fw*   fq *squeeze(data(:,i1,j1,:));
end
