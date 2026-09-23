function [xv, yv, Zc] = brt_project(V, gv, xdim, ydim)
% Slice 4-D BRT value function V onto the (xdim, ydim) plane, fixing the other
% two states at their deviation-zero node. Zc is oriented as Zc(y,x) for
% contour/contourc(xv, yv, Zc, ...). gv: 1x4 cell of deviation grid vectors.
    nd  = 4;
    fix = setdiff(1:nd, [xdim, ydim]);

    idx = repmat({':'}, 1, nd);
    for d = fix
        [~, i0] = min(abs(gv{d}));
        idx{d}  = i0;
    end

    Vs = squeeze(V(idx{:}));
    xv = gv{xdim};
    yv = gv{ydim};
    if xdim < ydim
        Zc = Vs.';
    else
        Zc = Vs;
    end
end
