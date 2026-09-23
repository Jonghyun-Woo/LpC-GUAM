function segs = brt_zero_contour(xv, yv, Zc)
% Extract V=0 BRT-boundary points from slice Zc(y,x) as a cell of 2xM segments
% (row 1 = x, row 2 = y). Returns {} when there is no zero crossing.
    segs = {};
    if all(Zc(:) > 0) || all(Zc(:) < 0)
        return;
    end

    C = contourc(double(xv), double(yv), double(Zc), [0 0]);
    k = 1;
    while k < size(C, 2)
        n = C(2, k);
        segs{end+1} = C(:, k+1:k+n); %#ok<AGROW>
        k = k + n + 1;
    end
end
