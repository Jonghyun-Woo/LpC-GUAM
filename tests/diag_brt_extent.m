% Check the actual BRT tube extent in each dimension from existing npy stacks.
% For each UH, finds max |dev| where V(x)<=0 exists anywhere in the other dims.
% Helps decide which grid bounds can be safely shrunk.
%
% Grid bounds below must match what the npy files were actually computed with.
% LON: brt_grid_config.m is correct.
% LAT: grid was already shrunk (HANDOFF.md); bounds differ from brt_grid_config.m.

clear; clc;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

stackDir = fullfile(root, 'reachable_data', 'GUAM_timestack_');
UH_LIST  = [1 5 10 15 20];
WH       = 3;

% Grid that the npy files were ACTUALLY computed with
% LON: from brt_grid_config.m (correct)
lon_gmin = [-35.0; -50.0; -1.8; -1.35];
lon_gmax = [ 35.0;  50.0;  1.8;  1.35];
% LAT: from the shrunk config (HANDOFF.md / old guam_analysis_config.yml)
lat_gmin = [-13.5; -3.4; -0.15; -1.3];
lat_gmax = [ 13.5;  3.4;  0.15;  1.3];

axes_ = struct( ...
    'name',   {'lon',  'lat'}, ...
    'gmin',   {lon_gmin, lat_gmin}, ...
    'gmax',   {lon_gmax, lat_gmax}, ...
    'labels', {{'u','w','q','theta'}, {'v','p','r','phi'}}, ...
    'units',  {{'ft/s','ft/s','rad/s','rad'}, {'ft/s','rad/s','rad/s','rad'}});

fprintf('BRT tube extent  (max |dev| projected over all other dims, max over UH {1,5,10,15,20})\n\n');

for ai = 1:2
    ax   = axes_(ai).name;
    gmin = axes_(ai).gmin;
    gmax = axes_(ai).gmax;
    half = (gmax - gmin) / 2;

    tube_max = zeros(4, numel(UH_LIST));

    for iu = 1:numel(UH_LIST)
        uh    = UH_LIST(iu);
        fname = fullfile(stackDir, sprintf('GUAM_%s_BRT_UH%d_WH%d.npy', upper(ax), uh, WH));
        V     = readNPY(fname);
        sz    = size(V);   % actual npy dimensions

        % build gv from actual npy size (not brt_grid_config gnum, which may be stale)
        gv = arrayfun(@(a,b,n) linspace(a,b,n), gmin, gmax, sz(:), 'UniformOutput', false);

        if iu == 1
            fprintf('== %s ==  npy size: [%s]\n', upper(ax), num2str(sz, '%d '));
        end

        for d = 1:4
            perm = [d, setdiff(1:4, d)];
            Vp   = permute(V, perm);
            Vp   = reshape(Vp, sz(d), []);
            Vp   = min(Vp, [], 2);
            inside = find(Vp <= 0);
            if isempty(inside)
                tube_max(d, iu) = 0;
            else
                tube_max(d, iu) = max(abs(gv{d}(inside)));
            end
        end
    end

    fprintf('  %-8s  grid_half  tube_max  tube%%   note\n', 'dim');
    for d = 1:4
        tm  = max(tube_max(d, :));
        gh  = half(d);
        lbl = axes_(ai).labels{d};
        unt = axes_(ai).units{d};
        pct = 100 * tm / gh;
        % flag if tube is >80% of grid (grid may be too tight already)
        flag = '';
        if pct > 80, flag = '  *** near edge'; end
        if pct < 50, flag = '  (can shrink)'; end
        fprintf('  %-8s  %8.3f   %7.3f  %5.1f%%  [%s]%s\n', lbl, gh, tm, pct, unt, flag);
    end
    fprintf('\n');
end
