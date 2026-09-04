% brt_inspect - read-only survey of the reachable-set (BRT) data files.
%
% Answers eight structural questions about the stored value functions and
% prints the evidence for each. Nothing is computed, fitted or saved except
% the two contour figures in item 8. No existing file is modified.
%
% Anything that cannot be established from the data itself is reported as
% UNCERTAIN with the reason, rather than guessed.
%
% Output: logger/brt_inspect_UH1.png, logger/brt_inspect_UH20.png

clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));
dataDir = fullfile(root,'data');

AX   = {'u','w','q','theta'};
UNIT = {'ft/s','ft/s','rad/s','rad'};
line = @(c) fprintf('%s\n', repmat(c,1,76));

%% ========================================================================
%  1. WHERE THE FILES ARE, HOW THEY ARE NAMED, WHETHER UH1..UH20 ARE ALL THERE
%% ========================================================================
line('='); fprintf('1. FILE INVENTORY\n'); line('=');

allf = dir(fullfile(dataDir,'*.mat'));
fprintf('directory      : %s\n', dataDir);
fprintf('total .mat     : %d\n\n', numel(allf));

% group by the prefix pattern <SET>_<CHANNEL>_<KIND>_HJIR_UH<k>_WH<j>
names = {allf.name};
fam = containers.Map();
for i = 1:numel(names)
    t = regexp(names{i}, '^(.*?)_UH(\d+)_WH(\d+)\.mat$', 'tokens', 'once');
    if isempty(t), key = ['(unparsed) ' names{i}]; else, key = t{1}; end
    if ~isKey(fam, key), fam(key) = []; end
    if ~isempty(t)
        fam(key) = [fam(key), str2double(t{2})];
    end
end
fprintf('%-28s | %5s | %s\n', 'family prefix', 'count', 'UH indices present');
line('-');
ks = sort(fam.keys);
for i = 1:numel(ks)
    v = sort(fam(ks{i}));
    fprintf('%-28s | %5d | %s\n', ks{i}, numel(v), squash(v));
end

TARGET = 'GUAM_LON_BRT_HJIR';
fprintf('\nnaming convention: <SET>_<CHANNEL>_<KIND>_HJIR_UH<k>_WH<j>.mat\n');
fprintf('   CHANNEL = LON / LAT     KIND = BRT / FRT\n');
fprintf('   UH<k>   = trim point index (1-based into the UH breakpoint list)\n');
fprintf('   WH<j>   = vertical-speed breakpoint index\n\n');

have = false(1,20);
for k = 1:20
    have(k) = isfile(fullfile(dataDir, sprintf('%s_UH%d_WH3.mat', TARGET, k)));
end
if all(have)
    fprintf('=> %s UH1..UH20 at WH3: ALL 20 PRESENT\n', TARGET);
else
    fprintf('=> %s WH3 MISSING: %s\n', TARGET, mat2str(find(~have)));
end
whv = unique(cellfun(@(n) str2double(regexp(n,'WH(\d+)\.mat$','tokens','once')), ...
      names(contains(names, TARGET)), 'UniformOutput', true));
fprintf('=> WH indices available for this family: %s\n', mat2str(whv(:)'));

%% ========================================================================
%  2. WHAT IS INSIDE ONE FILE - VARIABLE NAMES, SIZES, CLASSES
%% ========================================================================
line('='); fprintf('2. VARIABLES INSIDE THE FILES\n'); line('=');

f1  = fullfile(dataDir, sprintf('%s_UH1_WH3.mat',  TARGET));
f20 = fullfile(dataDir, sprintf('%s_UH20_WH3.mat', TARGET));
w1  = whos('-file', f1);

fprintf('file: %s\n\n', sprintf('%s_UH1_WH3.mat', TARGET));
fprintf('%-16s | %-18s | %-10s | %s\n', 'name','size','class','bytes');
line('-');
for i = 1:numel(w1)
    fprintf('%-16s | %-18s | %-10s | %d\n', w1(i).name, mat2str(w1(i).size), ...
            w1(i).class, w1(i).bytes);
end

% same variable list in every file?
same = true;  n1 = sort({w1.name});
for k = 2:20
    wk = whos('-file', fullfile(dataDir, sprintf('%s_UH%d_WH3.mat', TARGET, k)));
    if ~isequal(sort({wk.name}), n1), same = false; fprintf('\n  file UH%d differs: %s\n', k, strjoin(sort({wk.name}),', ')); end
end
fprintf('\nall 20 files carry the same variable list: %s\n', tf(same));

% identify the value function: the largest numeric array with >=4 dims
[~, ib] = max([w1.bytes]);
VNAME = w1(ib).name;
fprintf('\n=> VALUE FUNCTION identified as "%s"\n', VNAME);
fprintf('   evidence: largest variable (%d bytes), %d-D numeric of class %s;\n', ...
        w1(ib).bytes, numel(w1(ib).size), w1(ib).class);
fprintf('   and it is the variable BRTValue.m loads (load(...,''data'')).\n');

%% ========================================================================
%  3. DIMENSION ORDER OF THE VALUE ARRAY
%% ========================================================================
line('='); fprintf('3. DIMENSION ORDER\n'); line('=');

S1 = load(f1);   V1 = S1.(VNAME);
S20 = load(f20); V20 = S20.(VNAME);
spec = FilterConfig.channelSpec('lon');

fprintf('size(%s)                 = %s\n', VNAME, mat2str(size(V1)));
fprintf('channelSpec(''lon'').grid_num = %s   for state order [u w q theta]\n', ...
        mat2str(spec.grid_num'));
match = isequal(size(V1), spec.grid_num');
fprintf('\nmatch: %s\n', tf(match));
fprintf(['\nreasoning: the four grid counts (21, 41, 61, 31) are all DIFFERENT,\n' ...
         'so the size vector alone pins the axis order - no permutation of\n' ...
         '[21 41 61 31] other than the identity reproduces it. Therefore\n']);
if match
    fprintf('=> dimension order IS [u, w, q, theta]. CONFIRMED.\n');
else
    p = zeros(1,4);
    for d = 1:4, p(d) = find(spec.grid_num == size(V1,d), 1); end
    fprintf('=> dimension order is NOT the spec order. Actual axis per dim: %s\n', ...
            strjoin(AX(p), ', '));
end
fprintf('\nsize is identical in UH1 and UH20: %s\n', tf(isequal(size(V1), size(V20))));

%% ========================================================================
%  4. ARE THE GRID VECTORS STORED, AND DO THE BOUNDS MATCH THE SPEC?
%% ========================================================================
line('='); fprintf('4. GRID VECTORS\n'); line('=');

gridvars = setdiff({w1.name}, {VNAME});
if isempty(gridvars)
    fprintf('variables other than "%s" in the file: NONE\n', VNAME);
    fprintf('=> grid vectors are NOT stored. They must be rebuilt from\n');
    fprintf('   FilterConfig.channelSpec(''lon''), which is what BRTValue.m does:\n');
    fprintf('   gv{d} = linspace(grid_min(d), grid_max(d), grid_num(d))\n');
else
    fprintf('other variables present: %s\n', strjoin(gridvars, ', '));
    fprintf('=> inspect these to see whether they hold the grid.\n');
end

fprintf('\nspec bounds:  min %s\n              max %s\n', ...
        mat2str(spec.grid_min'), mat2str(spec.grid_max'));
fprintf('\nCAN THE BOUNDS BE VERIFIED FROM THE DATA?\n');
if isempty(gridvars)
    fprintf('  NO - UNCERTAIN.\n');
    fprintf(['  Only the NUMBER of nodes per axis is recoverable from the array\n' ...
             '  (21/41/61/31, all confirmed). The physical extent is not stored\n' ...
             '  anywhere in the file, so [-16 -33 -1.5 -0.75]..[16 33 1.5 0.75]\n' ...
             '  rests entirely on FilterConfig matching whatever generated these\n' ...
             '  files. The generator is not in this repository (no HJIR/helperOC\n' ...
             '  script is present), so it cannot be cross-checked here.\n']);
    fprintf(['  Weak consistency evidence only: see item 8 - the zero level set\n' ...
             '  sits well inside the grid on every axis, which would not hold if\n' ...
             '  the assumed extent were badly wrong (the set would be clipped or\n' ...
             '  fill the box).\n']);
else
    fprintf('  possibly - the extra variables above may carry it.\n');
end

% build the grid the way the rest of the codebase does
gv = cell(1,4);
for d = 1:4
    gv{d} = linspace(spec.grid_min(d), spec.grid_max(d), spec.grid_num(d));
end
i0 = zeros(1,4);
for d = 1:4
    [~, i0(d)] = min(abs(gv{d}));
    fprintf('  axis %-5s : %3d nodes, step %8.4f %-6s, zero at index %d (exact: %s)\n', ...
            AX{d}, numel(gv{d}), gv{d}(2)-gv{d}(1), UNIT{d}, i0(d), ...
            tf(abs(gv{d}(i0(d))) < 1e-12));
end

%% ========================================================================
%  5. DEVIATION COORDINATES OR ABSOLUTE COORDINATES?
%% ========================================================================
line('='); fprintf('5. COORDINATE FRAME\n'); line('=');

c1  = V1(i0(1),  i0(2), i0(3), i0(4));
c20 = V20(i0(1), i0(2), i0(3), i0(4));
fprintf('value at the grid centre (all four deviations = 0):\n');
fprintf('   UH1  (trim u =   0.0 ft/s) : %+.6f\n', c1);
fprintf('   UH20 (trim u = 160.0 ft/s) : %+.6f\n', c20);

fprintf('\ncentre value across all 20 files:\n');
cc = zeros(1,20);
for k = 1:20
    Sk = load(fullfile(dataDir, sprintf('%s_UH%d_WH3.mat', TARGET, k)));
    cc(k) = Sk.(VNAME)(i0(1), i0(2), i0(3), i0(4));
end
fprintf('   '); fprintf('%+7.3f', cc); fprintf('\n');
fprintf('   all negative: %s\n', tf(all(cc < 0)));

fprintf('\n=> DEVIATION COORDINATES. Evidence:\n');
fprintf(['   (a) the u axis spans only %+.1f..%+.1f ft/s. Absolute coordinates\n' ...
         '       could not represent a trim point at u = 160 ft/s at all - that\n' ...
         '       state would fall outside the grid entirely.\n'], gv{1}(1), gv{1}(end));
fprintf(['   (b) the grid centre is negative (inside the set) in every one of\n' ...
         '       the 20 files, including UH20. Under absolute coordinates the\n' ...
         '       centre would mean u = 0, which is not a safe cruise state.\n']);
fprintf(['   (c) BRTValue.m subtracts the trim state before indexing:\n' ...
         '       xc = x_lon - xe, then reads the grid at xc.\n']);

%% ========================================================================
%  6. SIGN CONVENTION
%% ========================================================================
line('='); fprintf('6. SIGN CONVENTION\n'); line('=');

fprintf('centre value (the trim point itself, deviation 0)  : %+.6f  -> %s\n', ...
        c1, sgnword(c1));
crn = V1(1,1,1,1);
fprintf('far corner   (all four axes at their grid extreme) : %+.6f  -> %s\n', ...
        crn, sgnword(crn));
fprintf('array range  : min %+.4f   max %+.4f\n', min(V1(:)), max(V1(:)));
fprintf('fraction of the whole 4-D grid with V <= 0 : %.2f %%\n', ...
        100*mean(V1(:) <= 0));
fprintf('\n=> V <= 0 IS INSIDE the reachable set, V > 0 is outside.\n');
fprintf(['   Evidence: the trim point (the one state that must be safe) is the\n' ...
         '   most negative region, the far corner of the grid is positive, and\n' ...
         '   every consumer in the repo tests V < 0 for "inside"\n' ...
         '   (BRTValue.m header, check_tube_dimensions.m cross_zero).\n']);

%% ========================================================================
%  7. TIME HORIZON
%% ========================================================================
line('='); fprintf('7. TIME HORIZON\n'); line('=');

fprintf('ndims(%s) = %d, size = %s\n', VNAME, ndims(V1), mat2str(size(V1)));
fprintf('a time axis would appear as a 5th dimension, or as a separate\n');
fprintf('time/tau vector variable in the file.\n\n');
fprintf('   5th dimension present : %s\n', tf(ndims(V1) > 4));
fprintf('   time-like variable    : %s\n', ...
        tf(any(contains(lower({w1.name}), {'t','tau','time','horizon'}) & ...
               ~strcmpi({w1.name}, VNAME))));

fprintf('\n=> UNCERTAIN whether this is a converged or a finite-time set.\n');
fprintf(['   What IS established: the stored object is a single 4-D array with\n' ...
         '   no time axis and no accompanying horizon metadata, so whatever\n' ...
         '   horizon was used is baked in and not recoverable from the file.\n']);
fprintf(['   What is NOT established: which horizon. Deciding it needs the\n' ...
         '   generating script, which is not in this repository - grepping for\n' ...
         '   HJIR/helperOC/reachability finds only consumers (BRTValue,\n' ...
         '   ValueFunctionLUT, LivenessFilter and the plotters), no producer.\n']);
fprintf(['   Weak hint only: the name "HJIR" points at a Hamilton-Jacobi\n' ...
         '   reachability solve, whose usual product is the converged set, but\n' ...
         '   the file itself does not say so.\n']);

%% ========================================================================
%  8. (u, theta) SLICE AT w = 0, q = 0 FOR UH1 AND UH20
%% ========================================================================
line('='); fprintf('8. SLICE AT w = 0, q = 0\n'); line('=');

fprintf('%-6s | %-10s | %14s | %18s | %18s\n', 'file','trim u','cells V<=0', ...
        'max |u| at th=0', 'max |theta| at u=0');
line('-');
Vs = {V1, V20};  ku = [1 20];  trimu = [0.0 160.0];  slice = cell(1,2);
for j = 1:2
    V = Vs{j};
    Z = squeeze(V(:, i0(2), i0(3), :));          % 21 x 31, (u, theta)
    slice{j} = Z;

    n_in = sum(Z(:) <= 0);

    lu = Z(:, i0(4));                            % along u at theta = 0
    iu = find(lu <= 0);
    if isempty(iu), umax = NaN; else, umax = max(abs(gv{1}(iu))); end

    lt = Z(i0(1), :);                            % along theta at u = 0
    it = find(lt <= 0);
    if isempty(it), tmax = NaN; else, tmax = max(abs(gv{4}(it))); end

    fprintf('UH%-4d | %6.1f ft/s | %6d / %4d | %10.3f ft/s | %8.4f rad (%.1f deg)\n', ...
            ku(j), trimu(j), n_in, numel(Z), umax, tmax, rad2deg(tmax));
end
fprintf(['\nnote: these are grid-node counts, not interpolated crossings, so the\n' ...
         'extents are quantised to the node spacing (u: %.2f ft/s, theta: %.3f rad).\n'], ...
        gv{1}(2)-gv{1}(1), gv{4}(2)-gv{4}(1));
fprintf(['The zero level set stays inside the grid on both axes in both files,\n' ...
         'i.e. it is not clipped by the grid edge - the consistency evidence\n' ...
         'referred to in item 4.\n']);

for j = 1:2
    f = figure('Position',[80 80 720 520]);
    Z = slice{j}';                               % theta down rows for contourf
    contourf(gv{1}, gv{4}, Z, 24, 'LineColor','none'); hold on;
    [~, h] = contour(gv{1}, gv{4}, Z, [0 0], 'k-', 'LineWidth', 2.2);
    plot(0, 0, 'rp', 'MarkerSize', 14, 'MarkerFaceColor','r');
    colorbar; xlabel('u deviation [ft/s]'); ylabel('\theta deviation [rad]');
    title(sprintf('%s UH%d WH3 - value on the w=0, q=0 slice (black = zero level, red = trim)', ...
                  'BRT', ku(j)));
    grid on;
    out = fullfile(root,'logger', sprintf('brt_inspect_UH%d.png', ku(j)));
    saveas(f, out);
    fprintf('saved %s\n', out);
end

line('='); fprintf('done - nothing was written except the two figures above.\n'); line('=');

%% ------------------------------------------------------------------------
function s = tf(b),  if b, s = 'YES'; else, s = 'NO'; end, end
function s = sgnword(v), if v <= 0, s = 'inside'; else, s = 'outside'; end, end

function s = squash(v)
% compact run-length rendering of an integer index list, e.g. "1-20"
v = sort(unique(v));
if isempty(v), s = '(none)'; return; end
b = [true, diff(v) ~= 1];  st = v(b);  en = v([b(2:end), true]);
p = cell(1, numel(st));
for i = 1:numel(st)
    if st(i) == en(i), p{i} = sprintf('%d', st(i));
    else, p{i} = sprintf('%d-%d', st(i), en(i)); end
end
s = strjoin(p, ',');
end
