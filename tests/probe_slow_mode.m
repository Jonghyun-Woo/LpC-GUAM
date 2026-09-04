% probe_slow_mode - is the near-unit mode the reason the prediction is wrong,
% and can it simply be left out?
%
% ---------------------------------------------------------------------------
% THE QUESTION
% ---------------------------------------------------------------------------
% Phi carries a mode at |lambda| ~ 0.999, right under the stability boundary.
% Doubling the perturbation used to measure Phi pushes it past 1 (0.99907 ->
% 1.00057 at the hover anchor), which flips a 43 s transition from decaying by
% 50x to growing by 11x. The finite difference cannot pin this mode down, and
% no perturbation size fixes that - see project-secant-fails.
%
% But the mode lives in the position error and the longitudinal integrator,
% and POSITION IS NOT IN THE REACHABLE SETS. The safety test only ever reads
% u, w, q and theta. So the question is not whether the mode is well modelled;
% it is whether it reaches the four states that matter.
%
% ---------------------------------------------------------------------------
% THE TEST
% ---------------------------------------------------------------------------
% Two things are measured, both against simulator trajectories already
% recorded by tests/measure_model_validity:
%
%   1. PARTICIPATION - how much of each slow eigenvector sits in the four
%      rows the sets read. A mode that lives entirely in position has nothing
%      to say about safety however badly it is modelled.
%
%   2. DOES DROPPING IT HELP - propagate the error with the slow modes removed
%      from Phi and compare the (u,w,q,theta) prediction against the simulator.
%      If the prediction gets BETTER, the mode was hurting and can go. If it
%      gets worse, the mode is carrying real signal and the 19 s stays.
%
% Output: logger/slow_mode.mat
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

SEL = [4 6 11 8];
LAM_SLOW = 0.995;          % above this counts as the marginal family
line = @(c) fprintf('%s\n', repmat(c, 1, 92));

Sm = load(fullfile(root,'logger','closed_loop_model_descending.mat'));
M  = Sm.M;   dXE = Sm.dXE;   dt = M.dt;
Sv = load(fullfile(root,'logger','model_validity.mat'));   R = Sv.R;

T   = load('trim_table_Poly_ConcatVer4p0.mat');
UH  = T.UH(1:20);  UH = UH(:)';
trim_lon = T.XU0_interp([1 3 5 11], 1:20, 3);

%% 1. what do the slow modes consist of? -------------------------------
line('='); fprintf('WHAT THE MARGINAL MODES ARE MADE OF  (|lambda| > %.3f)\n', LAM_SLOW); line('=');
grp = {1:3,'pos err'; 4:6,'u v w'; 7:9,'euler'; 10:12,'p q r'; ...
       13:21,'rotors'; 22:26,'surf'; 27:32,'integr'; 33:34,'virt att'};
fprintf('   k  u_trim  n_slow  max|lam|   share of the slow eigenvectors, by block\n');
fprintf('%34s', '');
for g = 1:size(grp,1), fprintf('%9s', grp{g,2}); end
fprintf('\n');
share = zeros(20, size(grp,1));   part4 = zeros(1,20);
for k = 1:20
    [Ve, D] = eig(M.Phi{k});   lam = diag(D);
    sl = abs(lam) > LAM_SLOW;
    W  = abs(Ve(:, sl));
    W  = W ./ max(sum(W.^2, 1), eps);          % energy-normalise each mode
    e  = sum(W.^2, 2);   e = e / sum(e);       % pooled energy per state
    for g = 1:size(grp,1), share(k,g) = sum(e(grp{g,1})); end
    part4(k) = sum(e(SEL));
    fprintf('%4d %7.1f %7d %9.5f', k, UH(k), sum(sl), max(abs(lam(sl))));
    fprintf('%34s', '');
    for g = 1:size(grp,1), fprintf('%8.1f%%', 100*share(k,g)); end
    fprintf('\n');
end
fprintf('\n  energy of the marginal modes sitting in [u w q theta] : %.1f%% .. %.1f%%\n', ...
        100*min(part4), 100*max(part4));

%% 2. does dropping them improve the prediction? -----------------------
line('='); fprintf('PREDICTION WITH THE MARGINAL MODES REMOVED\n'); line('=');
fprintf('  error against the simulator, max over the ramp, as a fraction of\n');
fprintf('  the true excursion. Lower is better.\n\n');

Phi_cut = cell(1,20);
for k = 1:20
    [Ve, D] = eig(M.Phi{k});   lam = diag(D);
    lam(abs(lam) > LAM_SLOW) = 0;               % the mode contributes nothing
    Phi_cut{k} = real(Ve * diag(lam) / Ve);
end

nm = {'u','w','q','th'};  sc = [1 1 180/pi 180/pi];
fprintf('    a  | channel | full Phi | slow removed | better?\n');
res = struct();
for ia = 1:numel(R)
    a = R(ia).a;   t = R(ia).t;   v = R(ia).v;
    ERs = R(ia).ER;                                   % simulator, heading frame
    ELf = propagate(M, dXE, M.Phi,  v, dt, SEL, trim_lon, UH);
    ELc = propagate(M, dXE, Phi_cut, v, dt, SEL, trim_lon, UH);
    m = t >= 45;
    for c = 1:4
        ef = max(abs(ELf(c,m) - ERs(c,m)));
        ec = max(abs(ELc(c,m) - ERs(c,m)));
        tag = ' ';
        if ec < ef*0.95, tag = 'YES'; elseif ec > ef*1.05, tag = 'no'; end
        fprintf(' %5.1f | %-7s | %8.2f | %12.2f | %s\n', ...
                a, nm{c}, sc(c)*ef, sc(c)*ec, tag);
    end
    res(ia).a = a;  res(ia).ELf = ELf;  res(ia).ELc = ELc;  res(ia).ER = ERs;
end

save(fullfile(root,'logger','slow_mode.mat'), 'res', 'share', 'part4', 'LAM_SLOW');
fprintf('\nsaved logger/slow_mode.mat\n');

%% ---------------------------------------------------------------------
function EL = propagate(M, dXE, Phi, v, dt, SEL, trim_lon, UH)
% Same linear propagation measure_model_validity uses, but with whichever Phi
% is handed in. Returned as the heading-frame deviation from the trim curve,
% which is what a set lookup sees.
N = numel(v);   XL = zeros(4, N);   d = zeros(34,1);
ua = M.u_anchor;
for i = 1:N
    vv = v(i);
    if vv <= ua(1),        k = 1;            al = 0;
    elseif vv >= ua(end),  k = numel(ua)-1;  al = 1;
    else
        k = find(ua <= vv, 1, 'last');  k = min(k, numel(ua)-1);
        al = (vv - ua(k))/(ua(k+1) - ua(k));
    end
    P  = (1-al)*Phi{k}      + al*Phi{k+1};
    xe = (1-al)*M.xi_e(:,k) + al*M.xi_e(:,k+1);
    XL(:, i) = xe(SEL) + d(SEL);
    adot = 0;
    if i < N, adot = (v(i+1)-v(i))/dt; end
    d = P*d - dXE(:, k)*adot*dt;
end
th = XL(4,:);
HL = [ XL(1,:).*cos(th) + XL(2,:).*sin(th) ;
      -XL(1,:).*sin(th) + XL(2,:).*cos(th) ;
       XL(3,:) ; th ];
XT = zeros(4, N);   vv = min(max(v, UH(1)), UH(end));
for c = 1:4, XT(c,:) = interp1(UH, trim_lon(c,:), vv); end
EL = HL - XT;
end
