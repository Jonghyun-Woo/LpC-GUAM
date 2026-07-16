% verify_sim_logger  -  SimLogger core buffer / exportTrace check.
%
% Runs a short closed-loop transition with the two-call logging contract and
% asserts buffer sizing, non-NaN state capture, and trace export shape.
% Prints "verify_sim_logger PASS" on success; errors otherwise.

here = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(here));

steps = 200;
cfg = Config('althold', struct('steps', steps));
g   = LpC_GUAM(cfg);
L   = SimLogger(cfg.logger, cfg.sim);

rt = g.refTraj;
N  = size(rt.pos, 2);
g.reset();
for k = 1:N
    ref.pos = rt.pos(:, k);  ref.vel = rt.vel(:, k);
    ref.chi = rt.chi(k);     ref.chi_dot = rt.chidot(k);

    L.logState(k, g, ref);
    [engine, surface] = g.step(ref);
    L.logInputs(g, engine, surface);
end
L.finalize();

assert(L.count == N, 'count %d ~= N %d', L.count, N);
assert(size(L.buf.state, 1) == N, 'state rows %d ~= N', size(L.buf.state, 1));
assert(~any(isnan(L.buf.state(:, 4))), 'uBody column has NaN');
assert(~any(isnan(L.buf.engine(:, 1))), 'engine column has NaN');

tr = L.exportTrace();
assert(numel(tr.uBody) == N && numel(tr.k) == N, 'trace length mismatch');
assert(all(isfinite(tr.uBody)) && all(isfinite(tr.thetaDeg)), 'trace has non-finite coords');

disp('verify_sim_logger PASS');
