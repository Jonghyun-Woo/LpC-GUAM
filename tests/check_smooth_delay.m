% check_smooth_delay - does rounding the speed schedule delay the mission?
%
% No vehicle in the loop: this compares the reference schedules only, which
% is where a delay would have to come from.
clear;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt = 0.01;
params = struct('filter_mode','off');  params.T_seg = 2.0;
ref = ReferenceTrajectory.build('trim_schedule', dt, params);
t = ref.time;  v0 = ref.vel(1,:);

fprintf('\n%7s | %9s %9s %9s | %10s %10s\n', 'tau [s]', ...
        'start','reach 95%','reach end', 'distance','peak du/dt');
fprintf('%s\n', repmat('-',1,64));
for tau = [0 1 2 3 4]
    if tau == 0
        v = v0;
    else
        n = max(3, 2*floor(tau/dt/2)+1);
        pad = [repmat(v0(1),1,(n-1)/2), v0, repmat(v0(end),1,(n-1)/2)];
        v = conv(pad, ones(1,n)/n, 'valid');
    end
    vf = v0(end);
    i_s  = find(v > 0.01*vf, 1);
    i_95 = find(v >= 0.95*vf, 1);
    i_e  = find(v >= 0.999*vf, 1);
    fprintf('%7.1f | %9.2f %9.2f %9.2f | %10.1f %10.3f\n', tau, ...
            t(i_s), t(i_95), t(i_e), trapz(t, v), max(diff(v))/dt);
end
fprintf(['\n"start" is when the command first leaves zero, "reach end" when it\n' ...
         'first gets within 0.1%% of the final trim speed. "distance" is the\n' ...
         'integral of the whole schedule.\n']);
