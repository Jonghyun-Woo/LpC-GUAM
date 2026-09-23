function results = load_sim_results(path)
    % Load the results struct saved by run_transition_sim.
    if nargin < 1 || isempty(path)
        path = fullfile('reachable_data', 'transition_results.mat');
    end
    S = load(path);
    results = S.results;
end
