%% Analyze leak durations
% Analyzes the results saved by run_leak_simulations.m:
%   1) leaky periods: continuous periods with at least one leak through
%      the junction, from the number of leaks at each time step
%   2) individual leaks: lifetimes of the separately tracked leaks
% Durations are in time steps.
%
% Run this script from the folder that contains it.

%% Settings
[fileName, folder] = uigetfile(fullfile('results', '*.mat'), 'Select a simulation results file');
if isequal(fileName, 0)
    disp('No results file selected.');
    return
end
resultsFile = fullfile(folder, fileName);
minDuration = 2;         % shorter periods/leaks (in time steps) are left out
binWidth    = 5;         % histogram bin width (time steps)
maxDuration = 200;       % histogram range; longer durations go into the last bin

%% Load the results
results = load(resultsFile);
params  = results.params;
nSim    = numel(results.nLeaksInTime);
fprintf('Loaded %d simulations of %d time steps.\n', nSim, params.nSteps);

% Histogram bins: the last bin also collects everything above maxDuration
edges   = [0:binWidth:maxDuration-binWidth, Inf];
centers = (0:binWidth:maxDuration-binWidth) + binWidth/2;

%% 1) Leaky periods
periodLengths = cell(1,nSim);  % durations of the leaky periods, per simulation
nPeriods      = zeros(1,nSim); % number of leaky periods, per simulation

for i = 1:nSim
    % Mark the steps with at least one leak, pad with a zero at both ends
    % and take the difference: +1 where a leaky period starts and -1 on
    % the step after it ends
    diffs = diff([0; results.nLeaksInTime{i}; 0] > 0);

    % Duration of each leaky period in time steps
    periodLengths{i} = find(diffs == -1) - find(diffs == 1);
    nPeriods(i) = numel(periodLengths{i});
end

% Durations from all simulations, leaving out the shortest ones
allPeriods = cell2mat(periodLengths(:));
allPeriods(allPeriods < minDuration) = [];

fprintf('\nLeaky periods (at least %d steps):\n', minDuration);
fprintf('  number:          %d (%.1f per simulation)\n', numel(allPeriods), numel(allPeriods)/nSim);
fprintf('  median duration: %g time steps\n', median(allPeriods));
fprintf('  fraction of time leaky: %.3f\n', ...
        mean(cellfun(@(n) mean(n > 0), results.nLeaksInTime)));

%% 2) Individual leaks
% Leaks that were already open at the first time step or still open at
% the end of the simulation have an unknown lifetime and are left out
allStarts = cell2mat(results.leakStarts(:));
allEnds   = cell2mat(results.leakEnds(:));
complete  = allStarts > 1 & ~isnan(allEnds);

leakLifetimes = allEnds(complete) - allStarts(complete) + 1;
leakLifetimes(leakLifetimes < minDuration) = [];

fprintf('\nIndividual leaks (at least %d steps, complete lifetime):\n', minDuration);
fprintf('  number:          %d (%.1f per simulation)\n', numel(leakLifetimes), numel(leakLifetimes)/nSim);
fprintf('  median lifetime: %g time steps\n', median(leakLifetimes));

%% Plot the duration distributions
figure('Name', 'Leak durations', 'NumberTitle', 'off');
hold on;
if ~isempty(allPeriods)
    counts = histcounts(allPeriods, edges);
    plot(centers, counts/sum(counts), '-o', 'DisplayName', 'Leaky periods');
end
if ~isempty(leakLifetimes)
    counts = histcounts(leakLifetimes, edges);
    plot(centers, counts/sum(counts), '-s', 'DisplayName', 'Individual leaks');
end
hold off;
xlabel('Duration (time steps)');
ylabel('Fraction');
title(sprintf('Leak durations, %d simulations', nSim));
legend('show');
box on;
