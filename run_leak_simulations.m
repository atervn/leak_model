%% Run leak model simulations
% Runs nSim independent simulations of the tight junction leak model and
% saves the results to results/leak_simulation_results.mat, which is read
% by analyze_leak_durations.m.
%
% The simulations run in parallel with the Parallel Computing Toolbox. 
% Without the toolbox, parfor runs as an ordinary for loop.
%
% Run this script from the folder that contains it.

% Export file name
exportFile = 'leak_simulation_results.mat';

%% Parameters
% Simulation
nSteps = 7200;           % number of time steps per simulation
nSim   = 100;            % number of independent simulations (2000 in the original runs)
seed   = 1;              % random seed, for reproducible results

% Network geometry
nStrands    = 4;         % number of horizontal strands
nComp       = 63;        % number of compartments in an odd compartment row
lStrandSame = 100e-9;    % strand length between side-by-side compartments in the same row (m)
lStrandDiff = 50e-9;     % strand length between compartments in neighbouring rows (m)

% Strand dynamics
pBreak = 0.013;          % breaking formation rate per µm of strand per time step
pSeal  = 0.012;          % sealing rate of a broken strand per time step

%% Build the network
% The network is the same for all simulations, so it is built only once
[nNodes, intRows, intCols, intLengths] = interaction_matrix(nStrands, nComp, lStrandSame, lStrandDiff);

%% Run the simulations
nLeaksInTime = cell(1,nSim);   % number of leaks at each time step
leakStarts   = cell(1,nSim);   % first time step of each leak
leakEnds     = cell(1,nSim);   % last time step of each leak (NaN = still open at the end)

% The random number stream is set inside the loop; keep the current one
% so it can be restored afterwards
originalStream = RandStream.getGlobalStream;

tic
parfor i = 1:nSim
    % Every simulation gets its own substream of the same random stream,
    % so the results are reproducible and do not depend on the number of
    % parallel workers
    stream = RandStream('mrg32k3a', 'Seed', seed);
    stream.Substream = i;
    RandStream.setGlobalStream(stream);

    [nLeaksInTime{i}, leakStarts{i}, leakEnds{i}] = leak_model(nSteps, nStrands, nNodes, nComp, ...
        intRows, intCols, intLengths, pBreak, pSeal, lStrandDiff);
end
fprintf('Ran %d simulations in %.1f s.\n', nSim, toc);

RandStream.setGlobalStream(originalStream);

%% Save the results
params = struct('nSteps', nSteps, 'nSim', nSim, 'seed', seed, ...
                'nStrands', nStrands, 'nComp', nComp, ...
                'lStrandSame', lStrandSame, 'lStrandDiff', lStrandDiff, ...
                'pBreak', pBreak, 'pSeal', pSeal);

if ~exist('results', 'dir')
    mkdir('results');
end
resultsFile = fullfile('results', exportFile);
save(resultsFile, 'params', 'nLeaksInTime', 'leakStarts', 'leakEnds');
fprintf('Results saved to %s.\n', resultsFile);
