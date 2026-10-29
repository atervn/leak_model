%% Animate the leak model
% Runs one simulation of a small network and shows it as an animation.
%   grey lines - broken (open) strands
%   blue lines - shortest open paths across the junction
%   red lines  - bottleneck edges of each distinct leak
%
% To save the animation as a video, set videoFile below.

%% Parameters
% Simulation
nSteps    = 800;         % number of time steps
pauseTime = 0;           % pause between frames in seconds (animation speed)
videoFile = '';          % e.g. 'leak_animation.mp4' ('' = no video)
seed      = 1;           % random seed, for a reproducible animation

% Network geometry (smaller than in run_leak_simulations.m so the
% animation stays readable)
nStrands    = 4;         % number of horizontal strands
nComp       = 20;        % number of compartments in an odd compartment row
lStrandSame = 100e-9;    % strand length between side-by-side compartments in the same row (m)
lStrandDiff = 50e-9;     % strand length between compartments in neighbouring rows (m)

% Strand dynamics (faster than in run_leak_simulations.m, so that leaks
% appear and disappear within a short animation)
pBreak = 0.08;           % breaking probability per µm of strand per time step
pSeal  = 0.033;          % sealing probability of a broken strand per time step

%% Build the network and run the animation
[nNodes, intRows, intCols, intLengths] = interaction_matrix(nStrands, nComp, lStrandSame, lStrandDiff);

rng(seed);
nLeaks = leak_model(nSteps, nStrands, nNodes, nComp, intRows, intCols, intLengths, ...
    pBreak, pSeal, lStrandDiff, 'Plot', true, 'PauseTime', pauseTime, 'VideoFile', videoFile);

fprintf('Fraction of time steps with at least one leak: %.2f\n', mean(nLeaks > 0));
