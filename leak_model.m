function [nLeaks, leakStarts, leakEnds] = leak_model(nSteps, nStrands, nNodes, nComp, ...
        intRows, intCols, intLengths, pBreak, pSeal, lStrandDiff, varargin)
%LEAK_MODEL Simulate strand breaking and sealing and track the leaks.
%
%   [nLeaks, leakStarts, leakEnds] = LEAK_MODEL(nSteps, nStrands, nNodes,
%       nComp, intRows, intCols, intLengths, pBreak, pSeal, lStrandDiff)
%   simulates the tight junction strand network for nSteps time steps.
%   The network (nNodes, intRows, intCols, intLengths) is built with
%   INTERACTION_MATRIX using the same nStrands, nComp and lStrandDiff.
%
%   [...] = LEAK_MODEL(..., Name, Value) sets options:
%     'Plot'      - true to animate the simulation (default: false)
%     'PauseTime' - pause between animation frames in seconds (default: 0)
%     'VideoFile' - file to record the animation to, e.g. 'leaks.mp4' or
%                   'leaks.avi' (default: '' = no recording). Only used
%                   when 'Plot' is true.
%
%   Model:
%     Every strand segment (edge of the network) is either intact or
%     broken. In each time step an intact segment breaks with rate
%     pBreak per µm of its length, and a broken segment seals with
%     rate pSeal. A leak exists when broken segments form a
%     continuous path from a source node (one side of the junction) to a
%     sink node (the other side).
%
%     In each time step the shortest open path between every connected
%     source–sink pair is found. Paths that share edges are grouped into
%     one leak, represented by the edges common to its paths (the
%     bottleneck). A leak continues to the next time step if its
%     bottleneck shares an edge with a leak in that step. When two leaks
%     merge, the older one continues and the younger one ends.
%
%   Inputs:
%     nSteps      - number of time steps to simulate
%     nStrands    - number of horizontal strands
%     nNodes      - number of nodes (from INTERACTION_MATRIX)
%     nComp       - number of compartments in an odd compartment row
%     intRows     - first node of every edge (from INTERACTION_MATRIX)
%     intCols     - second node of every edge (from INTERACTION_MATRIX)
%     intLengths  - strand length of every edge in m (from INTERACTION_MATRIX)
%     pBreak      - breaking formation rate per µm of strand per time step
%     pSeal       - sealing rate of a broken strand per time step
%     lStrandDiff - strand length between compartments in neighbouring
%                   rows in m, used to set the initial state
%
%   Outputs:
%     nLeaks      - nSteps x 1, number of distinct leaks at each time step
%     leakStarts  - nLeakEvents x 1, first time step of each leak
%     leakEnds    - nLeakEvents x 1, last time step of each leak (NaN if
%                   the leak is still open at the end of the simulation)
%
%   See also INTERACTION_MATRIX.

%% -----------------------------------------------------------------------
%  Options.
% ------------------------------------------------------------------------
parser = inputParser;
parser.addParameter('Plot', false, @(x) islogical(x) || isnumeric(x));
parser.addParameter('PauseTime', 0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
parser.addParameter('VideoFile', '', @(x) ischar(x) || isstring(x));
parser.parse(varargin{:});
opts = parser.Results;

%% -----------------------------------------------------------------------
%  Initial state of the strands.
%  The probability that a strand of length lStrandDiff is broken is taken
%  from the steady state of its two-state (intact/broken) Markov chain:
%  pBreakDiff / (pBreakDiff + pSeal), where pBreakDiff is the per-step
%  breaking rate of a strand of length lStrandDiff. The chance that
%  a given strand starts out broken scales with its length.
% ------------------------------------------------------------------------
pBreakDiff = pBreak*lStrandDiff/1e-6;
pInitOpen  = pBreakDiff/(pBreakDiff + pSeal);

nInts = length(intRows);                     % total number of edges
breaks = rand(nInts,1) < pInitOpen .* intLengths ./ lStrandDiff;

%% -----------------------------------------------------------------------
%  Source and sink nodes.
%  INTERACTION_MATRIX places nComp source nodes first and one sink node per
%  top-row compartment last. The top row has nComp compartments when
%  nStrands is even and nComp+1 compartments when nStrands is odd.
% ------------------------------------------------------------------------
sources = 1:nComp;

if mod(nStrands,2) == 0
    % Even number of strands → nComp sinks
    sinks = (nNodes-nComp+1):nNodes;
else
    % Odd number of strands → nComp+1 sinks
    sinks = (nNodes-nComp):nNodes;
end

%% -----------------------------------------------------------------------
%  Output containers.
% ------------------------------------------------------------------------
nLeaks      = zeros(nSteps,1);   % number of distinct leaks at each step
leakStarts  = zeros(0,1);        % first time step of each leak
leakEnds    = zeros(0,1);        % last time step of each leak (NaN = still open)

pastLeakPathEdges = {};          % bottleneck edges of the leaks in the previous step
leakIds    = [];                 % IDs of those leaks (same order), index into leakStarts/leakEnds

%% -----------------------------------------------------------------------
%  Animation setup.
% ------------------------------------------------------------------------
if opts.Plot
    plotState = init_plot(nStrands, nComp, opts.VideoFile);
end

%% -----------------------------------------------------------------------
%  Main simulation loop – iterate over time steps.
% ------------------------------------------------------------------------
for t = 1:nSteps

    % ---------------------------------------------------------------
    %  1) Randomly decide which intact strands break and which broken
    %     strands seal, based on the supplied rates.
    % ---------------------------------------------------------------
    probs = rand(nInts,1);                     % uniform random numbers
    sealMask  = breaks & (probs < pSeal);                     % broken → intact
    breakMask = ~breaks & (probs < pBreak .* intLengths/1e-6);% intact → broken

    breaks(sealMask)  = false;   % apply sealing
    breaks(breakMask) = true;    % apply breaking

    % ---------------------------------------------------------------
    %  2) Build the graph of the currently broken (open) strands.
    % ---------------------------------------------------------------
    g = graph(intRows(breaks), intCols(breaks), [], nNodes);

    % ---------------------------------------------------------------
    %  3) Find the source–sink pairs that are connected by open strands.
    % ---------------------------------------------------------------
    distanceMatrix = distances(g, sources, sinks);
    reachable = isfinite(distanceMatrix);           % true where a path exists

    [srcIdx, snkIdx] = find(reachable);
    nRoutes = numel(srcIdx);                         % number of connected pairs

    LeakPathEdges = {};      % bottleneck edges of each distinct leak
    allOpenPaths  = {};      % edges of every shortest path (for plotting)

    % ---------------------------------------------------------------
    %  4) Find the shortest path of each connected pair and group the
    %     paths into distinct leaks.
    % ---------------------------------------------------------------
    for k = 1:nRoutes
        srcNode = sources(srcIdx(k));
        snkNode = sinks(snkIdx(k));

        % Edge indices of the shortest path
        [~,~,edges] = shortestpath(g, srcNode, snkNode);
        % The edge indices refer to g, which is rebuilt every step, so the
        % same strand gets a different index from step to step. Convert
        % them to stable edge keys (linear index of the node pair) so the
        % leaks can be compared between time steps.
        edgeNodes = g.Edges.EndNodes(edges,:);
        edges = sub2ind([nNodes nNodes], edgeNodes(:,1), edgeNodes(:,2))';
        allOpenPaths{end+1} = edges; %#ok<AGROW>

        % The first path starts the first leak. Each further path is
        % added to the first leak it shares edges with, which reduces
        % that leak to the shared edges; a path that shares no edges with
        % any leak starts a new leak.
        noMatch = true;
        for kk = 1:numel(LeakPathEdges)
            mask = ismember(LeakPathEdges{kk}, edges);
            if any(mask)
                LeakPathEdges{kk} = LeakPathEdges{kk}(mask);
                noMatch = false;
                break
            end
        end
        if noMatch
            LeakPathEdges{end+1} = edges; %#ok<AGROW>
        end
    end

    % ---------------------------------------------------------------
    %  5) Match the leaks to those of the previous time step to decide
    %     which leaks continue, which are new and which have ended.
    % ---------------------------------------------------------------
    % ID of each current leak, in the same order as LeakPathEdges
    % (0 = not matched to a past leak yet)
    currentIds = zeros(1,numel(LeakPathEdges));

    for k = 1:numel(pastLeakPathEdges)
        noMatch = true;
        for kk = 1:numel(LeakPathEdges)
            if any(ismember(pastLeakPathEdges{k}, LeakPathEdges{kk}))
                noMatch = false;
                if currentIds(kk) == 0
                    % Past leak still present → the current leak keeps its ID.
                    currentIds(kk) = leakIds(k);
                else
                    % Two past leaks merged into the same current leak:
                    % keep the older ID (smaller = started earlier) and
                    % end the younger leak at the previous step.
                    olderId   = min(currentIds(kk), leakIds(k));
                    youngerId = max(currentIds(kk), leakIds(k));
                    leakEnds(youngerId) = t-1;
                    currentIds(kk) = olderId;
                end
                break
            end
        end
        if noMatch
            % Past leak vanished – its last step was the previous one.
            leakEnds(leakIds(k)) = t-1;
        end
    end

    % Current leaks that were not matched to a past leak are new leaks.
    for kk = find(currentIds == 0)
        leakStarts(end+1,1) = t;                   % first step of this leak
        leakEnds(end+1,1)   = NaN;                 % still open
        currentIds(kk)      = numel(leakStarts);   % new unique ID
    end

    % Keep the current leaks and their IDs for the next time step.
    pastLeakPathEdges = LeakPathEdges;
    leakIds = currentIds;

    nLeaks(t) = numel(LeakPathEdges);

    % ---------------------------------------------------------------
    %  6) Update the animation.
    % ---------------------------------------------------------------
    if opts.Plot
        plotState = update_plot(plotState, t, nSteps, nNodes, intRows, intCols, ...
            breaks, allOpenPaths, LeakPathEdges, opts.PauseTime);
    end
end

if opts.Plot
    close_plot(plotState);
end
end


%% =======================================================================
%  Local functions for the animation
% ========================================================================

function plotState = init_plot(nStrands, nComp, videoFile)
% Create the figure, draw the nodes and open the video file (optional).

[plotState.nodeX, plotState.nodeY] = node_coordinates(nStrands, nComp);

plotState.fig = figure('Name','Tight junction leak model','NumberTitle','off');
hold on;
% Nodes (sources, compartments and sinks) as grey circles, these stay the
% same during the whole animation
plot(plotState.nodeX, plotState.nodeY, 'ko', 'MarkerSize', 4, ...
     'MarkerFaceColor', [0.9 0.9 0.9]);
axis equal tight off;

plotState.hOpenEdges       = gobjects(0);   % all open strands (grey)
plotState.hRouteEdges      = gobjects(0);   % shortest source–sink paths (blue)
plotState.hBottleneckEdges = gobjects(0);   % bottleneck edges of the leaks (red)

% Title for the statistics, aligned with the left edge of the axes
hTitle = title('', 'FontWeight', 'normal', 'HorizontalAlignment', 'left', ...
               'Units', 'normalized');
hTitle.Position(1) = 0;
plotState.hTitle = hTitle;

% Video recording
plotState.video = [];
if ~isempty(char(videoFile))
    [~,~,ext] = fileparts(char(videoFile));
    if strcmpi(ext, '.mp4')
        profile = 'MPEG-4';             % not available on all platforms
    else
        profile = 'Motion JPEG AVI';
    end
    plotState.video = VideoWriter(char(videoFile), profile);
    plotState.video.FrameRate = 15;     % frames per second
    plotState.video.Quality   = 95;     % compression quality, 0-100
    open(plotState.video);
end
end


function plotState = update_plot(plotState, t, nSteps, nNodes, intRows, intCols, ...
        breaks, allOpenPaths, LeakPathEdges, pauseTime)
% Redraw the open strands, routes and bottlenecks for time step t.

nodeX = plotState.nodeX;
nodeY = plotState.nodeY;

% Remove the edges drawn in the previous step
delete(plotState.hOpenEdges);
delete(plotState.hRouteEdges);
delete(plotState.hBottleneckEdges);
plotState.hRouteEdges = gobjects(0);
plotState.hBottleneckEdges = gobjects(0);

% All open strands as thin grey lines
iOpen = intRows(breaks);
jOpen = intCols(breaks);
plotState.hOpenEdges = draw_edges(nodeX, nodeY, iOpen, jOpen, [0.7 0.7 0.7], 1);

nRoutes = numel(allOpenPaths);
if nRoutes > 0
    % Edges of all shortest source–sink paths as blue lines (each edge
    % drawn once), decoded from the edge keys back to their end nodes
    [ri, rj] = ind2sub([nNodes nNodes], unique([allOpenPaths{:}]));
    plotState.hRouteEdges = draw_edges(nodeX, nodeY, ri, rj, [0.2 0.4 0.9], 2);

    % Bottleneck edges of each leak as thick red lines on top
    [bi, bj] = ind2sub([nNodes nNodes], unique([LeakPathEdges{:}]));
    plotState.hBottleneckEdges = draw_edges(nodeX, nodeY, bi, bj, 'r', 4);

    % Path length in edges, averaged over the connected pairs
    avgLen = mean(cellfun(@numel, allOpenPaths));
    txt = sprintf('Step %d/%d – %d open routes, avg len = %.1f, leaks = %d', ...
                  t, nSteps, nRoutes, avgLen, numel(LeakPathEdges));
else
    txt = sprintf('Step %d/%d – no open route across the junction', t, nSteps);
end
plotState.hTitle.String = txt;

drawnow limitrate;
pause(pauseTime);

% Add the frame to the video
if ~isempty(plotState.video)
    writeVideo(plotState.video, getframe(plotState.fig));
end
end


function close_plot(plotState)
% Finalize the video file.
if ~isempty(plotState.video)
    close(plotState.video);
end
end


function h = draw_edges(nodeX, nodeY, i, j, color, width)
% Draw edges between nodes i and j as lines, one line object per edge.
if isempty(i)
    h = gobjects(0);
    return
end
i = i(:)'; j = j(:)';
h = line([nodeX(i); nodeX(j)], [nodeY(i); nodeY(j)], 'Color', color, 'LineWidth', width);
end


function [nodeX, nodeY] = node_coordinates(nStrands, nComp)
% Plotting coordinates of all nodes, in the same order as the node
% indices from INTERACTION_MATRIX: sources, compartment rows 1 to
% nStrands-1, sinks. Each node row is drawn at its own height, sources at
% the bottom (y = 1).

% Sources (y = 1) and the first compartment row (y = 2), nComp each
nodeX = [1:nComp, 1:nComp];
nodeY = [ones(1,nComp), 2*ones(1,nComp)];

% Remaining compartment rows
for ii = 2:nStrands-1
    if mod(ii,2) == 0
        % Even row: nComp+1 compartments shifted by half a compartment,
        % the first and last ones are half-width
        nodeX = [nodeX, 0.75, 1.5:nComp, nComp+0.25]; %#ok<AGROW>
        nodeY = [nodeY, (1+ii)*ones(1,nComp+1)];      %#ok<AGROW>
    else
        % Odd row: nComp compartments aligned with the first row
        nodeX = [nodeX, 1:nComp];                     %#ok<AGROW>
        nodeY = [nodeY, (1+ii)*ones(1,nComp)];        %#ok<AGROW>
    end
end

% Sinks on top, laid out like the last compartment row
if mod(nStrands,2) == 0
    nodeX = [nodeX, 1:nComp];
    nodeY = [nodeY, (nStrands+1)*ones(1,nComp)];
else
    nodeX = [nodeX, 0.75, 1.5:nComp, nComp+0.25];
    nodeY = [nodeY, (nStrands+1)*ones(1,nComp+1)];
end
end
