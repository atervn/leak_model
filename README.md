# Tight junction leak model

A stochastic model of leaks through the tight junction strand network. The
junction is modelled as a set of horizontal strands that divide the space
between them into rows of compartments. Strand segments break and seal at
random, and a leak exists when broken segments open a continuous path from
one side of the junction to the other. The model tracks the number of
distinct leaks over time and the lifetime of each leak.

## Files

| File | Description |
| --- | --- |
| `interaction_matrix.m` | Builds the strand network as a graph: compartments are nodes, strand segments between them are edges. |
| `leak_model.m` | Runs one simulation and tracks the leaks. Can optionally animate the simulation and record it as a video. |
| `run_leak_simulations.m` | Example: runs many independent simulations and saves the results to `results/`. |
| `analyze_leak_durations.m` | Example: analyzes and plots the durations of the leaks from the saved results. |
| `animate_leak_model.m` | Example: animates a single simulation of a small network. |

## Requirements

- MATLAB R2016b or later (uses the `graph` functions `distances` and `shortestpath`).
- Optional: Parallel Computing Toolbox, to run the simulations in parallel.
  Without it, `run_leak_simulations.m` runs the simulations one after another.

## Quick start

From the folder containing the files, run in MATLAB:

```matlab
animate_leak_model       % watch a single simulation
run_leak_simulations     % run the simulations (a few minutes with the default settings)
analyze_leak_durations   % analyze and plot the leak durations
```

The parameters are set at the top of each script. `run_leak_simulations.m`
uses `nSim = 100` by default so that it finishes quickly; the original runs
used `nSim = 2000`.

## Using the model directly

```matlab
[nNodes, intRows, intCols, intLengths] = interaction_matrix(nStrands, nComp, lStrandSame, lStrandDiff);
[nLeaks, leakStarts, leakEnds] = leak_model(nSteps, nStrands, nNodes, nComp, ...
    intRows, intCols, intLengths, pBreak, pSeal, lStrandDiff);
```

Add `'Plot', true` to animate the simulation, and `'VideoFile', 'leaks.mp4'`
to record it. See `help leak_model` and `help interaction_matrix` for details.

### Parameters

| Parameter | Description | Unit |
| --- | --- | --- |
| `nSteps` | Number of time steps | – |
| `nStrands` | Number of horizontal strands | – |
| `nComp` | Number of compartments in an odd compartment row (even rows have `nComp+1`) | – |
| `lStrandSame` | Strand length between side-by-side compartments in the same row | m |
| `lStrandDiff` | Strand length between compartments in neighbouring rows | m |
| `pBreak` | Breaking formation rate per µm of strand per time step | – |
| `pSeal` | Sealing rate of a broken strand per time step | – |

### Outputs

| Output | Description |
| --- | --- |
| `nLeaks` | Number of distinct leaks at each time step |
| `leakStarts` | First time step of each leak |
| `leakEnds` | Last time step of each leak (`NaN` if still open at the end of the simulation) |

## How leaks are defined

In each time step, the shortest open path is found for every pair of a
source node (one side of the junction) and a sink node (the other side)
that are connected through broken strands. Paths that share strand segments
belong to the same leak, which is represented by the segments common to its
paths (its bottleneck). A leak continues into the next time step if its
bottleneck shares a segment with a leak in that step. When two leaks merge,
the older one continues and the younger one ends.
