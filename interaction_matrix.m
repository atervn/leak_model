function [nNodes, intRows, intCols, intLengths] = interaction_matrix(nStrands,nComp,lStrandSame,lStrandDiff)
%INTERACTION_MATRIX Build the graph of the tight junction strand network.
%
%   [nNodes, intRows, intCols, intLengths] = INTERACTION_MATRIX(nStrands, nComp, lStrandSame, lStrandDiff)
%
%   The junction is modelled as nStrands horizontal strands. The space
%   between two neighbouring strands forms a row of compartments, so there
%   are nStrands-1 compartment rows. The rows alternate in size:
%     - odd rows contain nComp compartments,
%     - even rows contain nComp+1 compartments and are shifted by half a
%       compartment (half-width compartments at both ends), so that every
%       compartment touches two compartments in each neighbouring row.
%
%   Every compartment is a node in the graph and every piece of strand
%   that separates two compartments is an edge. Two extra sets of nodes
%   are added on the outside of the network:
%     - nComp source nodes (indices 1:nComp), one below each compartment
%       of the first row, and
%     - one sink node above each compartment of the last row.
%
%   Inputs:
%     nStrands    - number of horizontal strands
%     nComp       - number of compartments in an odd compartment row
%     lStrandSame - length of the strand segment separating two
%                   side-by-side compartments in the same row
%     lStrandDiff - length of the strand segment separating two diagonally
%                   adjacent compartments in neighbouring rows
%
%   Outputs:
%     nNodes      - total number of nodes (sources + compartments + sinks)
%     intRows     - first node of every edge
%     intCols     - second node of every edge (intRows < intCols)
%     intLengths  - strand length of every edge

% Total number of compartments: nComp per row, plus one extra compartment
% for every even row
nFunctions = (nStrands - 1)*nComp + floor((nStrands-1)/2);

% Creates a vector that will form the diagonal of the interaction matrix
% (interactions with self)
selfInteractions = zeros(nFunctions,1);

% Initialize the vector for the interactions and interaction lengths with
% the next compartments in the same row
nextInteractions = [];
nextInteractionLengths = [];

% Go through the strand rows
for i = 1:nStrands-1
    % Even strand rows
    if mod(i,2) == 0
        % Add the interactions for the row and zero in to end for the last
        % compartment
        nextInteractions = [nextInteractions ones(1,nComp) 0];
        % Add the interaction lenghts similarly to the interactions
        nextInteractionLengths = [nextInteractionLengths ones(1,nComp).*lStrandSame 0];
    % Odd strand rows
    else
        % Add the interactions for the row and zero in to end for the last
        % compartment
        nextInteractions = [nextInteractions ones(1,nComp-1) 0];
        % Add the interaction lenghts similarly to the interactions
        nextInteractionLengths = [nextInteractionLengths ones(1,nComp-1).*lStrandSame 0];
    end
end

% Remove the final extra zeros from the vectors
nextInteractions(end) = [];
nextInteractionLengths(end) = [];

% Create the interactionMatrix including the self and next interactions
% using the diagonal vectors
interactionMatrix = diag(selfInteractions) + diag(nextInteractions,1);
% Create the interactionLenghts matrix including the next interaction
% lengths using the diagonal vector (no length for self interaction)
interactionLengths = diag(nextInteractionLengths,1);

% If there are more than 2 strands
if nStrands > 2

    % Each compartment (excluding the top row) has two compartments above
    % them, initializing a vector for each and for the lengths
    aboveFirstInteractions = [];
    aboveSecondInteractions = [];
    aboveFirstInteractionLengths = [];
    aboveSecondInteractionLengths = [];

    % Go through the compartmentrows excluding the top one
    for i = 1:nStrands-2
        % Even row
        if mod(i,2) == 0

            % Left above interactions and lengths
            aboveFirstInteractions = [aboveFirstInteractions ones(1,nComp)];
            aboveFirstInteractionLengths = [aboveFirstInteractionLengths ones(1,nComp).*lStrandDiff];
            % Rigth above interactions and lengths
            aboveSecondInteractions = [aboveSecondInteractions ones(1,nComp) 0];
            aboveSecondInteractionLengths = [aboveSecondInteractionLengths ones(1,nComp).*lStrandDiff 0];

            % Odd row 
        else
            % Left above interactions and lengths
            aboveFirstInteractions = [aboveFirstInteractions ones(1,nComp) 0];
            aboveFirstInteractionLengths = [aboveFirstInteractionLengths ones(1,nComp).*lStrandDiff 0];
            % Rigth above interactions and lengths
            aboveSecondInteractions = [aboveSecondInteractions ones(1,nComp)];
            aboveSecondInteractionLengths = [aboveSecondInteractionLengths ones(1,nComp).*lStrandDiff];
        end
    end

    % If last compartment row is odd, remove the last zeros of the second
    % interactions and lengths
    if mod(nStrands-1,2) ~= 0
        aboveSecondInteractions(end) = [];
        aboveSecondInteractionLengths(end) = [];
    end

    % Add the above interactions to the interactionMatrix
    interactionMatrix = interactionMatrix + diag(aboveFirstInteractions,nComp) + diag(aboveSecondInteractions,nComp+1);
    % Add the above interaction lengths to the interactionLengths matrix
    interactionLengths = interactionLengths + diag(aboveFirstInteractionLengths,nComp) + diag(aboveSecondInteractionLengths,nComp+1);
end

% Add the source nodes in front of the compartments. Each iteration
% prepends one new node (a new first row and column). The 1 is placed at
% column nComp+1; because every later prepend shifts the indices by one,
% source node j ends up connected to compartment j of the first row.
for i = 1:nComp
    % New first row of the interaction matrix: connection to one
    % compartment of the first compartment row
    newRow = zeros(1,size(interactionMatrix,2)+1);
    newRow(nComp+1) = 1;
    % New first column is all zeros (the matrix is kept upper triangular,
    % so the connection is only stored in the row)
    newCol = zeros(size(interactionMatrix,2),1);
    interactionMatrix = [newRow; newCol interactionMatrix];

    % Same for the lengths: the strand between a source and a first-row
    % compartment spans the whole compartment width (2*lStrandDiff)
    newLRow = zeros(1,size(interactionLengths,2)+1);
    newLRow(nComp+1) = 2*lStrandDiff;
    interactionLengths = [newLRow; newCol interactionLengths];
end

% Number of compartments in the last (top) compartment row. With an even
% number of strands the last row (nStrands-1) is odd and has nComp
% compartments, otherwise it is even and has nComp+1 compartments.
if mod(nStrands,2) == 0
    nTopComp = nComp;
else
    nTopComp = nComp + 1;
end


% Add the sink nodes after the compartments, one for each top row
% compartment. Each iteration appends one new node (a new last column and
% row). As the matrix grows by one in each iteration, the 1 at row
% size-nTopComp+1 connects sink node i to compartment i of the top row.
for i = 1:nTopComp

    % New last column: connection from one top-row compartment to the sink
    newCol = zeros(size(interactionMatrix,1),1);
    newCol(size(interactionMatrix,1)-nTopComp+1) = 1;
    % New last row is all zeros (upper triangular matrix)
    newRow = zeros(1,size(interactionMatrix,2)+1);
    interactionMatrix = [interactionMatrix newCol; newRow];

    % Length of the strand between the top-row compartment and the sink
    newLCol = zeros(size(interactionLengths,1),1);
    % If the top row is a shifted (even) row, its first and last
    % compartments are only half as wide, so their strand is shorter
    if mod(nStrands,2) == 1 && (i == 1 || i == nTopComp)
        newLCol(size(interactionLengths,1)-nTopComp+1) = lStrandDiff;
    else
        % Full-width compartment
        newLCol(size(interactionLengths,1)-nTopComp+1) = 2*lStrandDiff;
    end
    interactionLengths = [interactionLengths newLCol; newRow];

end

% Collect all the connections (edges) from the finished matrix
idx = find(interactionMatrix);                 % 1‑based linear indices

% Total number of nodes: sources + compartments + sinks
nNodes = size(interactionMatrix,1);

% Convert those linear indices to row‑ and column‑subscripts
[intRows, intCols] = ind2sub(size(interactionMatrix), idx);

% Strand lengths of the edges, in the same order as intRows and intCols
intLengths = interactionLengths(idx);

end