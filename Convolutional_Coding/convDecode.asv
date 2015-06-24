function [ output_args ] = convDecode( rxFrame, nG )
% Reshape frame to brinary stream
% -------------------------------
rxFrame = reshape(rxFrame,1,[]);
nBits = length(rxFrame);

% Define Trellis map for Decoder
% ------------------------------
stages  = nBits+nG;     % Number of stages in trellis
states = 2^nG;          % Number of states in trellis, dependent on encoding rate i.e 1/3

block_st= 1;

% Define struct current state
% ---------------------------
for i = 1:nG
    currentState.m{i} = 0;
end

currentState.st    = 1;
currentState.in    = 0;

% Define node for trellis map
% ---------------------------
for ii = 1:states
    node{ii}{ii}.p{1}       = NaN;
    node{ii}{ii}.p{2}       = NaN;
    node{ii}{ii}.f{1}       = NaN;
    node{ii}{ii}.f{2}       = NaN;
    node{ii}{ii}.cost       = -100000;
    node{ii}{ii}.visit      = 0;
    node{ii}{ii}.surviver   = NaN;
end
node{1}{1}.visit    = 1;
node{1}{1}.surviver = 1;
node{1}{1}.cost     = 0;
end