function [ inputValue, numInput, PreviousState, Value, LogicValue, Decode ] = NodeSetup( Decision, numStates )
% NodeSetup
%
% Initalises Node Values
%
% Usage :
%               [ inputValue, numInput, currentStates, PreviousState,
%               Value, LogicValue, Decode, Surviver ] = NodeSetup( Decision, statesDec, numStates )
%
% Where         Decision = HardDecision (0,1) or SoftDecision (LLR)
%               numStates = Number of States

inputValue = [0 1];                         % Logic Values
numInput = length(inputValue);              % Length Logic Values           
PreviousState = zeros(numStates,1);
LogicValue = zeros(numStates,1);
Decode = 0;
switch Decision
    case 'HardDecision'
        Value(1:numStates,1) = 100000;      % Inital HardDecision Value
    case 'SoftDecision'
        Value(1:numStates,1) = -100000;     % Inital SoftDecision Value
end
end

