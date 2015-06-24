function [ strucStages, inputValue, numInput, currentStates ] = InitialiseNodes( numStates, statesBinary, statesDec, genPoly, NumberOfGenPoly, Memory )
% Define Nodes
% ------------
inputValue = [0 1];
numInput = length(inputValue);
currentStates = statesDec + 1;
[ BranchNum, NextState ] = BranchLogicValues( numStates, inputValue, numInput, statesBinary, genPoly, NumberOfGenPoly, Memory, zeros(numStates,NumberOfGenPoly,numInput ), zeros(numStates,Memory,numInput ) );

hammDistance = NaN;
branchValues = BranchNum;
nextState = NextState;

PreviousState = NaN(numStates,1);
Value(1:numStates,1) = -1000;
Logic = NaN(numStates,1);
Surviver = NaN;
Decode = NaN;
Visit = zeros(numStates,1);
% hammDistance = zeros(numStates,numInput,stages-1 );
strucStages = struct('States',{currentStates},'BranchLogic',{branchValues},'NextState',{nextState},'PreviousState',{PreviousState},'Visit',{Visit},'Value',{Value},'CodeWord',{[NaN NaN]},'SurviverPath',{Surviver},'Decode',{Logic},'SurviverDecode',{Decode},'HammDistance',{hammDistance});
% strucStages = repmat(strucStages,1,stages);


end

