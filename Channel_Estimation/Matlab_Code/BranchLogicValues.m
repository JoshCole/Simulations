function [ BranchNum, nextState ] = BranchLogicValues( numStates, inputValue, numInput, statesBinary, genPoly, NumberOfGenPoly, Memory, branchValues, nextState )
% BranchLogicValues
%
% Calculates for each state the Binary value of the Trellis branch for a
% Logic Zero and One, Then determines the row Number associated with a
% branch for example [1 1 1] is row number four in BranchLogic from
% function GetGenPolyData
%
% Usage :
%               [ BranchNum, nextState ] = BranchLogicValues( numStates,
%               inputValue, numInput, statesBinary, genPoly, NumberOfGenPoly, Memory, branchValues, nextState )
%
% Where         numStates = is number of States
%               inputValue = Logic [0,1] 
%               numInput = is 2 for Logic [0,1] Loop
%               statesBinary = State values in Binary
%               genPoly = The Encoding Polynomials
%               NumberOfGenPoly = Number of Polynomials
%               Memory = Memory for state register
%               branchValues = Empty Matrix
%               nextState = Empty Matrix

for i = 1:numStates
    for ii = 1:numInput
        register = [inputValue(ii), statesBinary(i,:)];         % State plus binary input i.e 1|001 where 001 is in the polynomial register
        [ branch, nState ] = EncodeBranch( register, genPoly, NumberOfGenPoly, Memory);
        branchValues(i,:,ii) = branch;      % Branch and Next state values
        nextState(i,:,ii) = nState;
    end
end
[ nextState ] = viterbiBi2Dec(nextState);       % Conversion to Decimal values for row numbering
[ BranchNum ] = viterbiBi2Dec(branchValues);

end

