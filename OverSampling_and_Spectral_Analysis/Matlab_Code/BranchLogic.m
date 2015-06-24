function [ output_args ] = BranchLogic( input_args )
% Compute Branch Logic
% --------------------
inputValue = [0 1];
numInput = length(inputValue);
branchValues = zeros(states,numGenPoly,numInput );
nextState = zeros(states,mem,numInput );
for i = 1:states
    for ii = 1:numInput
        register = [inputValue(ii), statesBinary(i,:)];
        [ branchLogic, branchLogicBlock, nBits, State ] = encodeDataStream( register, newLLR, xorIndex, numGenPoly, mem, 'Trellis Map' );
        branchValues(i,:,ii) = branchLogic;
        nextState(i,:,ii) = State;
    end
end
[m,n]=size(newLLR);


end

