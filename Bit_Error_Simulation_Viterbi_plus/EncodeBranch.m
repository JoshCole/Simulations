function [ codeWord, nextState ] = EncodeBranch( register, genPoly, numGenPoly, Memory)
% EncodeBranch
%
% Encodes each state to determine branch metric values
%
% Usage :
%               [ codeWord, nextState ] = EncodeBranch( register, genPoly,
%               numGenPoly, Memory)
%
% Where         register = Input Logic|States Binary i.e 1|001
%               genPoly = Polynomials for encoding
%               Memory = Memory in register

xorIndex = cell(1,numGenPoly);              % Deterimes which elements of the register are Xored
genPolyElements = Memory+1;
for i = 1:numGenPoly
    arrayIndex = 1;
    for ii = 1:genPolyElements
        if (genPoly(i,ii))
            xorIndex{i}(arrayIndex) = ii;   % the Index at which a Xor is preformed in register
            arrayIndex = arrayIndex+1;
        end
    end
end

% Preforms Encoding of Branch
% ---------------------------
genOutput = zeros(1,numGenPoly);
for i = 1:numGenPoly
    Bit = register(xorIndex{i}(1));      % Loops through register preforming Xor for each value
    for ii = 1: length(xorIndex{i})-1
        Bit2 = register(xorIndex{i}(ii+1));
        Bit = xor(Bit,Bit2);
    end
    genOutput(1,i)= Bit;
end
numElementsCode = length(genOutput(:,i));
codeWord = zeros(1,numGenPoly*numElementsCode);

% Orders Code Word into Stream of Data
% ------------------------------------
for i = 0:numGenPoly-1
    codeWord(i+1:numGenPoly:end) = genOutput(:,i+1);
end
nextState = register(1:Memory);
end

