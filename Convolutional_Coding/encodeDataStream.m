function [ codeWord, genOutput, nBits, nextState ] = encodeDataStream( register, dataStream, xorIndex, numGenPoly, memory, method )

if strcmp(method,'Encode Data')
    processLength = length(dataStream);
elseif strcmp(method,'Trellis Map')
    processLength = 1;
end

genOutput = zeros(processLength,numGenPoly);
for i = 1:processLength
    if strcmp(method,'Encode Data')
        register(1) = dataStream(i);
    end
    for ii = 1:numGenPoly
        % Exicuite 1/3 decoding
        % ---------------------
        value1 = register(xorIndex{ii}(1));
        for j = 1: length(xorIndex{ii})-1
            value2 = register(xorIndex{ii}(j+1));
            value1 = xor(value1,value2);
        end
        genOutput(i,ii)= value1;
    end
    if strcmp(method,'Encode Data')
        for jj = length(register):-1:1
            if jj ~=1
                register(jj) = register(jj-1);
            end
        end
    end
end
    numElementsCode = length(genOutput(:,ii));
    codeWord = zeros(1,numGenPoly*numElementsCode);
for jj = 0:numGenPoly-1
    codeWord(jj+1:numGenPoly:end) = genOutput(:,jj+1);
end
numCodeWord = length(codeWord);
nBits = numCodeWord/numGenPoly;
nextState = register(1:memory);

end

