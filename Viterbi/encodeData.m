function [ codeWord ] = encodeData( dataStream, genPoly, numGenPoly )
N = length(dataStream);
code = zeros(numGenPoly,N);
codeWord = zeros(1,numGenPoly*N);

for i = 1:numGenPoly
    temp = mod(conv(dataStream,genPoly(i,:)),2);
    code(i,:) = temp(1:N);
end


%Assemble code word from n outputs
for i = 0:numGenPoly-1
    codeWord(1+i:numGenPoly:end) = code(i+1,:);
end

end

