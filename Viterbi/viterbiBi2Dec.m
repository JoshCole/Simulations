function [ Dec ] = viterbiBi2Dec(data)
[numData,nBits,numOfArrays] = size(data);
Dec = zeros(numData,numOfArrays);

for i = 1:numOfArrays
    bitIndex = 1;
    numBits = nBits;
    while numBits>0
        powOf2 = 2.^(numBits-1);
        numPow = powOf2*data(:,bitIndex,i);
        Dec(:,i) = Dec(:,i) + numPow;
        numBits = numBits -1;
        bitIndex = bitIndex + 1;
    end

end
Dec = Dec +1;
end

