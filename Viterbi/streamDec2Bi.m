function [ Bi ] = streamDec2Bi(data, nBits)
m =length(data);
Bi = zeros(1,m*nBits);
powOf2 = 2.^(0:nBits-1);
bitIndex = 1;
for i = 1:m
    powIndex = nBits;
    while powIndex>0
        if data(i) >= powOf2(powIndex)
            data(i) = data(i)-powOf2(powIndex);
            Bi(bitIndex) = 1;
        end
        powIndex = powIndex - 1;
        bitIndex = bitIndex+1;
    end
end