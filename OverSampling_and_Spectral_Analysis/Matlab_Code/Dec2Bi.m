function [ Bi ] = viterbiDec2Bi(data, nBits)
m =length(data);
Bi = zeros(m,nBits);
powOf2 = 2.^(0:nBits-1);
for i = 1:m
    bitIndex = 1;
    powIndex = nBits;
    while powIndex>0
        if data(i) >= powOf2(powIndex)
            data(i) = data(i)-powOf2(powIndex);
            Bi(i,bitIndex) = 1;
        end
        powIndex = powIndex - 1;
        bitIndex = bitIndex+1;
    end
end

