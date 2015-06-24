function [ Bi ] = Dec2Bi(data, nBits)
lenData = length(data);
Bi = zeros(1,nBits,lenData);
powOf2 = fliplr(2.^(0:nBits-1));
for i = 1:lenData
    for ii = 1:nBits
        if data(i) >= powOf2(ii)
            data(i) = data(i)-powOf2(ii);
            Bi(1,ii,i) = 1;
        end
    end
end
end


