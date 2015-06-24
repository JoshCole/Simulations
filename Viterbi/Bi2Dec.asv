function [ Dec ] = Bi2Dec(data, nBits)
[m, n, numberOfDec]=size(data);
Dec = zeros(numberOfDec,1);
powOf2 = fliplr(2.^(0:nBits-1));

for i = 1:nBits
    tempDec = powOf2(i)*data(1,i,:);
    Dec = squeeze(tempDec) + Dec;
end

end

