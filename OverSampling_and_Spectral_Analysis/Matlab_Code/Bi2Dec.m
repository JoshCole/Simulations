function [ Dec ] = Bi2Dec(data, nBits)
% Bi2Dec
%
% This function calculates the binary to decimal values for a given packet
% Usage :
%
% [ Dec ] = Bi2Dec(data, nBits)
%
% Where         data        = Data to be converted to Decimal
%
%				nBits       = Number of Bits for symbol definition

[m, n, numberOfDec]=size(data);
Dec = zeros(numberOfDec,1);         % Define array size
powOf2 = fliplr(2.^(0:nBits-1));    % Calculates 2^N-1....2^2, 2^1, 2^0

for i = 1:nBits
    tempDec = powOf2(i)*data(1,i,:);    % Column wise multication for each bit
    Dec = squeeze(tempDec) + Dec;       % Decimal value
end

end

