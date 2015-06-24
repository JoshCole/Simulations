function [ Dec ] = Bi2Dec(data, nBits)
Dec = 0;
cnt = nBits;
cnt_2 = 1;

while cnt>0
    powOf2 = 2.^(cnt-1);
    temp = powOf2*data(cnt_2,:);
    Dec = Dec + temp;
    cnt = cnt -1;
    cnt_2 = cnt_2 + 1;
end

end

