function [ Bi ] = Dec2Bi(data, nBits)
Bi = false(nBits,1);
powOf2 = 2.^(0:nBits-1);
cnt = nBits;
cnt_2 = 1;
while cnt>0
    if data >= powOf2(cnt)
        data = data-powOf2(cnt);
        Bi(cnt_2) = true;
    end
    cnt = cnt - 1;
    cnt_2 = cnt_2+1;
end
end


