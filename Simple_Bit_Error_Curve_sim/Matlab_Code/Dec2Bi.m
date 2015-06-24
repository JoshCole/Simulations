function [ bin ] = Dec2Bi(data, nBits)
[m,n] = size(data);
Bi = false(nBits, n);
bin = cell(m,1);
for i = 1:m
    for ii = 1:n
        powOf2 = 2.^(0:nBits-1);
        cnt = nBits;
        cnt_2 = 1;
    
        while sum(data(i,1:n))>0
            if data(i,ii) >= powOf2(cnt)
                data(i,ii) = data(i,ii)-powOf2(cnt);
                Bi(cnt_2,ii) = true;
                
            end
            cnt = cnt - 1;
            cnt_2 = cnt_2+1;
        end
    end
    bin{i} = Bi;
end


