function [ Bi ] = grayCode2bi( grayCode, N, numBits )

Bi = false(numBits,N);

for ii = 1:N
    
    cnt = 1;
    
    while cnt<numBits
        Bi(1,ii) = grayCode(1,ii);
        a = Bi(cnt,ii);
        b = grayCode(cnt+1,ii);
        Bi(cnt+1,ii) = mod(a + b,2);
        cnt= cnt+1;
    end
end

end



