function [ di ] = readOutBits( yi, Xi, Fi )
di = zeros(1,Xi);
for i = 0:Fi-1
    di(1+i:Fi:end) = yi(:,i+1);
end

end

