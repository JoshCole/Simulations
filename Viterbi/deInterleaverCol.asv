function [ yi ] = deInterleaverCol( ti,ttiIndex)
TTI = ['10 ms' '20 ms' '40 ms' '80 ms'];
c1Array= [1 2 4 8];

C1 = c1Array(ttiIndex);
Xi = length(ti);
R1 = Xi/C1;
yi = zeros(R1,C1); % Interleaver

for i = 0:C1-1
    Index = i*R1;
    for ii = 1:R1
        yi(ii,i+1) = ti(Index+ii);
    end
end

end

