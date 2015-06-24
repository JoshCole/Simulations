function [ yi, Xi, DeltaNij, Fi, Nij, rateChangeArray] = interLeaver1st( ti,ttiIndex, newPacketSize )
TTI = ['10 ms' '20 ms' '40 ms' '80 ms'];

c1Array= [1 2 4 8];

InterColumnPermutation{1} = 0;
InterColumnPermutation{2} = [0 1];
InterColumnPermutation{3} = [0 2 1 3];
InterColumnPermutation{4} = [0 4 2 6 1 5 3 7];

C1 = c1Array(ttiIndex);
Xi = length(ti);
R1 = Xi/C1;
yi = zeros(R1,C1); % Interleaver

% Create Interleaving Matrix
for i = 0:R1-1
    c1Index = i*C1;
    for ii = 1:C1
        yi(i+1,ii) = ti(c1Index+ii);
    end
end

% Use inter column permutation pattern
yi = yi(:,InterColumnPermutation{ttiIndex}+1);
[Nij,Fi] = size(yi);     % [Number of radio frames, number of bits per frame]

newRow = newPacketSize/C1;
if (newRow)<Nij
    DeltaNij = newRow - Nij;
elseif (newRow)>Nij 
    DeltaNij = newRow -Nij;
else
    DeltaNij = 0;
end
rateChangeArray = zeros(Nij+DeltaNij,1);
end

