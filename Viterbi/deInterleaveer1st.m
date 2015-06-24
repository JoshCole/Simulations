function [ codeWord, codeBlock] = deInterleaveer1st( ti,ttiIndex, numGenPoly )
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

% for i = 0:R1-1
%     Index = i*C1;
%     for ii = 1:C1
%         yi(i+1,ii) = ti(Index+ii);
%     end
% end
for i = 0:C1-1
    Index = i*R1;
    for ii = 1:R1
        yi(ii,i+1) = ti(Index+ii);
    end
end

% Use inter column permutation pattern
yi = yi(:,InterColumnPermutation{ttiIndex}+1);
[Nij,Fi] = size(yi);     % [Number of radio frames, number of bits per frame]
codeWord = zeros(1,Xi);
for i = 0:Fi-1
    codeWord(1+i:Fi:end) = yi(:,i+1);
end

rows = Xi/numGenPoly;
codeBlock = zeros(rows,numGenPoly);
for i = 0:rows-1
    Index = i*numGenPoly;
    for ii = 1:numGenPoly
        codeBlock(i+1,ii) = codeWord(Index+ii);
    end
end


end

